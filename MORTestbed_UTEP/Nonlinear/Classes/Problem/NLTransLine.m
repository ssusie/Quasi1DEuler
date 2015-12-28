classdef NLTransLine
    
    properties
        config;
        staticFlag = false;
        
        mesh;
        ndim;
        Jstruct;
        ic;
        B;
        G;
        C;
        D;
    end
    
    properties (Hidden=true)
        caseNum;
        vInd;
        numInd;
        p;
    end
    
    methods
        function  [obj] = NLTransLine(cfgobj,oldobj)
            %This is the constructor for the Nonlinear Transmission Line.
            %If there is 1 input, the constructor will create the class
            %instance from the variables defined in this functino.  If
            %there are 2 inputs, the first input will be unused and the
            %output will be a deep copy of the oldobj.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %cfgobj - CONFIG object
            %oldobj - NLTransLine object whose properties will be copied
            %         into the new object.
            %
            %Outputs:
            %--------
            %obj    - NLTransLine object
            %--------------------------------------------------------------
            
            if nargin == 0
                %Keep everything at defaults if no input arguments
                return;
            end
            
            % Construct a new object based on a deep copy of an old object
            % of this class by copying properties over.
            if nargin == 2 && isa(oldobj,'NLTransLine')
                props = properties(oldobj);
                for i = 1:length(props)
                    % Use Dynamic Expressions to copy the required property.
                    % For more info on usage of Dynamic Expressions, refer to
                    % the section "Creating Field Names Dynamically" in:
                    % web([docroot '/techdoc/matlab_prog/br04bw6-38.html#br1v5a9-1'])
                    obj.(props{i}) = oldobj.(props{i});
                end
                return;
            end
            
             %Store config object
            obj.config = cfgobj;
            
            %Determine mesh
            obj.mesh.node = [1:cfgobj.nNodes(1)]';
            obj.mesh.conn = [1:cfgobj.nNodes(1)-1;2:cfgobj.nNodes(1)]';
            obj.mesh.dbc = [];
            
            %Number of unknowns per node
            obj.ndim = 1;
            
            %Determine grid spacing
            coord = obj.mesh.node(2:end,1); %want coordinate vector to only include nodes that remain after removing dbc
            
            %Setup Jacobian Structure
            obj.Jstruct = JacobianStructure(obj);
            
            %Determine source term vector from the source function and the
            %coordinate vector
            if isa(cfgobj.GFunc,'function_handle')
                obj.G = cfgobj.GFunc(coord);
            else
                obj.G = feval(cfgobj.GFunc,coord);
            end
            
            %Determine ic term vector from the ic function and the
            %coordinate vector
            if isa(cfgobj.icFunc,'function_handle')
                obj.ic = cfgobj.icFunc(coord);
            else
                obj.ic = feval(cfgobj.icFunc,coord);
            end
            
            %Determine input matrix from the input matrix function
            if isa(cfgobj.BFunc,'function_handle')
                obj.B = feval(cfgobj.BFunc);
            else
                obj.B = feval(cfgobj.BFunc);
            end
            
            %Determine output matrix from the output matrix function
            if isa(cfgobj.CFunc,'function_handle')
                obj.C = feval(cfgobj.CFunc);
            else
                obj.C = feval(cfgobj.CFunc);
            end
            
            %Determine feedthrough matrix from the feedthrough matrix function
            if isa(cfgobj.DFunc,'function_handle')
                obj.D = feval(cfgobj.DFunc);
            else
                obj.D = feval(cfgobj.DFunc);
            end
        end
                
        function  [R,J] =  ResJac(obj,v,t)
        %This function calcuates the residual and jacobian of the nonlinear
        %transmission line problem
        %------------------------------------------------------------------
        %Inputs:
        %-------
        %obj - NLTrans object
        %v   - N x 1 vector containing the "solution" (or guess or
        %      approximation) at the current time
        %t   - scalar indicating the time
        %
        %Outputs:
        %--------
        %R   - N x 1 vector containing the nonlinear residual of the
        %      solution v at the current time t
        %J   - N x N vector containing the nonlinear jacobian of the
        %      solution v at the current time t
        %------------------------------------------------------------------
        
        %Determine the number of DOFs in the FOM
        N = length(v);
        
        %Comopute voltage difference between adjacent nodes and evaluate
        %the source term (and its derivative) at these potential
        %differences
        Delv = v(1:N-1,1) - v(2:N,1);
        [gDelv,dgDelv] = NLTransLine.nonlinTerm(Delv);
        [g1,dg1] = NLTransLine.nonlinTerm(v(1,1));
        
        %Set up residual (see Reweinski Thesis).  The first three terms are
        %from the tridiagonal matrix A times v, the fourth term is the G
        %vector, the fifth term is the source and the fourth is the input.
        R = [-2*v(1:end-1);-v(end)] + [0;v(1:end-1)] + [v(2:end);0] + ...
            [2-g1 - gDelv(1,1);gDelv(1:N-2,1) - gDelv(2:N-1,1);gDelv(N-1,1) - 1] + ...
            obj.G + obj.B*feval(obj.config.inFunc,t);
        
        %Set up and fill the remainder of the jacobian matrix (from
        %Reweinski Thesis).
        diagJG = [-dg1-dgDelv(1,1);...
                  -dgDelv(1:N-2,1)-dgDelv(2:N-1,1);...
                  -dgDelv(N-1,1)];
        
        J = spdiags([diagJG-2,[dgDelv+1;0],[0;dgDelv+1]],[0,-1,1],N,N);
        J(end,end) = J(end,end)+1;
        
        end
        
        function  [R,J] =  ResJacGNAT(obj,v,t)
            %This function calcuates the residual and jacobian of the
            %nonlinear transmission line
            %------------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - NLTrans object
            %U       - length(jrow) x 1 vector containing the "solution" (or
            %          guess or approximation) at the current time
            %t       - scalar indicating the time
            %
            %Outputs:
            %--------
            %R   - indN x 1 vector containing the nonlinear residual of the
            %      solution U at the current time t
            %J   - vector containing the nonlinear jacobian of the solution
            %      U at the current time t at the rows indN
            %--------------------------------------------------------------
            
            switch obj.caseNum
                case 1
                    v1 = v(obj.vInd); v2 = v(obj.vInd+1); v3 = v(obj.vInd+2);
                    [g1,dg1] = NLTransLine.nonlinTerm(v(1));
                    [g1Del,dg1Del] = NLTransLine.nonlinTerm( v(1) - v(2) );
                    [gnDel,dgnDel] = NLTransLine.nonlinTerm(v(end-1) - v(end));
                    [gDelv1,dgDelv1] = NLTransLine.nonlinTerm( v1 - v2 );
                    [gDelv2,dgDelv2] = NLTransLine.nonlinTerm( v2 - v3 );
                    R =  [-2*v(1) + v(2) + 2 - g1 - g1Del; v1 - 2*v2 + v3 + gDelv1 - gDelv2; v(end-1) - v(end) + gnDel - 1] + obj.G + obj.B*feval(obj.config.inFunc,t);
                    J = [-2;1;repmat([1; -2; 1],obj.numInd-2,1);1;-1] +[ -dg1 - dg1Del; dg1Del ; reshape([dgDelv1, -dgDelv1 - dgDelv2, dgDelv2]',3*(obj.numInd-2),1); dgnDel; -dgnDel];
                case 2
                    v1 = v(obj.vInd); v2 = v(obj.vInd+1); v3 = v(obj.vInd+2);
                    [g1,dg1] = NLTransLine.nonlinTerm(v(1));
                    [g1Del,dg1Del] = NLTransLine.nonlinTerm( v(1) - v(2) );
                    [gDelv1,dgDelv1] = NLTransLine.nonlinTerm( v1 - v2 );
                    [gDelv2,dgDelv2] = NLTransLine.nonlinTerm( v2 - v3 );
                    R =  [-2*v(1) + v(2) + 2 - g1 - g1Del; v1 - 2*v2 + v3 + gDelv1 - gDelv2] + obj.G + obj.B*feval(obj.config.inFunc,t);
                    J = [-2;1;repmat([1; -2; 1],obj.numInd-1,1)] +[ -dg1 - dg1Del; dg1Del ; reshape([dgDelv1, -dgDelv1 - dgDelv2, dgDelv2],3*(obj.numInd-1),1)];
                case 3
                    v1 = v(obj.vInd); v2 = v(obj.vInd+1); v3 = v(obj.vInd+2);
                    [gnDel,dgnDel] = NLTransLine.nonlinTerm(v(end-1) - v(end));
                    [gDelv1,dgDelv1] = NLTransLine.nonlinTerm( v1 - v2 );
                    [gDelv2,dgDelv2] = NLTransLine.nonlinTerm( v2 - v3 );
                    R =  [v1 - 2*v2 + v3 + gDelv1 - gDelv2; v(end-1) - v(end) + gnDel - 1] + obj.G + obj.B*feval(obj.config.inFunc,t);
                    J = [repmat([1; -2; 1],obj.numInd-1,1);1;-1] +[reshape([dgDelv1, -dgDelv1 - dgDelv2, dgDelv2],3*(obj.numInd-1),1); dgnDel; -dgnDel];
                case 4
                    v1 = v(obj.vInd); v2 = v(obj.vInd+1); v3 = v(obj.vInd+2);
                    [gDelv1,dgDelv1] = NLTransLine.nonlinTerm( v1 - v2 );
                    [gDelv2,dgDelv2] = NLTransLine.nonlinTerm( v2 - v3 );
                    R = v1 - 2*v2 + v3 + gDelv1 - gDelv2 + obj.G + obj.B*feval(obj.config.inFunc,t);
                    J = repmat([1; -2; 1],obj.numInd,1) + reshape([dgDelv1, -dgDelv1 - dgDelv2, dgDelv2],3*obj.numInd,1);
            end
        end
        
        function  [R,J] = ResJacTimeDer(obj,U,~)
            R = U;
            J = eye(size(U));
        end
        
        function  [R,J] = ResJacGNATTimeDer(obj,U,~,~,~,jrow,jdiag)
            R = U(jrow(logical(jdiag)));
            J = jdiag;
        end
        
        function  [R,J] = ResJacWithoutInput(obj,U,t)
            %This function returns the residual and jacobian without the
            %input term contributions.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj - NLTransLine object
            %
            %Outputs:
            %--------
            %Radj - residual without input term contribution
            %Jadj - jacobian without input term contribution
            %--------------------------------------------------------------
            
            [R,J] = obj.ResJac(U,t);
            R = R - obj.B*feval(obj.config.inFunc,t);
            
        end
                
        function  [ind] = node2ind(obj,nodenum)
            %This function takes a node in the mesh and maps it to an index
            %in the full state vector for the nonlinear transmission line
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - NLTransLine object
            %nodenum - scalar specifying node numbers
            %
            %Outputs:
            %--------
            %ind     - vector specifying the corresponding indices in the
            %          full state vector
            %--------------------------------------------------------------
            
            if nodenum == 1
                ind = nan;
            else
                ind = nodenum - 1;
            end
        end
        
        function  [JacStruct] = JacobianStructure(obj)
            %This function forms a boolean matrix whose true entries
            %correspond to nonzero entries in the jacobian for the
            %nonlinear transmission line
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj       - NLTransLine object
            %
            %Outputs:
            %--------
            %JacStruct - ndof x ndof boolean matrix whose true values
            %            correspond to nonzero entries in the jacobian
            %--------------------------------------------------------------
            
            N = obj.config.ndof;
            JacStruct = spdiags(ones(N,3),-1:1,N,N);            
        end
        
        function  [probGNAT] = createCopy4GNAT(obj,gnat)
            %This function makes a copy of the current object and only
            %copies the source term and input matrix evaluated at indN
            %rows.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - NLTransLine object
            %gnat    - GNAT object
            %
            %Outputs:
            %--------
            %probGNAT - new "reduced" NLTransLine object
            %--------------------------------------------------------------
            
            probGNAT = NLTransLine([],obj);
            
            temp = probGNAT.config.inFunc;
            temp1 = probGNAT.config.ndof;
            temp2 = probGNAT.config.form;
            probGNAT.config = [];
            probGNAT.config.inFunc = temp;
            probGNAT.config.ndof = temp1;
            probGNAT.config.form = temp2;
            probGNAT.mesh   = [];
            probGNAT.ndim = [];
            probGNAT.Jstruct = [];
            
            probGNAT.numInd = length(gnat.sampleInd);
            if gnat.sampleInd(1) == 1 && gnat.sampleInd(end) == obj.config.ndof
                probGNAT.caseNum = 1;
                probGNAT.vInd = gnat.irstart(2:end-2);
            elseif gnat.sampleInd(1) == 1
                probGNAT.caseNum = 2;
                probGNAT.vInd = gnat.irstart(2:end-1);
            elseif gnat.sampleInd(end) == obj.config.ndof
                probGNAT.caseNum = 3;
                probGNAT.vInd = gnat.irstart(2:end-2);
            else
                probGNAT.caseNum = 4;
                probGNAT.vInd = gnat.irstart(1:end-1);
            end
            
            if length(obj.B) > 1
                probGNAT.B = obj.B(gnat.sampleInd,1);
            end
            
            if length(obj.G) > 1
                probGNAT.G = obj.G(gnat.sampleInd,1);
            end
            
            if length(obj.ic) > 1
                probGNAT.ic = obj.ic(unique(gnat.jrow),1);
            end
        end
        
        function  [probGNAT] = createCopy4GNATloc(obj,gnat)
            %This function makes a copy of the current object and only
            %copies the source term and input matrix evaluated at indN
            %rows.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - NLTransLine object
            %gnat    - locGNAT object
            %
            %Outputs:
            %--------
            %probGNAT - new "reduced" NLTransLine object
            %--------------------------------------------------------------
            
            for kk = 1:gnat.nBases
                probGNAT(kk) = NLTransLine([],obj);
                
                temp = probGNAT(kk).config.inFunc;
                temp1 = probGNAT(kk).config.ndof;
                temp2 = probGNAT(kk).config.form;
                probGNAT(kk).config = [];
                probGNAT(kk).config.inFunc = temp;
                probGNAT(kk).config.ndof = temp1;
                probGNAT(kk).config.form = temp2;
                probGNAT(kk).mesh   = [];
                probGNAT(kk).ndim = [];
                probGNAT(kk).Jstruct = [];
                
                probGNAT(kk).numInd = length(gnat.sampleInd{kk,1});
                if gnat.sampleInd{kk,1}(1) == 1 && gnat.sampleInd{kk,1}(end) == obj.config.ndof
                    probGNAT(kk).caseNum = 1;
                    probGNAT(kk).vInd = gnat.irstart{kk,1}(2:end-2);
                elseif gnat.sampleInd{kk,1}(1) == 1
                    probGNAT(kk).caseNum = 2;
                    probGNAT(kk).vInd = gnat.irstart{kk,1}(2:end-2);
                elseif gnat.sampleInd{kk,1}(end) == obj.config.ndof
                    probGNAT(kk).caseNum = 3;
                    probGNAT(kk).vInd = gnat.irstart{kk,1}(2:end-2);
                else
                    probGNAT(kk).caseNum = 4;
                    probGNAT(kk).vInd = gnat.irstart{kk,1}(1:end-1);
                end
                
                if length(obj.B) > 1
                    probGNAT(kk).B = obj.B(gnat.sampleInd{kk,1},1);
                end
                
                if length(obj.G) > 1
                    probGNAT(kk).G = obj.G(gnat.sampleInd{kk,1},1);
                end
                
                if length(obj.ic) > 1
                    probGNAT(kk).ic = obj.ic(unique(gnat.jrow{kk,1}),1);
                end
            end
        end
        
        function   [axOUT] = problemPaperPlot(obj,modelobj,modelAuxobj,pstr,axIN)
            %This function generates the Nonlinear Transmission line plot
            %from the Rewienski thesis.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj         - Problem object
            %modelobj    - FOM, ROM, GNAT, or TPWL object
            %modelAuxobj - ROM object if modelobj is of class GNAT or
            %              locGNAT (unused otherwise)
            %pstr        - cell array indicating plot object formatting
            %axIN        - axis handle to place new plots
            %
            %Outputs:
            %--------
            %axOUT       - axis handle that new plots were placed
            %--------------------------------------------------------------
            
            %Extract reconstructed state vectors (with DBCs)
            svF = obj.returnPosSV(modelobj,modelAuxobj);
            
            %Extract mesh spacing
            x = obj.mesh.node;
            tvec = modelobj.time.T(1):modelobj.time.dt:modelobj.time.T(2);
            
            %Set up axes
            if isempty(axIN)
                ax = axes;
                title('Voltage Time History at Node 1','fontsize',16,'interpreter','latex');
                xlabel('Time [s]','fontsize',14,'interpreter','latex');
                ylabel('Voltage at Node 1 [V]','fontsize',14,'interpreter','latex');
            else
                ax = axIN;
            end
            
            plot(ax,tvec,obj.C*svF,pstr{:});
            axOUT = ax;
        end
        
        function  [X,Y,svDBC] = returnPosSV(obj,modelobj,modelAuxobj,~,flag,basenum)
            %This function takes a state vector from a simulation and
            %returns the state vector with the Dirichlet boundary
            %conditions.  A vector of positions (including those with DBCs)
            %is also returned.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj         - Problem object
            %modelobj    - FOM, ROM, GNAT, or TPWL object
            %modelAuxobj - ROM object if modelobj is of class GNAT, empty
            %              otherwise
            %flag        - string ('pod','res','jac','clusCenter','sv')
            %basenum     - basis number to use (only needed for phi)
            %
            %Outputs:
            %--------
            %X           - x coordinate positions
            %Y           - y coordinate positions (empty for this problem)
            %svDBC       - state vector matrix with DBCs
            %--------------------------------------------------------------
            
            %Extract data (state vectors, cluster centers,
            %residual/jacobians, pod vectors)
            switch flag
                case 'pod'
                    svF = modelobj.phi(:,:,basenum);
                case 'res'
                    svF = modelobj.res;
                case 'jac'
                    svF = modelobj.jac;
                case 'clusCenter'
                    svF = modelobj.clusCenter;
                case 'sv'
                    svF = modelobj.reconstructFullState(modelAuxobj);
            end
            
            X = obj.mesh.node;
            Y = [];
        end
        
        function   [] = problemAnimate(obj,modelobj,modelAuxobj,pstr,~,freq)
            %This function generates the animation for the Nonlinear
            %Transmission line.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj      - Problem object
            %modelobj - FOM, ROM, GNAT, or TPWL object
            %modelAuxobj - ROM object if modelobj is of class GNAT or
            %              locGNAT (unused otherwise)
            %pstr        - cell array indicating plot object formatting
            %freq        - scalar indicating how often to plot sv in
            %              animation
            %
            %Outputs:
            %--------
            %There are no outputs for this function.  It generates an
            %animiation.
            %--------------------------------------------------------------
            
            if isempty(pstr)
                pstr = {'k.','MarkerSize',6};
            end
            
            %Extract reconstructed state vectors (with DBCs)
            svF = obj.returnPosSV(modelobj,modelAuxobj);
            
            %Extract mesh spacing
            tvec = modelobj.time.T(1):modelobj.time.dt:modelobj.time.T(2);
            
            %Set up axes
            set(gcf,'position',[439   566   714   540]);
            ax = axes;
            title('Voltage Time History at Node 1','fontsize',16,'fontname','georgia');
            xlabel('Time [s]','fontsize',14,'interpreter','latex');
            ylabel('Voltage at Node 1 [V]','fontsize',14,'fontname','georgia');
            set(ax,'nextplot','replacechildren','ylim',[ymin,ymax]);
            
            output = obj.C*svF;
            for i = 1:freq:length(tvec)
                plot(ax,tvec(i),output(1,i),pstr{:});
            end
        end
        
        function  [varargout] = setPrevSV(varargin)
        end
    end
    
    methods (Static)
        function [g,dg] = nonlinTerm(U)
            %This function defines the nonlinear term of the nonlinear ODE.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %U   - state vector
            %
            %Outputs:
            %--------
            %g   - vector containing the evaluation of the nonlinear term
            %      at each component of the state vector
            %dg  - vector containing the evaluation of the derivative of
            %      the nonlinear term at each component of the state vector
            %--------------------------------------------------------------
            
            g = exp(40*U);
            dg = 40*exp(40*U);
        end
    end
end
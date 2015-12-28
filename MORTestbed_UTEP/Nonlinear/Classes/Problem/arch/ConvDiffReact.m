classdef ConvDiffReact < handle
    
    properties (SetAccess = private, GetAccess = public)
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
        dx;
        
        Pe_m = 5;
        Pe_h = 5;
        gamma = 25;
        BB = 0.5;
        DD = 0.165;
        beta = 2.5;
        th0 = 1;
        thIn;
        yIn;
    end
    
    properties (Hidden=true)
        Uind;
        ind1;
        lengthJROW;
    end
    
    methods
        function  [obj] = ConvDiffReact(cfgobj,oldobj)
            %This is the constructor for the 1D Convection-Diffusion
            %Reaction Equation.  If there is 1 input, the constructor will
            %create the class instance from the variables defined in this
            %function.  If there are 2 inputs, the first input will be
            %unused and the output will be a deep copy of the oldobj.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %cfgobj - CONFIG object
            %oldobj - ConvDiffReact object whose properties will be copied
            %         into the new object.
            %
            %Outputs:
            %--------
            %obj    - ConvDiffReact object
            %--------------------------------------------------------------
            
            if nargin == 0
                %Keep everything at defaults if no input arguments
                return;
            end
            
            % Construct a new object based on a deep copy of an old object
            % of this class by copying properties over.
            if nargin == 2 && isa(oldobj,'ConvDiffReact')
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
            
            %Store number of DOFs (2 Neumann BCs - analytically removed from unknowns)
            ndof = 2*(cfgobj.nNodes(1) - 2);
            
            %Read altFile if it exists
            if isempty(obj.config.altFile)
                [obj.Pe_m,obj.Pe_h,obj.gamma,obj.BB,obj.DD,obj.beta,obj.th0] = feval(obj.config.altFile,cfgobj.nNodes(1));
            end
            
            %Determine grid spacing
            obj.dx = (cfgobj.DLim(1,2) - cfgobj.DLim(1,1))/(ndof+1);
            
            %Determine mesh
            obj.mesh.node = linspace(cfgobj.DLim(1,1),cfgobj.DLim(1,2),cfgobj.nNodes(1))';
            obj.mesh.conn = [1:cfgobj.nNodes(1)-1;2:cfgobj.nNodes(1)]';
            obj.mesh.dbc = [1,true;length(obj.mesh.node),true];
            
            %Number of unknowns per node
            obj.ndim = 2;
            
            %Determine grid spacing
            obj.dx = obj.mesh.node(2,1) - obj.mesh.node(1,1);
            coord = obj.mesh.node(2:end-1,1); %want coordinate vector to only include nodes that remain after removing dbc
            
            %Setup Jacobian Structure
            JacobianStructure(obj);
            
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
            %obj.B = (0.5/obj.dx)*[1;zeros(ndof-1,1)];
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
        
        function  [R] = ResidualOnly(obj,U,t)
            %This function calcuates the residual of the
            %Convection-Diffusion Equation
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj - ConvDiffReact object
            %U   - N x 1 vector containing the "solution" (or guess or
            %      approximation) at the current time
            %t   - scalar indicating the time
            %
            %Outputs:
            %--------
            %R   - N x 1 vector containing the nonlinear residual of the
            %      solution
            %      U at the current time t
            %--------------------------------------------------------------
            
            %Determine the number of DOFs in the FOM
            R = zeros(obj.config.ndof,1);
            
            %Define constants
            c1 = 1/(obj.dx*obj.dx*obj.Pe_m) - 1/(2*obj.dx);
            c2 = 1/(obj.dx*obj.dx*obj.Pe_m) + 1/(2*obj.dx);
            d1 = 1/(obj.dx*obj.dx*obj.Pe_h) - 1/(2*obj.dx);
            d2 = 1/(obj.dx*obj.dx*obj.Pe_h) + 1/(2*obj.dx);
            
            %Determine nonlinear contributions to residual
            NLterm = obj.nonlinFunc(U);
            
            %Extract y and theta
            y = U(1:2:end);
            th = U(2:2:end);
            
            %Fill the entries of the residual appropriately
            R(1,1) = c1*y(2) + (-2/(obj.Pe_m*obj.dx*obj.dx) + 1/(obj.dx^3*obj.Pe_m*(obj.Pe_m+1/obj.dx)) + 1/(2*obj.dx*obj.dx*(obj.Pe_m+1/obj.dx)))*y(1) - obj.DD*NLterm(1) + c2*obj.Pe_m/(obj.Pe_m+1/obj.dx);
            R(2,1) = d1*th(2) + (-2/(obj.Pe_h*obj.dx*obj.dx) + 1/(obj.dx^3*obj.Pe_h*(obj.Pe_h+1/obj.dx)) + 1/(2*obj.dx*obj.dx*(obj.Pe_h+1/obj.dx))-obj.beta)*y(1) + obj.BB*obj.DD*NLterm(1) + d2*obj.Pe_h/(obj.Pe_h+1/obj.dx) + obj.beta*obj.th0;
            R(3:2:end-3,1) = c1*y(3:end) + c2*y(1:end-2) - 2/(obj.dx*obj.dx*obj.Pe_m)*y(2:end-1) - obj.DD*NLterm(2:end-1);
            R(4:2:end-2,1) = d1*th(3:end) + d2*th(1:end-2) - (2/(obj.dx*obj.dx*obj.Pe_h) + obj.beta)*th(2:end-1) + obj.BB*obj.DD*NLterm(2:end-1) + obj.beta*obj.th0;
            R(end-1,1) = c1*y(end-1) - c2*y(end) - obj.DD*NLterm(end);
            R(end,1)   = d1*th(end-1) - (d2+obj.beta)*th(end) + obj.BB*obj.DD*NLterm(end) + obj.beta*obj.th0;
        end
        
        function  [R,J] = ResJac(obj,U,t)
            %This function calcuates the residual and jacobian of
            %Convection-Diffusion Equation
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj - ConvDiffReact object
            %U   - N x 1 vector containing the "solution" (or guess or
            %      approximation) at the current time
            %t   - scalar indicating the time
            %
            %Outputs:
            %--------
            %R   - N x 1 vector containing the nonlinear residual of the solution
            %      U at the current time t
            %J   - N x N vector containing the nonlinear jacobian of the solution
            %      U at the current time t
            %--------------------------------------------------------------
            
            %Compute residual
            R = obj.ResidualOnly(U,t);
            
            %Fill the entries of the jacobian using finite differences
            n = length(R);
            J = zeros(n,n);
            
            for j = 1:n
                vec = U;
                vec(j) = vec(j) + eps;
                J(:,j) = (1/eps)*(obj.ResidualOnly(vec,t) - R);
            end
        end
        
        function  [R,J] = ResJacGNAT(obj,U,t)
            %This function calcuates the residual and jacobian of Burger's
            %Equation
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - OneDBurgers object
            %U       - N x 1 vector containing the "solution" (or guess or
            %          approximation) at the current time
            %t       - scalar indicating the time
            %indN    - vector containing the sample indicies
            %irstart - irstart vector from GNAT class
            %
            %Outputs:
            %--------
            %R   - indN x 1 vector containing the nonlinear residual of the
            %      solution U at the current time t
            %J   - vector containing the nonlinear jacobian of the solution
            %      U at the current time t at the rows indN
            %--------------------------------------------------------------
            
            %             M = irstart(end)-1; %This is equal to numel(jrow)
            if obj.ind1 == 1
                %This is done first so the entire array is initialized
                R = [-U(1,1)*U(1,1)/(2*obj.dx);...
                    (U(obj.Uind,1).*U(obj.Uind,1)-U(obj.Uind+1,1).*U(obj.Uind+1,1))/(2*obj.dx)] + ...
                    obj.G + obj.B*feval(obj.config.inFunc,t);
                
                J = [-U(1,1)/obj.dx; reshape([U(obj.Uind,1)'; -U(obj.Uind+1,1)'],obj.lengthJROW-1,1)/(obj.dx)];
                
                %                 %This is done first so the entire array is initialized
                %                 R(2:numel(indN),1) = (U(irstart(2:end-1),1).*U(irstart(2:end-1),1)-U(irstart(2:end-1)+1,1).*U(irstart(2:end-1)+1,1))/(2*obj.dx);
                %                 J(2:M,1) = reshape([U(irstart(2:end-1),1)'; -U(irstart(2:end-1)+1,1)'],M-1,1)/(obj.dx);
                %
                %                 %This is hard-coded as such because I know that the indN will
                %                 %be sorted!
                %                 R(1,1) = -U(1,1)*U(1,1)/(2*obj.dx);
                %                 J(1,1) = -U(1,1)/obj.dx;
            else
                R = (U(obj.Uind,1).*U(obj.Uind,1)-U(obj.Uind+1,1).*U(obj.Uind+1,1))/(2*obj.dx) + obj.G + obj.B*feval(obj.config.inFunc,t);
                J = reshape([U(obj.Uind,1)'; -U(obj.Uind+1,1)'],obj.lengthJROW,1)/(obj.dx);
                return;
            end
            
            %             R = R + obj.G + obj.B*feval(obj.config.inFunc,t);
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
            %obj - ConvDiffReact object
            %
            %Outputs:
            %--------
            %Radj - residual without input term contribution
            %Jadj - jacobian without input term contribution
            %--------------------------------------------------------------
            
            [R,J] = obj.ResJac(U,t);
            R = R - obj.B*feval(obj.config.inFunc,t);
            
        end
        
        function  [f] = nonlinFunc(obj,U)
            %This function returns the nonlinear term of the nonlinear
            %convection-diffusion reaction equation evaluated at U.
            %---------
            %Inputs:
            %-------
            %obj - ConvDiffReact object
            %U   - N x 1 vector containing the "solution" (or guess or
            %      approximation) at the current time
            %
            %Outputs:
            %--------
            %f   - vector containing the nonlinear term evaluated at U
            %--------------------------------------------------------------
            
            f = U(1:2:end).*exp(obj.gamma - obj.gamma./U(2:2:end));
        end
        
        function  [ind] = node2ind(obj,nodenum)
            %This function takes a node in the mesh and maps it to an index
            %in the full state vector
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - OneDBurgers object
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
        
        function  [] = JacobianStructure(obj)
            %This function forms a boolean matrix whose true entries
            %correspond to nonzero entries in the jacobian
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - OneDBurgers object
            %
            %Outputs:
            %--------
            %There are no outputs.  It store a ndof x ndof boolean matrix
            %whose true values correspond to nonzero entries in the
            %jacobian in the class instance.
            %--------------------------------------------------------------
            
            N = obj.config.ndof;
            temp = diag(ones(N,1),0) + diag(ones(N-1,1),-1);
            obj.Jstruct = sparse(temp);
            
        end
        
        function  [probGNAT] = createCopy4GNAT(obj,gnat)
            %This function makes a copy of the current object and only
            %copies the source term and input matrix evaluated at indN
            %rows.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - OneDBurgers object
            %gnat    - gnat object
            %
            %Outputs:
            %--------
            %probGNAT - new "reduced" OneDBurgers object
            %--------------------------------------------------------------
            
            probGNAT = OneDBurgers([],obj);
            
            temp = probGNAT.config.inFunc;
            probGNAT.config = [];
            probGNAT.config.inFunc = temp;
            probGNAT.mesh   = [];
            probGNAT.ndim = [];
            probGNAT.Jstruct = [];
            
            probGNAT.ind1 = gnat.sampleInd(1);
            if probGNAT.ind1 == 1
                probGNAT.Uind = gnat.irstart(2:end-1);
            else
                probGNAT.Uind = gnat.irstart(1:end-1);
            end
            probGNAT.lengthJROW = gnat.irstart(end)-1;
            
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
            %obj     - OneDBurgers object
            %gnat    - gnat object
            %
            %Outputs:
            %--------
            %probGNAT - array of new "reduced" OneDBurgers object
            %--------------------------------------------------------------
            
            for kk = 1:gnat.nBases
                if isempty(gnat.sampleInd{kk,1})
                    continue;
                end
                probGNAT(kk) = OneDBurgers([],obj);
                
                temp = probGNAT(kk).config.inFunc;
                %                 probGNAT(kk).ic = obj.ic(jrowUnion);
                probGNAT(kk).config = [];
                probGNAT(kk).config.inFunc = temp;
                probGNAT(kk).mesh   = [];
                probGNAT(kk).ndim = [];
                probGNAT(kk).Jstruct = [];
                
                probGNAT(kk).ind1 = gnat.sampleInd{kk,1}(1);
                if probGNAT(kk).ind1 == 1
                    probGNAT(kk).Uind = gnat.irstart{kk,1}(2:end-1);
                else
                    probGNAT(kk).Uind = gnat.irstart{kk,1}(1:end-1);
                end
                probGNAT(kk).lengthJROW = gnat.irstart{kk,1}(end)-1;
                
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
            %This function generates the 1D Burger's Equation plot from the
            %Rewienski thesis.
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
            %axOUT    - axis handle that new plots were placed
            %--------------------------------------------------------------
            
            %Store the times instances to plot the solution (based on
            %Reweinski thesis)
            Tplots = [2.5,10,20,30,40,50];
            
            %Extract reconstructed state vectors (with DBCs)
            svF = obj.returnPosSV(modelobj,modelAuxobj);
            
            %Extract mesh spacing
            x = obj.mesh.node;
            n = length(Tplots);
            tvec = modelobj.time.T(1):modelobj.time.dt:modelobj.time.T(2);
            ind = zeros(1,2);
            
            %Set up axes
            if isempty(axIN)
                ax = axes;
                title('Shock Propagation Through Fluid','fontsize',16,'interpreter','latex');
                xlabel('Position','fontsize',14,'interpreter','latex');
                ylabel('Conserved Quantity','fontsize',14,'interpreter','latex');
            else
                ax = axIN;
            end
            set(ax,'nextplot','add');
            
            for i = 1:n
                %Determine indices in svF that correspond to entries in
                %Tplots
                if isempty(Tplots(i) == tvec)
                    [~,temp] = sort(abs(tvec - Tplots(i)));
                    ind = temp(1:2);
                else
                    ind(1,:) = find(tvec == Tplots(i));
                end
                U = 0.5*(svF(:,ind(1,1)) + svF(:,ind(1,2)));
                
                plot(ax,x,U,pstr{:});
            end
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
            n = size(svF,2);
            
            t = modelobj.time.T(1):modelobj.time.dt:modelobj.time.T(2);
            u = feval(obj.config.inFunc,t);
            
            svDBC = zeros(obj.config.ndof+1,n);
            svDBC(1,:) = sqrt(u);
            svDBC(2:end,:) = svF;
            
            X = obj.mesh.node;
            Y = [];
            
        end
        
        function   [] = problemAnimate(obj,modelobj,modelAuxobj,pstr,~,freq)
            %This function generates a problem specific 1D Burger's
            %Equation plot.
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
            
            %Extract reconstructed state vectors (with DBCs) and mesh
            %spacing
            [x,~,svF] = obj.returnPosSV(modelobj,modelAuxobj);
            
            tvec = modelobj.time.T(1):modelobj.time.dt:modelobj.time.T(2);
            
            %Determine ylimits
            ymin = min(min(svF)); ymax = max(max(svF));
            
            set(gcf,'position',[439   566   714   540]);
            ax = axes;
            title(ax,'Shock Propagation Through Fluid','fontsize',16,'fontname','georgia');
            xlabel(ax,'Position','fontsize',14,'fontname','georgia');
            ylabel(ax,'Conserved Quantity','fontsize',14,'fontname','georgia');
            set(ax,'nextplot','replacechildren','ylim',[ymin,ymax]);
            for i = 1:freq:length(tvec)
                plot(ax,x,svF(:,i),pstr{:}); pause(0.001);
            end
        end
    end
end

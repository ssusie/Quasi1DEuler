classdef FHN
    
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
        dx;
        
        b       = 0.5;
        epsilon = 0.015;
        gamma   = 2;
    end
    
    properties (Hidden = true)
        caseNum;
        Uind;
        indV;
        numVeq;
        numVvar;
        numWeq;
        indW;
        p;
    end
    
    methods
        function  [obj] = FHN(cfgobj,oldobj)
            %This is the constructor for the FitzHugh-Nagumo problem.
            %If there is 1 input, the constructor will create the class
            %instance from the variables defined in this functino.  If
            %there are 2 inputs, the first input will be unused and the
            %output will be a deep copy of the oldobj.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %cfgobj - CONFIG object
            %oldobj - FHN object whose properties will be copied
            %         into the new object.
            %
            %Outputs:
            %--------
            %obj    - FHN object
            %--------------------------------------------------------------
            
            if nargin == 0
                %Keep everything at defaults if no input arguments
                return;
            end
            
            % Construct a new object based on a deep copy of an old object
            % of this class by copying properties over.
            if nargin == 2 && isa(oldobj,'FHN')
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
            
            %Extract parameters (if they were specified)
            if length(cfgobj.param) == 1
                obj.b = cfgobj.param;
            elseif length(cfgobj.param) == 2
                obj.b = cfgobj.param(1);
                obj.epsilon = cfgobj.param(2);
            elseif length(cfgobj.param) == 3
                obj.b = cfgobj.param(1);
                obj.epsilon = cfgobj.param(2);
                obj.gamma = cfgobj.param(3);
            elseif length(cfgobj.param) > 3
                error('param in the .cfg file have at most 3 columns, where column 1 is b, column 2 is epsilon, and column 3 is gamma');
            end
                        
            %Determine grid spacing
            obj.dx = (cfgobj.DLim(1,2) - cfgobj.DLim(1,1))/(cfgobj.nNodes(1)-1);
            
            %Determine mesh
            obj.mesh.node = linspace(cfgobj.DLim(1,1),cfgobj.DLim(1,2),cfgobj.nNodes(1))';
            obj.mesh.conn = [1:cfgobj.nNodes(1)-1;2:cfgobj.nNodes(1)]';
            obj.mesh.dbc = [];
            
            %Number of unknowns per node
            obj.ndim = 2;
            
            %Determine grid spacing
            obj.dx = obj.mesh.node(2,1) - obj.mesh.node(1,1);
            coord = obj.mesh.node; %want coordinate vector to only include nodes that remain after removing dbc
            
            %Setup Jacobian Structure
            obj.Jstruct = JacobianStructure(obj);
            
            %Determine source term vector from the source function and the
            %coordinate vector
            if isa(cfgobj.GFunc,'function_handle')
                obj.G = [(1/obj.epsilon)*cfgobj.GFunc(coord);cfgobj.GFunc(coord)];
            else
                obj.G = [(1/obj.epsilon)*feval(cfgobj.GFunc,coord);feval(cfgobj.GFunc,coord)];
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
                obj.B = (obj.epsilon/obj.dx)*feval(cfgobj.BFunc);
            else
                obj.B = (obj.epsilon/obj.dx)*feval(cfgobj.BFunc);
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
        
        function  [R,J] =  ResJac(obj,U,t)
        %This function calcuates the residual and jacobian of Burger's Equation
        %--------------------------------------------------------------------------
        %Inputs:
        %-------
        %obj - FOM object
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
        %--------------------------------------------------------------------------
        
        %Determine the number of DOFs in the FOM
        N = obj.config.ndof/2;
        v = U(1:N,1); w = U(N+1:end,1);
        e_dx2 = obj.epsilon/(obj.dx^2);
        einv = 1/obj.epsilon;
        
        R = [e_dx2*[0.5*(v(3)-v(1));v(1:end-2) - 2*v(2:end-1) + v(3:end);0.5*(v(end-2)-v(end))];obj.b*v] + ...
            [(-einv)*w;-obj.gamma*w] + [einv*v.*(v-0.1).*(1-v);zeros(N,1)] + obj.G +...
            obj.B*feval(obj.config.inFunc,t);
        
%J = [diag(e_dx2*[-0.5;-2*ones(N-2,1);-0.5],0) + diag(e_dx2*[0;ones(N-2,1)],1) + diag(e_dx2*[ones(N-2,1);0],-1),diag((-einv)*ones(N,1),0);...
%     diag(obj.b*ones(N,1),0),diag((-obj.gamma)*ones(N,1),0)] + einv*diag([-3*v.^2 + 2.2*v - 0.1;zeros(N,1)]);
        J = [spdiags([e_dx2*[ones(N-2,1);0;0],e_dx2*[-0.5;-2*ones(N-2,1);-0.5] + einv*(-3*v.*v + 2.2*v - 0.1),e_dx2*[0;0;ones(N-2,1)]],-1:1,N,N),spdiags(-einv*ones(N,1),0,N,N); ...
            spdiags(obj.b*ones(N,1),0,N,N),spdiags(-obj.gamma*ones(N,1),0,N,N)];
        J(1,3) = 0.5*e_dx2;
        J(N,N-2) = 0.5*e_dx2;
        end
        
        function  [R,J] =  ResJacGNAT(obj,U,t)
            %This function calcuates the residual and jacobian of FHN
            %problem at the specified indices.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - FHN object
            %U       - N x 1 vector containing the "solution" (or guess or
            %          approximation) at the current time
            %t       - scalar indicating the time
            %
            %Outputs:
            %--------
            %R   - vector containing the nonlinear residual of the solution
            %      U at the current time t at the sample entries requested
            %J   - vector containing the nonlinear jacobian of the solution
            %      U at the current time t at the rows requested
            %--------------------------------------------------------------
            
            e_dx2 = obj.epsilon/(obj.dx^2);
            einv = 1/obj.epsilon;
            
            switch obj.caseNum
                case 1
                    v = U(obj.indV); w = U(~obj.indV);
                    
                    R = [0.5*e_dx2*(-v(1)+v(2)) + einv*(-w(1) + FHN.f(v(1)));...
                        e_dx2*(U(obj.Uind) - 2*U(obj.Uind+1) + U(obj.Uind+2)) + einv*(FHN.f(U(obj.Uind+1))-U(obj.Uind+3));...
                        0.5*e_dx2*(v(obj.numVvar-1)-v(obj.numVvar)) + einv*(-w(obj.indW(1)-1) + FHN.f(v(obj.numVvar)));...
                        obj.b*v(obj.numVvar+1:end) - obj.gamma*w(obj.indW)] + obj.G + obj.B*feval(obj.config.inFunc,t);
                    
                    J = [-0.5*e_dx2+einv*FHN.df(v(1));0.5*e_dx2;-einv;reshape([e_dx2*ones(obj.numVeq,1),-2*e_dx2+einv*FHN.df(U(obj.Uind+1)),e_dx2*ones(obj.numVeq,1),-einv*ones(obj.numVeq,1)]',4*obj.numVeq,1);0.5*e_dx2;-0.5*e_dx2+einv*FHN.df(v(obj.numVvar));-einv;repmat([obj.b;-obj.gamma],obj.numWeq,1)];
                case 2
                    v = U(obj.indV); w = U(~obj.indV);
                    
                    R = [0.5*e_dx2*(-v(1)+v(2)) + einv*(-w(1) + FHN.f(v(1)));...
                        e_dx2*(U(obj.Uind) - 2*U(obj.Uind+1) + U(obj.Uind+2)) + einv*(FHN.f(U(obj.Uind+1))-U(obj.Uind+3));...
                        obj.b*v(obj.numVvar+1:end) - obj.gamma*w(obj.indW)] + obj.G + obj.B*feval(obj.config.inFunc,t);
                    
                    J = [-0.5*e_dx2+einv*FHN.df(v(1));0.5*e_dx2;-einv;reshape([e_dx2*ones(obj.numVeq,1),-2*e_dx2+einv*FHN.df(U(obj.Uind+1)),e_dx2*ones(obj.numVeq,1),-einv*ones(obj.numVeq,1)]',4*obj.numVeq,1);repmat([obj.b;-obj.gamma],obj.numWeq,1)];
                case 3
                    v = U(obj.indV); w = U(~obj.indV);
                    
                    R = [e_dx2*(U(obj.Uind) - 2*U(obj.Uind+1) + U(obj.Uind+2)) + einv*(FHN.f(U(obj.Uind+1))-U(obj.Uind+3));...
                        0.5*e_dx2*(v(obj.numVvar-1)-v(obj.numVvar)) + einv*(-w(obj.indW(1)-1) + FHN.f(v(obj.numVvar)));...
                        obj.b*v(obj.numVvar+1:end) - obj.gamma*w(obj.indW)] + obj.G + obj.B*feval(obj.config.inFunc,t);
                    
                    J = [reshape([e_dx2*ones(obj.numVeq,1),-2*e_dx2+einv*FHN.df(U(obj.Uind+1)),e_dx2*ones(obj.numVeq,1),-einv*ones(obj.numVeq,1)]',4*obj.numVeq,1);0.5*e_dx2;-0.5*e_dx2+einv*FHN.df(v(obj.numVvar));-einv;repmat([obj.b;-obj.gamma],obj.numWeq,1)];
                case 4
                    v = U(obj.indV); w = U(~obj.indV);
                    
                    R = [e_dx2*(U(obj.Uind) - 2*U(obj.Uind+1) + U(obj.Uind+2)) + einv*(FHN.f(U(obj.Uind+1))-U(obj.Uind+3));...
                        obj.b*v(obj.numVvar+1:end) - obj.gamma*w(obj.indW)] + obj.G + obj.B*feval(obj.config.inFunc,t);
                    
                    J = [reshape([e_dx2*ones(obj.numVeq,1),-2*e_dx2+einv*FHN.df(U(obj.Uind+1)),e_dx2*ones(obj.numVeq,1),-einv*ones(obj.numVeq,1)]',4*obj.numVeq,1);repmat([obj.b;-obj.gamma],obj.numWeq,1)];
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
            %obj - OneDBurgers object
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
            %in the full state vector
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - FHN object
            %nodenum - scalar specifying node numbers
            %
            %Outputs:
            %--------
            %ind     - vector specifying the corresponding indices in the
            %          full state vector
            %--------------------------------------------------------------
 
            ind = [nodenum; nodenum + obj.config.ndof/2];
        end
        
        function   [JacStruct] = JacobianStructure(obj)
            %This function forms a boolean matrix whose true entries
            %correspond to nonzero entries in the jacobian
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - FHN object
            %
            %Outputs:
            %--------
            %JacStruct - ndof x ndof boolean matrix whose true values
            %            correspond to nonzero entries in the jacobian 
            %--------------------------------------------------------------
            
            N = obj.config.ndof/2;
            
%             A1 = diag(ones(N,1),0) + diag([ones(N-2,1);0],-1) + diag([0;ones(N-2,1)],1);
            A = spdiags([[ones(N-2,1);0;0],ones(N,1),[0;0;ones(N-2,1)]],-1:1,N,N);
            A(1,3) = 1; A(end,end-2) = 1;
            I = speye(N,N);

            JacStruct = [A,I;I,I];           
        end 
        
        function  [probGNAT] = createCopy4GNAT(obj,gnat)
            %This function makes a copy of the current object and only
            %copies the source term and input matrix evaluated at indN
            %rows.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - FHN object
            %gnat    - GNAT object
            %
            %Outputs:
            %--------
            %probGNAT - new "reduced" FHN object
            %--------------------------------------------------------------
            
            probGNAT = FHN([],obj);
            
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
        
            if length(obj.B) > 1
                probGNAT.B = obj.B(gnat.sampleInd,:);
            end
            
            if length(obj.G) > 1
                probGNAT.G = obj.G(gnat.sampleInd,1);
            end
            
            if length(obj.ic) > 1
                probGNAT.ic = obj.ic(unique(gnat.jrow),1);
            end
            
            indLastV = find(gnat.sampleInd <= obj.config.ndof/2,1,'last');
            probGNAT.indV = (gnat.jrow <= obj.config.ndof/2);
            if gnat.sampleInd(1) == 1 && (gnat.sampleInd(indLastV) == obj.config.ndof/2)
                probGNAT.caseNum = 1;
                probGNAT.Uind = gnat.irstart(2:indLastV-1);
                probGNAT.numVeq = length(2:indLastV-1); %m
            elseif gnat.sampleInd(1) == 1
                probGNAT.caseNum = 2;
                probGNAT.Uind = gnat.irstart(2:indLastV);
                probGNAT.numVeq = length(2:indLastV); %m
                
%                 probGNAT.Uind = gnat.irstart(1:indLastV-1);
%                 probGNAT.numVeq = length(1:indLastV-1); %m
            elseif gnat.sampleInd(end) == obj.config.ndof
                probGNAT.caseNum = 3;
                probGNAT.Uind = gnat.irstart(1:indLastV-1);
                probGNAT.numVeq = length(1:indLastV-1); %m
                
%                 probGNAT.Uind = gnat.irstart(2:indLastV);
%                 probGNAT.numVeq = length(2:indLastV); %m
            else
                probGNAT.caseNum = 4;
                probGNAT.Uind = gnat.irstart(1:indLastV);
                probGNAT.numVeq = length(1:indLastV); %m
            end
            probGNAT.numWeq = length(gnat.sampleInd(indLastV+1:end));
            probGNAT.numVvar = sum(probGNAT.indV) - probGNAT.numWeq; %n
            probGNAT.indW = indLastV+1:sum(~probGNAT.indV);
            
        end
    
        function  [probGNAT] = createCopy4GNATloc(obj,gnat)
            %This function makes a copy of the current object and only
            %copies the source term and input matrix evaluated at sampleInd
            %rows.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - FHN object
            %gnat    - locGNAT object
            %
            %Outputs:
            %--------
            %probGNAT - new "reduced" FHN object
            %--------------------------------------------------------------
            
            for kk = 1:gnat.nBases
                probGNAT(kk) = FHN([],obj);
                
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
                
                if length(obj.B) > 1
                    probGNAT(kk).B = obj.B(gnat.sampleInd{kk,1},:);
                end
                
                if length(obj.G) > 1
                    probGNAT(kk).G = obj.G(gnat.sampleInd{kk,1},1);
                end
                
                if length(obj.ic) > 1
                    probGNAT(kk).ic = obj.ic(unique(gnat.jrow{kk,1}),1);
                end
                
                indLastV = find(gnat.sampleInd{kk,1} <= obj.config.ndof/2,1,'last');
                probGNAT(kk).indV = (gnat.jrow{kk,1} <= obj.config.ndof/2);
                if gnat.sampleInd{kk,1}(1) == 1 && (gnat.sampleInd{kk,1}(indLastV) == obj.config.ndof/2)
                    probGNAT(kk).caseNum = 1;
                    probGNAT(kk).Uind = gnat.irstart{kk,1}(2:indLastV-1);
                    probGNAT(kk).numVeq = length(2:indLastV-1); %m
                elseif gnat.sampleInd{kk,1}(1) == 1
                    probGNAT(kk).caseNum = 2;
                    probGNAT(kk).Uind = gnat.irstart{kk,1}(1:indLastV-1);
                    probGNAT(kk).numVeq = length(1:indLastV-1); %m
                elseif gnat.sampleInd{kk,1}(end) == obj.config.ndof
                    probGNAT(kk).caseNum = 3;
                    probGNAT(kk).Uind = gnat.irstart{kk,1}(2:indLastV);
                    probGNAT(kk).numVeq = length(2:indLastV); %m
                else
                    probGNAT(kk).caseNum = 4;
                    probGNAT(kk).Uind = gnat.irstart{kk,1}(1:indLastV);
                    probGNAT(kk).numVeq = length(1:indLastV); %m
                end
                probGNAT(kk).numWeq = length(gnat.sampleInd{kk,1}(indLastV+1:end));
                probGNAT(kk).numVvar = sum(probGNAT(kk).indV) - probGNAT(kk).numWeq; %n
                probGNAT(kk).indW = indLastV+1:sum(~probGNAT(kk).indV);
            end
        end   
        
        function   [axOUT] = problemPaperPlot(obj,modelobj,modelAuxobj,pstr,axIN)
            %This function generates the FHN plot from Sorensen paper.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj         - Problem object
            %modelobj    - FOM, ROM, GNAT, or TPWL object
            %modelAuxobj - ROM object if modelobj is of class GNAT or
            %              locGNAT.
            %pstr        - cell array containing plot information.
            %axIN        - axis handle to place new plots
            %
            %Outputs:
            %--------
            %axOUT    - axis handle that new plots were placed
            %--------------------------------------------------------------
            
            %Extract reconstructed state vectors 
            [~,~,v] = obj.returnPosSV(modelobj,modelAuxobj,1);
            [~,~,w] = obj.returnPosSV(modelobj,modelAuxobj,2);
            
            freq = 5;
            
            %Set up axes
            if isempty(axIN)
                ax = axes;
                title('Phase Space Diagram','fontsize',16,'interpreter','latex');
                xlabel('$v(x,t)$','fontsize',14,'interpreter','latex');
                ylabel('$w(x,t)$','fontsize',14,'interpreter','latex');
            else
                ax = axIN;
            end
            
            %Get the current next plot tag
            nptag = get(ax,'nextplot');
            
            for i = 2:freq:modelobj.time.nstep+1 %don't plot IC b/c it is boring (v=w=0)
                if i == 1
                    %on the first iteration, make sure that I'm not
                    %overwritting data that I just plotted
                    set(ax,'nextplot','add');
                end
                plot(ax,v(:,i),w(:,i),pstr{:});
            end
            
            %Reset nextplot tag to what it was orginially
            set(ax,'nextplot',nptag);
            axOUT = ax;
        end
        
        function   [X,Y,svDBC] = returnPosSV(obj,modelobj,modelAuxobj,svid,flag,basenum)
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
            %svid        - scalar indicating the component of the state
            %              vector to extract (1 -> v, 2 -> w).
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
            
            v = svF(1:obj.config.ndof/2,:); w = svF(obj.config.ndof/2+1:end,:);
            
            switch svid
                case 1
                    %Extract mesh spacing
                    X = obj.mesh.node;
                    Y = [];
                    
                    %Include DBC in state
                    svDBC = v;
                case 2
                    %Extract mesh spacing
                    X = obj.mesh.node;
                    Y = [];
                    
                    %Include DBC in state
                    svDBC = w;
                otherwise
                    error('svid must be 1 or 2 for the FHN problem');
            end
        end
        
        function   [] = problemAnimate(obj,modelobj,modelAuxobj,pstr,type,freq)
            %This function generates a FHN animation.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj         - Problem object
            %modelobj    - FOM, ROM, GNAT, or TPWL object
            %modelAuxobj - ROM object if modelobj is of class GNAT
            %pstr        - cell array containing the plotting options
            %type        - scalar indicating which animation to use
            %freq        - scalar indicating how often to plot sv in
            %              animation
            %
            %Outputs:
            %--------
            %There are no outputs for this function.  It generates an
            %animiation.
            %--------------------------------------------------------------
            
            if nargin < 6
                freq = 1;
            end
            
            %Extract mesh spacing
            x = obj.mesh.node;
            
            %Establish time vector
            tvec = modelobj.time.T(1):modelobj.time.dt:modelobj.time.T(2);
            
            set(gcf,'position',[439   566   714   540]);
            switch type
                case 1 %v
                    ax = axes;
                    [~,~,v] = obj.returnPosSV(modelobj,modelAuxobj,1);
                    
                    %Determine plot limits
                    xmin = min(x); xmax = max(x);
                    ymin = min(min(v)); ymax = max(max(v));
                                        
                    title(ax,'Voltage as a Function of Time','fontsize',16,'fontname','georgia');
                    xlabel(ax,'Position','fontsize',14,'fontname','georgia');
                    ylabel(ax,'Voltage, v(x)','fontsize',14,'fontname','georgia');
                    set(ax,'nextplot','replacechildren','ylim',[ymin,ymax],'xlim',[xmin,xmax]);
                    
                    for i = 1:freq:length(tvec)
                        plot(ax,x,v(:,i),pstr{:}); pause(0.05);
                    end
                case 2 %w
                    ax = axes;
                    [~,~,w] = obj.returnPosSV(modelobj,modelAuxobj,2);
                    
                    %Determine plot limits
                    xmin = min(x); xmax = max(x);
                    ymin = min(min(w)); ymax = max(max(w));
                                        
                    title(ax,'Recovery of Voltage as a Function of Time','fontsize',16,'fontname','georgia');
                    xlabel(ax,'Position','fontsize',14,'fontname','georgia');
                    ylabel(ax,'Recovery of Voltage, w(x)','fontsize',14,'fontname','georgia');
                    set(ax,'nextplot','replacechildren','ylim',[ymin,ymax],'xlim',[xmin,xmax]);
                    
                    for i = 1:freq:length(tvec)
                        plot(ax,x,w(:,i),pstr{:}); pause(0.05);
                    end
                case 3 %v and w
                    [~,~,v] = obj.returnPosSV(modelobj,modelAuxobj,1);
                    [~,~,w] = obj.returnPosSV(modelobj,modelAuxobj,2);
                    
                    %Determine plot limits
                    xmin1 = min(x); xmax1 = max(x);
                    ymin1 = min(min(v)); ymax1 = max(max(v));
                    
                    xmin2 = min(x); xmax2 = max(x);
                    ymin2 = min(min(w)); ymax2 = max(max(w));
                    
                    ax(1) = subplot(2,1,1);
                    ax(2) = subplot(2,1,2);
                    
                    title(ax,'Voltage as a Function of Time','fontsize',16,'fontname','georgia');
                    xlabel(ax,'Position','fontsize',14,'fontname','georgia');
                    ylabel(ax,'Voltage, v(x)','fontsize',14,'fontname','georgia');
                    set(ax,'nextplot','replacechildren','ylim',[ymin,ymax],'xlim',[xmin,xmax]);
                    
                    title(ax,'Recovery of Voltage as a Function of Time','fontsize',16,'fontname','georgia');
                    xlabel(ax,'Position','fontsize',14,'fontname','georgia');
                    ylabel(ax,'Recovery of Voltage, w(x)','fontsize',14,'fontname','georgia');
                    set(ax,'nextplot','replacechildren','ylim',[ymin,ymax],'xlim',[xmin,xmax]);

                    for i = 1:freq:length(tvec)
                        plot(ax,x,v(:,i),pstr{:});
                        plot(ax,x,w(:,i),pstr{:}); pause(0.05);
                    end
            end
        end
        
        function  [varargout] = setPrevSV(varargin)
        end
    end
    
    methods (Static)
        function [nl] = f(x)
            %This function evaluates the nonlinear FHN problem at x
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %x  - point to evalute the nonlinear function
            %
            %Outputs:
            %--------
            %nl - f(x)
            %--------------------------------------------------------------
            
            nl = x.*(x-0.1).*(1-x);
        end
        
        function  [nlder] = df(x)
            %This function evaluates the nonlinear FHN problem at x
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %x     - point to evalute the nonlinear function derivative
            %
            %Outputs:
            %--------
            %nlder - df(x)
            %--------------------------------------------------------------
            
            nlder = -3*(x.*x) + 2.2*x - 0.1;
        end
    end
end
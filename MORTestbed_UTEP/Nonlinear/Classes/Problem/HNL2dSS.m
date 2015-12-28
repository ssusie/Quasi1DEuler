classdef HNL2dSS
    
    properties (GetAccess=public, SetAccess=private)
        config;
        staticFlag = true;
        
        mesh;
        ndim;
        Jstruct;
        ic;
        B;
        G;
        C;
        D;
        dx;
        
        M;
        N;
        coord;
    end
    
    properties (Hidden=true, GetAccess=public, SetAccess=private)
        Jconst;
        jdiag;
        H1;
        H2;
        H3;
        H4;
        ind_i_j;
        ind_ip1_j;
        ind_im1_j;
        ind_i_jp1;
        ind_i_jm1;
        ResInd_ip1_j;
        ResInd_im1_j;
    end
    
    methods
        function  [obj] = HNL2dSS(cfgobj,oldobj)
            %This is the constructor for the Highly Nonlinear 2d Steady
            %State problem. If there is 1 input, the constructor will
            %create the class instance from the variables defined in this
            %functino.  If there are 2 inputs, the first input will be
            %unused and the output will be a deep copy of the oldobj.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %cfgobj - CONFIG object
            %oldobj - HNL2dSS object whose properties will be copied
            %         into the new object.
            %
            %Outputs:
            %--------
            %obj    - HNL2dSS object
            %--------------------------------------------------------------
            
            if nargin == 0
                %Keep everything at defaults if no input arguments
                return;
            end
            
            % Construct a new object based on a deep copy of an old object
            % of this class by copying properties over.
            if nargin == 2 && isa(oldobj,'HNL2dSS')
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
            
            %Establish the order of the problem
            obj.M = cfgobj.nNodes(1) - 2;
            obj.N = cfgobj.nNodes(2) - 2;
            
            %Determine grid spacing
            obj.dx = (cfgobj.DLim(:,2) - cfgobj.DLim(:,1))./(cfgobj.nNodes(:)-1);
            
            %Set up grid information
            obj.coord.x = [cfgobj.DLim(1,1):obj.dx(1):cfgobj.DLim(1,2)]';
            obj.coord.y = [cfgobj.DLim(2,1):obj.dx(2):cfgobj.DLim(2,2)]';
            x_svform = repmat(obj.coord.x(2:end-1,1),obj.N,1);
            y_svform = reshape(repmat(obj.coord.y(2:end-1,1)',obj.M,1),obj.M*obj.N,1);
            
            xFull = repmat(obj.coord.x,obj.N+2,1);
            yFull = reshape(repmat(obj.coord.y',obj.M+2,1),(obj.M+2)*(obj.N+2),1);
            
            %Determine mesh
            obj.mesh.node = [xFull, yFull];
            obj.mesh.conn = [1:cfgobj.nNodes(1)-1;2:cfgobj.nNodes(1)]';
            obj.mesh.conn = [];
            for j = 1:cfgobj.nNodes(2)-1
                obj.mesh.conn = [obj.mesh.conn;...
                    [((j-1)*cfgobj.nNodes(1) + 1):j*cfgobj.nNodes(1)-1]',...
                    [((j-1)*cfgobj.nNodes(1) + 2):j*cfgobj.nNodes(1)]',...
                    [(j*cfgobj.nNodes(1) + 1):(j+1)*cfgobj.nNodes(1)-1]',...
                    [(j*cfgobj.nNodes(1) + 2):(j+1)*cfgobj.nNodes(1)]'];
            end
            obj.mesh.conn = [[1:size(obj.mesh.conn,1)]',obj.mesh.conn];
            temp = find(xFull == cfgobj.DLim(1,1) | xFull == cfgobj.DLim(1,2) | yFull == cfgobj.DLim(2,1) | yFull == cfgobj.DLim(2,2));
            obj.mesh.dbc = [temp,true(length(temp),1)];
            
            %Number of unknowns per node
            obj.ndim = 1;
            
            %Setup Jacobian Structure
            obj.Jstruct = JacobianStructure(obj);
            
            %Determine source term vector from the source function and the
            %coordinate vector
            if isa(cfgobj.GFunc,'function_handle')
                obj.G = cfgobj.GFunc(x_svform,y_svform);
            else
                obj.G = feval(cfgobj.GFunc,x_svform,y_svform);
            end
            
            %Determine ic term vector from the ic function and the
            %coordinate vector
            if isa(cfgobj.icFunc,'function_handle')
                obj.ic = cfgobj.icFunc(x_svform,y_svform);
            else
                obj.ic = feval(cfgobj.icFunc,x_svform,y_svform);
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
        
        function  [R,J] =  ResJac(obj,U,~)
            %This function calcuates the residual and jacobian of Burger's Equation
            %--------------------------------------------------------------------------
            %Inputs:
            %-------
            %obj - HNL2dSS object
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
            
            U_i_jp1 = [U(obj.M+1:obj.M*obj.N,1);obj.h3(obj.coord.x(2:end-1))];
            U_i_jm1 = [obj.h1(obj.coord.x(2:end-1));U(1:obj.M*obj.N-obj.M)];
            
            U_ip1_j = zeros(obj.config.ndof,1);
            U_im1_j = zeros(obj.config.ndof,1);
            ind_ip1_j = repmat([true(obj.M-1,1);false],obj.N,1);
            ind_im1_j = repmat([false;true(obj.M-1,1)],obj.N,1);
            U_ip1_j(ind_ip1_j,1)  = U(ind_im1_j,1);
            U_ip1_j(~ind_ip1_j,1) = obj.h2(obj.coord.y(2:end-1));
            U_im1_j(ind_im1_j,1)  = U(ind_ip1_j,1);
            U_im1_j(~ind_im1_j,1) = obj.h4(obj.coord.y(2:end-1));
            
            R = (U_ip1_j + U_im1_j)/(obj.dx(1)^2) + (U_i_jp1 + U_i_jm1)/(obj.dx(2)^2) - (2/(obj.dx(1)^2+obj.dx(2)^2))*U + ...
                (obj.config.param(1)/obj.config.param(2))*(exp(obj.config.param(2)*U)-1) + obj.G;
            J = spdiags([[ones(obj.M*obj.N-obj.M,1)/(obj.dx(2)^2);zeros(obj.M,1)],[ind_im1_j(2:end,1)/(obj.dx(1)^2);0],...
                (-2/(obj.dx(1)^2+obj.dx(2)^2))*ones(obj.M*obj.N,1) + obj.config.param(1)*exp(obj.config.param(2)*U),...
                [0;ind_ip1_j(1:end-1,1)/(obj.dx(1)^2)],[zeros(obj.M,1);ones(obj.M*obj.N-obj.M,1)/(obj.dx(2)^2)]],...
                [-obj.M,-1:1,obj.M],obj.M*obj.N,obj.M*obj.N);
%             J1 = diag((-2/(obj.dx(1)^2+obj.dx(2)^2))*ones(obj.M*obj.N,1) + obj.config.param(1)*exp(obj.config.param(2)*U),0) + ...
%                 diag(ones(obj.M*obj.N-obj.M,1)/(obj.dx(2)^2),obj.M) + ...
%                 diag(ones(obj.M*obj.N-obj.M,1)/(obj.dx(2)^2),-obj.M) + ...
%                 diag(ind_ip1_j(1:end-1,1)/(obj.dx(1)^2),1) + diag(ind_im1_j(2:end,1)/(obj.dx(1)^2),-1);
        end
        
        function  [R,J] =  ResJacGNAT(obj,U,~)
            %This function calcuates the residual and jacobian of HNL2dSS
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - HNL2dSS object
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
            
            U_i_j = U(obj.ind_i_j);
            U_i_jp1 = [U(obj.ind_i_jp1);obj.H3];
            U_i_jm1 = [obj.H1; U(obj.ind_i_jm1)];
            
            U_ip1_j(obj.ResInd_ip1_j,1) = U(obj.ind_ip1_j);
            U_ip1_j(~obj.ResInd_ip1_j,1) = obj.H2;
            U_im1_j(obj.ResInd_im1_j,1) = U(obj.ind_im1_j);
            U_im1_j(~obj.ResInd_im1_j,1) = obj.H4;
            
            R = (1/(obj.dx(1)*obj.dx(1)))*(U_ip1_j + U_im1_j) + ...
                (1/(obj.dx(2)*obj.dx(2)))*(U_i_jp1 + U_i_jm1) - ...
                (2/(obj.dx(1)^2 + obj.dx(2)^2))*U_i_j + ...
                (obj.config.param(1)/obj.config.param(2))*(exp(obj.config.param(2)*U_i_j)-1) + obj.G;
            
            J = obj.Jconst;
            J(logical(obj.jdiag),1) = J(logical(obj.jdiag)) + obj.config.param(1)*exp(obj.config.param(2)*U_i_j);
            
        end  
                
        function  [ind] = node2ind(obj,nodenum)
            %This function takes a node in the mesh and maps it to an index
            %in the full state vector
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - HNL2dSS object
            %nodenum - scalar specifying node numbers
            %
            %Outputs:
            %--------
            %ind     - vector specifying the corresponding indices in the
            %          full state vector
            %--------------------------------------------------------------
            
            ind = zeros(length(nodenum),1);
            for i = 1:length(nodenum)
                if ~isempty(find(obj.mesh.dbc(:,1) == nodenum(i),1))
                    ind(i) = nan;
                else
                    temp = floor(nodenum(i)/obj.config.nNodes(1));
                    ind(i) = nodenum(i) - obj.config.nNodes(1) - 2*(temp-1) - 1;
                end
            end
        end
    
        function  [JacStruct] = JacobianStructure(obj)
            %This function forms a boolean matrix whose true entries
            %correspond to nonzero entries in the jacobian
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - HNL2dSS object
            %
            %Outputs:
            %--------
            %JacStruct - ndof x ndof boolean matrix whose true values
            %            correspond to nonzero entries in the jacobian 
            %--------------------------------------------------------------
            
            m = obj.M;
            n = obj.N;
            
            vec = repmat([ones(m-1,1);0],n,1); vec(end,:) = [];
            JacStruct = spdiags([[ones(m*n-m,1);zeros(m,1)],[vec;0],ones(m*n,1),[0;vec],[zeros(m,1);ones(m*n-m,1)]],[-m,-1,0,1,m],m*n,m*n);
%             temp = diag(ones(m*n,1),0) + diag(ones(m*n-m,1),m) + diag(ones(m*n-m,1),-m) + ...
%                 diag(vec,-1) + diag(vec,1);
%             JacStruct = sparse(temp);
        end
        
        function  [probGNAT] = createCopy4GNAT(obj,gnat)
            %This function makes a copy of the current object and only
            %copies the source term and input matrix evaluated at indN
            %rows.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - HNL2dSS object
            %indN    - scalar specifying node numbers
            %
            %Outputs:
            %--------
            %probGNAT - new "reduced" HNL2dSS object
            %--------------------------------------------------------------
            
            probGNAT = HNL2dSS([],obj);
            
            temp = probGNAT.config.param;
            probGNAT.config = [];
            probGNAT.config.param = temp;
            probGNAT.mesh   = [];
            probGNAT.ndim = [];
            probGNAT.Jstruct = [];
        
            if length(obj.B) > 1
                probGNAT.B = obj.B(gnat.sampleInd,1);
            end
            
            if length(obj.G) > 1
                probGNAT.G = obj.G(gnat.sampleInd,1);
            end
            
            if length(obj.ic) > 1
                probGNAT.ic = obj.ic(unique(gnat.jrow),1);
            end
            
            %Store the jdiag GNAT property in this object
            probGNAT.jdiag = gnat.jdiag;
            
            %Define the constant part of the jacobian
            [~,tempJ] = obj.ResJac(-inf*ones(obj.config.ndof,1),0);
            tempJ = tempJ(gnat.sampleInd,:);
            probGNAT.Jconst = [];
            for i = 1:size(tempJ,1)
                for j = 1:size(tempJ,2)
                    if abs(tempJ(i,j)) > 1e-16
                        probGNAT.Jconst = [probGNAT.Jconst; tempJ(i,j)];
                    end
                end
            end
            
            tempx = obj.coord.x(2:end-1);
            tempy = obj.coord.y(2:end-1);
            probGNAT.ind_i_j = false(length(gnat.jrow),1);
            probGNAT.ind_im1_j = false(length(gnat.jrow),1);
            probGNAT.ind_ip1_j = false(length(gnat.jrow),1);
            probGNAT.ind_i_jm1 = false(length(gnat.jrow),1);
            probGNAT.ind_i_jp1 = false(length(gnat.jrow),1);
            probGNAT.ResInd_im1_j = true(length(gnat.sampleInd),1);
            probGNAT.ResInd_ip1_j = true(length(gnat.sampleInd),1);

            probGNAT.H1 = []; probGNAT.H2 = []; probGNAT.H3 = []; probGNAT.H4 = [];
            
            for i = 1:length(gnat.sampleInd)
                x = tempx(mod(gnat.sampleInd(i)-1,obj.M)+1);
                y = tempy(ceil(gnat.sampleInd(i)/obj.M));
               if gnat.sampleInd(i) <= obj.M
                   if gnat.sampleInd(i) == 1
                       probGNAT.ind_i_j(gnat.irstart(i),1) = true;
                       probGNAT.ind_ip1_j(gnat.irstart(i)+1,1) = true;
                       probGNAT.ind_i_jp1(gnat.irstart(i)+2,1) = true;
                       
                       probGNAT.ResInd_im1_j(i) = false;
                       
                       probGNAT.H1 = [probGNAT.H1; obj.h1(x)];
                       probGNAT.H4 = [probGNAT.H4; obj.h4(y)];
                   elseif gnat.sampleInd(i) == obj.M;
                       probGNAT.ind_im1_j(gnat.irstart(i),1) = true;
                       probGNAT.ind_i_j(gnat.irstart(i)+1,1) = true;
                       probGNAT.ind_i_jp1(gnat.irstart(i)+2,1) = true;
                       
                       probGNAT.ResInd_ip1_j(i) = false;
                       
                       probGNAT.H1 = [probGNAT.H1; obj.h1(x)];
                       probGNAT.H2 = [probGNAT.H2; obj.h2(y)];
                   else
                       probGNAT.ind_im1_j(gnat.irstart(i),1) = true;
                       probGNAT.ind_i_j(gnat.irstart(i)+1,1) = true;
                       probGNAT.ind_ip1_j(gnat.irstart(i)+2,1) = true;
                       probGNAT.ind_i_jp1(gnat.irstart(i)+3,1) = true;
                       
                       probGNAT.H1 = [probGNAT.H1; obj.h1(x)];
                   end
               elseif gnat.sampleInd(i) >= obj.M*obj.N - obj.M + 1 
                   if gnat.sampleInd(i) == obj.M*obj.N - obj.M + 1
                       probGNAT.ind_i_jm1(gnat.irstart(i),1) = true;
                       probGNAT.ind_i_j(gnat.irstart(i)+1,1) = true;
                       probGNAT.ind_ip1_j(gnat.irstart(i)+2,1) = true;
                       
                       probGNAT.ResInd_im1_j(i) = false;
                       
                       probGNAT.H3 = [probGNAT.H3; obj.h3(x)];
                       probGNAT.H4 = [probGNAT.H4; obj.h4(y)];
                   elseif gnat.sampleInd(i) == obj.M*obj.N
                       probGNAT.ind_i_jm1(gnat.irstart(i),1) = true;
                       probGNAT.ind_im1_j(gnat.irstart(i)+1,1) = true;
                       probGNAT.ind_i_j(gnat.irstart(i)+2,1) = true;
                       
                       probGNAT.ResInd_ip1_j(i) = false;
                       
                       probGNAT.H2 = [probGNAT.H2; obj.h2(y)];
                       probGNAT.H3 = [probGNAT.H3; obj.h3(x)];
                   else
                       probGNAT.ind_i_jm1(gnat.irstart(i),1) = true;
                       probGNAT.ind_im1_j(gnat.irstart(i)+1,1) = true;
                       probGNAT.ind_i_j(gnat.irstart(i)+2,1) = true;
                       probGNAT.ind_ip1_j(gnat.irstart(i)+3,1) = true;
                       
                       probGNAT.H3 = [probGNAT.H3; obj.h3(x)];
                   end
               elseif mod(gnat.sampleInd(i),obj.M) == 0
                   probGNAT.ind_i_jm1(gnat.irstart(i),1) = true;
                   probGNAT.ind_im1_j(gnat.irstart(i)+1,1) = true;
                   probGNAT.ind_i_j(gnat.irstart(i)+2,1) = true;
                   probGNAT.ind_i_jp1(gnat.irstart(i)+3,1) = true;
                   
                   probGNAT.ResInd_ip1_j(i) = false;
                   
                   probGNAT.H2 = [probGNAT.H2; obj.h2(y)];
               elseif mod(gnat.sampleInd(i)-1,obj.M) == 0
                   probGNAT.ind_i_jm1(gnat.irstart(i),1) = true;
                   probGNAT.ind_i_j(gnat.irstart(i)+1,1) = true;
                   probGNAT.ind_ip1_j(gnat.irstart(i)+2,1) = true;
                   probGNAT.ind_i_jp1(gnat.irstart(i)+3,1) = true;
                   
                   probGNAT.ResInd_im1_j(i) = false;
                   
                   probGNAT.H4 = [probGNAT.H4; obj.h4(y)];
               else
                   probGNAT.ind_i_jm1(gnat.irstart(i),1) = true;
                   probGNAT.ind_im1_j(gnat.irstart(i)+1,1) = true;
                   probGNAT.ind_i_j(gnat.irstart(i)+2,1) = true;
                   probGNAT.ind_ip1_j(gnat.irstart(i)+3,1) = true;
                   probGNAT.ind_i_jp1(gnat.irstart(i)+4,1) = true;
               end
            end
        end
        
        function   [axOUT] = problemPaperPlot(obj,modelobj,modelAuxobj,pstr,axIN)
            %This function generates the HNL2dSS plot from Sorensen paper.
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
            [X,Y,U] = obj.returnPosSV(modelobj,modelAuxobj);
            
            %Find the indicies corresponding to the line x = y
            ind = (abs(X - Y) < 1e-10);
            temp = X(ind); 
            temp = unique(temp(:));
            s = sqrt(2)*(temp - temp(1));
            %Extract the indices from U along the line x = y
            u = U(ind);
            
            %Set up axes
            if isempty(axIN)
                ax(1) = axes;
                title('HNL2dSS Solution over Entire Domain','fontsize',16,'interpreter','latex');
                xlabel('X','fontsize',14,'interpreter','latex');
                ylabel('Y','fontsize',14,'interpreter','latex');
                zlabel('Z (solution)','fontsize',14,'interpreter','latex');
            else
                ax(1) = axIN;
            end
            
            figure;
            ax(2) = axes;
            title('HNL2dSS Solution along line $x = y','fontsize',16,'interpreter','latex');
            xlabel('Distance Along Line $x = y$','fontsize',14,'interpreter','latex');
            ylabel('Solution','fontsize',14,'interpreter','latex');
                
            surf(ax(1),X,Y,U);
            plot(ax(2),s,u,pstr{:});
            
            %Reset nextplot tag to what it was orginially
            axOUT = ax;
        end
        
        function [X,Y,svDBC] = returnPosSV(obj,modelobj,modelAuxobj,~,flag,~)
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
            %              otherwise.
            %flag        - string ('pod','res','jac','sv')
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
                    svF = modelobj.phi;
                case 'res'
                    svF = modelobj.res;
                case 'jac'
                    svF = modelobj.jac;
                case 'sv'
                    svF = modelobj.reconstructFullState(modelAuxobj);
                    svF(:,1) = []; %Don't care about initial guess
            end
            n = size(svF,2);
            
            %Extract and format mesh data and set up svDBC matrix that will
            %contain the state vector with its DBCs
            [X,Y] = meshgrid(obj.coord.x,obj.coord.y);
            svDBC = zeros(length(obj.coord.x),length(obj.coord.y),n);
            
            %Fill svDBC
            temp = reshape(svF,obj.M,obj.N,n);
            svDBC(2:end-1,2:end-1,:) = temp;
            svDBC(1,:,:) = repmat(obj.h1(obj.coord.x)',[1,1,n]);
            svDBC(end,:,:) = repmat(obj.h3(obj.coord.x)',[1,1,n]);
            svDBC(:,1,:) = repmat(obj.h4(obj.coord.y)',[1,1,n]);
            svDBC(:,end,:) = repmat(obj.h2(obj.coord.y)',[1,1,n]);
        end
        
        function   [] = problemAnimate(obj,modelobj,modelAuxobj,pstr,~,~)
            %This function generates the animation for the HNL2dSS
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj      - Problem object
            %modelobj - FOM, ROM, GNAT, or TPWL object
            %modelAuxobj - ROM object if modelobj is of class GNAT or
            %              locGNAT (unused otherwise)
            %pstr        - cell array indicating plot object formatting
            %
            %Outputs:
            %--------
            %There are no outputs for this function.  It generates an
            %animiation.
            %--------------------------------------------------------------
            
            %This is a steady state problem so the animation is just the
            %solution plot
            obj.problemPaperPlot(modelobj,modelAuxobj,pstr,[]);
        end
    end
    
    methods (Static)
        function [u] = h1(x)
            %This function defines the DBC along the bottom edge of the
            %domain.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %x - vector containing x coordinates
            %
            %Outputs:
            %--------
            %u - vector containing the DBC at each position x
            %--------------------------------------------------------------
            
            u = zeros(length(x),1);
        end
        
        function [u] = h2(y)
            %This function defines the DBC along the right edge of the
            %domain.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %y - vector containing y coordinates
            %
            %Outputs:
            %--------
            %u - vector containing the DBC at each position y
            %--------------------------------------------------------------
            
            u = zeros(length(y),1);
        end
        
        function [u] = h3(x)
            %This function defines the DBC along the top edge of the
            %domain.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %x - vector containing x coordinates
            %
            %Outputs:
            %--------
            %u - vector containing the DBC at each position x
            %--------------------------------------------------------------
            
            u = zeros(length(x),1);
        end
        
        function [u] = h4(y)
            %This function defines the DBC along the left edge of the
            %domain.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %y - vector containing y coordinates
            %
            %Outputs:
            %--------
            %u - vector containing the DBC at each position y
            %--------------------------------------------------------------
            
            u = zeros(length(y),1);
        end
    end
end
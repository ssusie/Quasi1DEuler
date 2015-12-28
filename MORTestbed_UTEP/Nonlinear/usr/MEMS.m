classdef MEMS < handle
    
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
        dy;
        coord;
        
        S0;
        E;
        In0;
        rho0;
        lambda;
        p0;
        z0;
        width;
        mu;
        height;
        b0;
        del;
        INT;
        D2;
        grady;
        gradx2;
        gradx;
        LAP;
        D4;
        M;
        N;
        
        indLTn;
        indEQnp1;
        indEQnp2;
        indEQ2ns1;
        indEQ2n;
        indLT2nElse;
        indLT2npm;
        indGTmnp2nsmp1;
        indGTmnElse;
        
        ResInd;
        JacInd;
        
        lenRes;
        lenJac;
        
        p;
    end
    
    methods
        function  [obj] = MEMS(cfgobj,oldobj)
            %This is the constructor for the Micromachined Switch.  If
            %there is 1 input, the constructor will create the class
            %instance from the variables defined in this functino.  If
            %there are 2 inputs, the first input will be unused and the
            %output will be a deep copy of the oldobj.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %cfgobj - CONFIG object
            %oldobj - MEMS object whose properties will be copied
            %         into the new object.
            %
            %Outputs:
            %--------
            %obj - MEMS object
            %--------------------------------------------------------------
            
            if nargin == 0
                %Keep everything at defaults if no input arguments
                return;
            end
            
            % Construct a new object based on a deep copy of an old object
            % of this class by copying properties over.
            if nargin == 2 && isa(oldobj,'MEMS')
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
            
            obj.config = cfgobj;
            obj.N = cfgobj.nNodes(1)-2;
            obj.M = cfgobj.nNodes(2)-2;
            
            if isempty(obj.config.altFile)
                error('For the MEMS problem, you must specify inputs in an alternate file due to the large number of parameters');
            end
            
            [~,consts,operators] = feval(obj.config.altFile,[obj.M obj.N],obj.config.DLim);
            obj.S0 = consts(1);
            obj.E = consts(2);
            obj.In0 = consts(3);
            obj.rho0 = consts(4);
            obj.lambda = consts(5);
            obj.p0 = consts(6);
            obj.z0 = consts(7);
            obj.width = consts(8);
            obj.dx = consts(11);
            obj.dy = consts(9);
            obj.mu = consts(10);
            obj.height = consts(12);
            obj.b0 = consts(13);
            obj.del = consts(14);
            obj.INT = operators{1};
            obj.D2 = operators{3};
            obj.grady = operators{4};
            obj.B = operators{5};
            obj.C = operators{6};
            obj.gradx2 = operators{7};
            obj.gradx = operators{8};
            obj.LAP = operators{9};
            obj.D4 = operators{10};
            
            obj.ic = operators{11};
            
            %Set up grid information
            obj.coord.x = [cfgobj.DLim(1,1):obj.dx:cfgobj.DLim(1,2)]';
            obj.coord.y = [cfgobj.DLim(2,1):obj.dy:cfgobj.DLim(2,2)]';
            
            xFull = reshape(repmat(obj.coord.x',obj.M+2,1),(obj.N+2)*(obj.M+2),1);
            yFull = repmat(obj.coord.y,obj.N+2,1);
            
            %Determine mesh
            obj.mesh.node = [obj.coord.x, 0.5*(obj.coord.y(1)+obj.coord.y(2))*ones(length(obj.coord.x),1);...
                xFull      , yFull];
            
            %             obj.mesh.conn = [];
            obj.mesh.conn = [(1:length(obj.coord.x)-1)',(2:length(obj.coord.x))',zeros(length(obj.coord.x)-1,1),zeros(length(obj.coord.x)-1,1)];
            k = length(obj.mesh.conn);
            for j = 1:cfgobj.nNodes(1)-1
                obj.mesh.conn = [obj.mesh.conn;...
                    (k+1)+[((j-1)*cfgobj.nNodes(2) + 1):j*cfgobj.nNodes(2)-1]',...
                    (k+1)+[((j-1)*cfgobj.nNodes(2) + 2):j*cfgobj.nNodes(2)]',...
                    (k+1)+[(j*cfgobj.nNodes(2) + 1):(j+1)*cfgobj.nNodes(2)-1]',...
                    (k+1)+[(j*cfgobj.nNodes(2) + 2):(j+1)*cfgobj.nNodes(2)]'];
            end
            %             obj.mesh.conn = [obj.mesh.conn;...
            %                              (length(xFull)+1:length(xFull)+length(obj.coord.x)-1)', (length(xFull)+2:length(xFull)+length(obj.coord.x))',zeros(length((length(xFull)+1:length(xFull)+length(obj.coord.x)-1)'),1),zeros(length((length(xFull)+1:length(xFull)+length(obj.coord.x)-1)'),1)];
            obj.mesh.conn = [[1:size(obj.mesh.conn,1)]',obj.mesh.conn];
            temp = find(xFull == cfgobj.DLim(1,1) | xFull == cfgobj.DLim(1,2) | yFull == cfgobj.DLim(2,1) | yFull == cfgobj.DLim(2,2));
            temp2 = find(obj.coord.x == cfgobj.DLim(1,1) | obj.coord.x == cfgobj.DLim(1,2));
            dbc1 = [temp+length(obj.coord.x),true(length(temp),1),false(length(temp),1)];
            dbc2 = [temp2,true(length(temp2),2)];
            obj.mesh.dbc = [dbc2;dbc1];
            %Number of unknowns per node
            obj.ndim = (2*obj.N + obj.M*obj.N)/(obj.N + obj.M);
            
            %Setup Jacobian Structure
            obj.Jstruct = JacobianStructure(obj);
        end
        
        function  [R,J] =  ResJac(obj,U,t)
            %This function calcuates the residual and jacobian of the MEMS
            %switch
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
            %------------------------------------------------------------------
            
            x1 = U(1:obj.N); x2 = U(obj.N+1:2*obj.N); x3 = U(2*obj.N+1:obj.N*(obj.M+2));
            rho = obj.rho0*obj.height*obj.width;
            In = obj.In0*obj.width*obj.height^3;
            
            f1 = x2./(3*x1.^2)*1/obj.del;
            f2 = (2*x2.^2)./(3*x1.^3)*1/obj.del + (3*x1.^2)/(rho).*( obj.INT*(x3-obj.p0)+obj.S0*obj.height*obj.width*(obj.D2*(x1-obj.z0))-obj.E*In*(obj.D4*(x1-obj.z0)))*obj.del;
            
            arg31 = (x2./(3*x1.^3))*1/obj.del;
            f31 = -DiagHat((arg31),obj.M,obj.N) * x3;
            arg32 = (x1+4*obj.lambda)./(4*obj.mu).*(obj.gradx2*(x1-obj.z0));
            f32 = DiagHat((arg32),obj.M,obj.N) * (x3.*(obj.gradx*(x3-obj.p0)));
            arg33 = ((x1.^2+6*obj.lambda*x1)/(12*obj.mu));
            f33 = DiagHat((arg33),obj.M,obj.N) * (x3.*(obj.LAP*(x3-obj.p0))+(obj.gradx*(x3-obj.p0)).^2+(obj.grady*(x3-obj.p0)).^2);
            f3 = f31 + f32 + f33;
            
            R = [f1 ; f2 ; f3] + obj.B*feval(obj.config.inFunc,t);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            J11 = -2/3*spdiags(x2./(x1.^3),0,obj.N,obj.N)*1/obj.del;
            J12 = 1/3*spdiags(1./(x1.^2),0,obj.N,obj.N)*1/obj.del;
            J13 = spalloc(obj.N,obj.N*obj.M,0);
%             J11 = -2/3*diag(x2./(x1.^3))*1/obj.del; J11=sparse(J11);
%             J12 = 1/3*diag(1./(x1.^2))*1/obj.del;   J12=sparse(J12);
%             J13 = zeros(obj.N,obj.N*obj.M);         J13=sparse(J13);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            J21 = -2*spdiags(x2.^2./(x1.^4),0,obj.N,obj.N)*1/obj.del + 6/rho*spdiags(x1,0,obj.N,obj.N)*spdiags( obj.INT*(x3-obj.p0) + obj.S0*obj.height*obj.width*(obj.D2*(x1-obj.z0)) ...
                - obj.E*In*(obj.D4*(x1-obj.z0)), 0, obj.N,obj.N )*obj.del + 3/rho*spdiags(x1.^2,0,obj.N,obj.N)*(obj.S0*obj.height*obj.width*obj.D2 - obj.E*In*obj.D4)*obj.del;
            J22 = 4/3*spdiags((x2/obj.del)./(x1.^3),0,obj.N,obj.N);
            J23 = 3/rho*spdiags(x1.^2,0,obj.N,obj.N)*obj.INT*obj.del;
%             J21 = -2*diag(x2.^2./(x1.^4))*1/obj.del + 6/rho*diag(x1)*diag( obj.INT*(x3-obj.p0) + obj.S0*obj.height*obj.width*(obj.D2*(x1-obj.z0)) ...
%                 - obj.E*In*(obj.D4*(x1-obj.z0)) )*obj.del + 3/rho*diag(x1.^2)*(obj.S0*obj.height*obj.width*obj.D2 - obj.E*In*obj.D4)*obj.del;
%             J22 = 4/3*diag((x2/obj.del)./(x1.^3));
%             J23 = 3/rho*diag(x1.^2)*obj.INT*obj.del;
%             J21=sparse(J21); J22=sparse(J22); J23=sparse(J23);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            J31a = spdiags(x3,0,obj.N*obj.M,obj.N*obj.M)*WeirdDiag( spdiags((x2/obj.del)./(x1.^4),0,obj.N,obj.N),obj.M,obj.N) ;
            J31b = spdiags(x3.*(obj.gradx*(x3-obj.p0)),0,obj.N*obj.M,obj.N*obj.M)*WeirdDiag((spdiags(obj.gradx2*(x1-obj.z0),0,obj.N,obj.N)+spdiags(x1,0,obj.N,obj.N)*obj.gradx2+4*obj.lambda*obj.gradx2)/(4*obj.mu),obj.M,obj.N);
%             J31b = diag(x3.*(obj.gradx*(x3-obj.p0)))*WeirdDiag((diag(obj.gradx2*x1-obj.z0)+diag(x1)*obj.gradx2+4*obj.lambda*obj.gradx2)/(4*obj.mu),obj.M,obj.N);
            J31c = spdiags(x3.*(obj.LAP*(x3-obj.p0))+(obj.gradx*(x3-obj.p0)).^2+(obj.grady*(x3-obj.p0)).^2,0,obj.N*obj.M,obj.N*obj.M)*WeirdDiag( spdiags((2*x1+6*obj.lambda)/(12*obj.mu),0,obj.N,obj.N),obj.M,obj.N);
            J31 = J31a + J31b + J31c;
            %%%
            J32 = -spdiags(x3,0,obj.N*obj.M,obj.N*obj.M)*WeirdDiag( spdiags(1./(3*x1.^3),0,obj.N,obj.N),obj.M,obj.N)*1/obj.del;
            %%%
            arg331b = (x1+4*obj.lambda)./(4*obj.mu).*(obj.gradx2*(x1-obj.z0));
            J331 = -DiagHat(((x2/obj.del)./(3*x1.^3)),obj.M,obj.N) + DiagHat(arg331b,obj.M,obj.N)*(spdiags(obj.gradx*(x3-obj.p0),0,obj.N*obj.M,obj.N*obj.M) + spdiags(x3,0,obj.N*obj.M,obj.N*obj.M)*obj.gradx);
%             J331 = -DiagHat(((x2/obj.del)./(3*x1.^3)),obj.M,obj.N) + DiagHat(arg331b,obj.M,obj.N)*(diag(obj.gradx*(x3-obj.p0))) + diag(x3)*obj.gradx;
            %%%
            J332 = DiagHat(((x1.^2+6*obj.lambda*x1)/(12*obj.mu)),obj.M,obj.N)*(spdiags(obj.LAP*(x3-obj.p0),0,obj.N*obj.M,obj.N*obj.M) + spdiags(x3,0,obj.N*obj.M,obj.N*obj.M)*obj.LAP + 2*spdiags(obj.gradx*(x3-obj.p0),0,obj.N*obj.M,obj.N*obj.M)*obj.gradx+2*spdiags(obj.grady*(x3-obj.p0),0,obj.N*obj.M,obj.N*obj.M)*obj.grady);
            %%%
            J33 = J331 + J332;
%             J31a = diag(x3)*WeirdDiag( diag((x2/obj.del)./(x1.^4)),obj.M,obj.N) ;
%             J31b = diag(x3.*(obj.gradx*(x3-obj.p0)))*WeirdDiag((diag(obj.gradx2*(x1-obj.z0))+diag(x1)*obj.gradx2+4*obj.lambda*obj.gradx2)/(4*obj.mu),obj.M,obj.N);
% %             J31b = diag(x3.*(obj.gradx*(x3-obj.p0)))*WeirdDiag((diag(obj.gradx2*x1-obj.z0)+diag(x1)*obj.gradx2+4*obj.lambda*obj.gradx2)/(4*obj.mu),obj.M,obj.N);
%             J31c = diag(x3.*(obj.LAP*(x3-obj.p0))+(obj.gradx*(x3-obj.p0)).^2+(obj.grady*(x3-obj.p0)).^2)*WeirdDiag( diag((2*x1+6*obj.lambda)/(12*obj.mu)),obj.M,obj.N);
%             J31 = J31a + J31b + J31c;
%             J31=sparse(J31);
%             %%%
%             J32 = -diag(x3)*WeirdDiag( diag(1./(3*x1.^3)),obj.M,obj.N)*1/obj.del;
%             J32=sparse(J32);
%             %%%
%             arg331b = (x1+4*obj.lambda)./(4*obj.mu).*(obj.gradx2*(x1-obj.z0));
%             J331 = -DiagHat(((x2/obj.del)./(3*x1.^3)),obj.M,obj.N) + DiagHat(arg331b,obj.M,obj.N)*(diag(obj.gradx*(x3-obj.p0)) + diag(x3)*obj.gradx);
% %             J331 = -DiagHat(((x2/obj.del)./(3*x1.^3)),obj.M,obj.N) + DiagHat(arg331b,obj.M,obj.N)*(diag(obj.gradx*(x3-obj.p0))) + diag(x3)*obj.gradx;
%             %%%
%             J332 = DiagHat(((x1.^2+6*obj.lambda*x1)/(12*obj.mu)),obj.M,obj.N)*(diag(obj.LAP*(x3-obj.p0)) + diag(x3)*obj.LAP + 2*diag(obj.gradx*(x3-obj.p0))*obj.gradx+2*diag(obj.grady*(x3-obj.p0))*obj.grady);
%             %%%
%             J33 = J331 + J332;
%             J33=sparse(J33);
            %%%%%%%%%%%%%%%
            J = [J11 J12 J13 ; J21 J22 J23; J31 J32 J33];
        end
        
        function  [R,J] =  ResJacGNAT(obj,U,t)
            %This function calcuates the residual and jacobian of MEMS
            %problem at the specified indices.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - MEMS object
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
            
            R = zeros(obj.lenRes,1);
            J = zeros(obj.lenJac,1);
            
            rho = obj.rho0*obj.height*obj.width;
            In = obj.In0*obj.width*obj.height^3;
            
            if ~isempty(obj.indLTn) %index <=n
                x1 = U(obj.indLTn,1);
                x2 = U(obj.indLTn+1,1);
                
                R(obj.ResInd.indLTn,1) = (1/3)*x2./(x1.^2)/obj.del;
                J(obj.JacInd.indLTn,1) = reshape([-(2/(3*obj.del))*x2./(x1.^3), 1./(3*obj.del*x1.^2)]',length(obj.JacInd.indLTn),1);
            end
            
            if ~isempty(obj.indEQnp1)
                x1 = U(obj.indEQnp1:obj.indEQnp1+2,1);
                x2 = U(obj.indEQnp1+3,1);
                x3 = U(obj.indEQnp1+4:obj.indEQnp1+(obj.M+3),1);
                
                R(obj.ResInd.indEQnp1,1) = (2/(3*obj.del))*x2*x2/(x1(1).^3) + (3/rho)*(x1(1)^2)*(obj.INT*(x3-obj.p0) + ...
                    obj.S0*obj.height*obj.width*(obj.D2(1,1:3)*(x1-obj.z0)) - ...
                    obj.E*In*(obj.D4(1,1:3)*(x1-obj.z0)))*obj.del;
                
                J21 = [-2*(x2.^2./(x1(1,:).^4))/obj.del + 6/rho*x1(1,:).*( obj.INT*(x3-obj.p0) + obj.S0*obj.height*obj.width*(obj.D2(1,1:2)*(x1(1:2,:)-obj.z0)) - ...
                    obj.E*In*(obj.D4(1,1:3)*(x1-obj.z0)) )*obj.del;0;0] + (3/rho*x1(1,:).^2).*(obj.S0*obj.height*obj.width*obj.D2(1,1:3)' - obj.E*In*obj.D4(1,1:3)')*obj.del;
                
                J22 = (4/3*(x2/obj.del)./(x1(1,:).^3));
                J23 = obj.INT'*(3*obj.del/rho*x1(1,:).^2);
                
                temp = [J21;J22;J23];
                J(obj.JacInd.indEQnp1,1) = temp(:);
            end
            
            if ~isempty(obj.indEQnp2)
                x1 = U(obj.indEQnp2:obj.indEQnp2+3,1);
                x2 = U(obj.indEQnp2+4,1);
                x3 = U(obj.indEQnp2+5:obj.indEQnp2+(obj.M+4),1);
                
                R(obj.ResInd.indEQnp2,1) = (2/(3*obj.del))*x2*x2/(x1(2).^3) + (3/rho)*(x1(2)^2)*(obj.INT*(x3-obj.p0) + ...
                    obj.S0*obj.height*obj.width*(obj.D2(2,:)*(x1(1:3,1)-obj.z0)) - ...
                    obj.E*In*(obj.D4(2,1:4)*(x1-obj.z0)))*obj.del;
                
                J21 = [0;-2*(x2.^2./(x1(2,:).^4))/obj.del + 6/rho*x1(2,:).*( obj.INT*(x3-obj.p0) + obj.S0*obj.height*obj.width*(obj.D2(2,:)*(x1(1:3,:)-obj.z0)) - ...
                    obj.E*In*(obj.D4(2,1:4)*(x1-obj.z0)) )*obj.del;0;0] + (3/rho*x1(2,:).^2).*(obj.S0*obj.height*obj.width*[obj.D2(2,:),0]' - obj.E*In*obj.D4(2,1:4)')*obj.del;
                
                J22 = (4/3*(x2/obj.del)./(x1(2,:).^3));
                J23 = obj.INT'*(3*obj.del/rho*x1(2,:).^2);
                
                temp = [J21;J22;J23];
                J(obj.JacInd.indEQnp2,1) = temp(:);
            end
            
            if ~isempty(obj.indLT2nElse)
                k = length(obj.indLT2nElse);
                x1 = [U(obj.indLT2nElse,1)';U(obj.indLT2nElse+1,1)';U(obj.indLT2nElse+2,1)';U(obj.indLT2nElse+3,1)';U(obj.indLT2nElse+4,1)'];
                x2 = U(obj.indLT2nElse+5,1)';
                ind = repmat(obj.indLT2nElse',obj.M,1) + repmat([6:(obj.M+5)]',1,k);
                x3 = U(ind);
                %                 x3 = [U(obj.indLT2nElse+6,1)';U(obj.indLT2nElse+7,1)';U(obj.indLT2nElse+8,1)';U(obj.indLT2nElse+9,1)';U(obj.indLT2nElse+10,1)'];
                
                R(obj.ResInd.indLT2nElse,1) = ((2/(3*obj.del))*(x2.^2)./(x1(3,:).^3) + ...
                    (3/rho)*(x1(3,:).^2).*((obj.INT*(x3-obj.p0)) + ...
                    obj.S0*obj.height*obj.width*(obj.D2(2,1:3)*(x1(2:4,:)-obj.z0)) - ...
                    obj.E*In*(obj.D4(3,:)*(x1-obj.z0)))*obj.del)';
                
                J21 = [zeros(2,k);-2*(x2.^2./(x1(3,:).^4))/obj.del + 6/rho*x1(3,:).*( obj.INT*(x3-obj.p0) + obj.S0*obj.height*obj.width*(obj.D2(2,:)*(x1(2:4,:)-obj.z0)) - ...
                    obj.E*In*(obj.D4(3,:)*(x1-obj.z0)) )*obj.del;zeros(2,k)] + (obj.S0*obj.height*obj.width*[0,obj.D2(2,:),0]' - obj.E*In*obj.D4(3,:)')*(3/rho*(x1(3,:)).^2)*obj.del;
                
                J22 = (4/3*(x2/obj.del)./(x1(3,:).^3));
                J23 = obj.INT'*(3*obj.del/rho*x1(3,:).^2);
                
                temp = [J21;J22;J23];
                J(obj.JacInd.indLT2nElse,1) = temp(:);
            end
            
            if ~isempty(obj.indEQ2ns1)
                x1 = U(obj.indEQ2ns1:obj.indEQ2ns1+3,1);
                x2 = U(obj.indEQ2ns1+4,1);
                x3 = U(obj.indEQ2ns1+5:obj.indEQ2ns1+(obj.M+4),1);
                
                R(obj.ResInd.indEQ2ns1,1) = (2/(3*obj.del))*x2*x2/(x1(3).^3) + (3/rho)*(x1(3)^2)*(obj.INT*(x3-obj.p0) + ...
                    obj.S0*obj.height*obj.width*(obj.D2(2,1:3)*(x1(2:4,1)-obj.z0)) - ...
                    obj.E*In*(obj.D4(4,2:5)*(x1-obj.z0)))*obj.del;
                
                J21 = [0;0;-2*(x2.^2./(x1(3,:).^4))/obj.del + 6/rho*x1(3,:).*( obj.INT*(x3-obj.p0) + obj.S0*obj.height*obj.width*(obj.D2(2,:)*(x1(2:4,:)-obj.z0)) - ...
                    obj.E*In*(obj.D4(4,2:5)*(x1-obj.z0)) )*obj.del;0] + (3/rho*x1(3,:).^2).*(obj.S0*obj.height*obj.width*[0,obj.D2(2,:)]' - obj.E*In*obj.D4(4,2:5)')*obj.del;
                
                J22 = (4/3*(x2/obj.del)./(x1(3,:).^3));
                J23 = obj.INT'*(3*obj.del/rho*x1(3,:).^2);
                
                temp = [J21;J22;J23];
                J(obj.JacInd.indEQ2ns1,1) = temp(:);
            end
            
            if ~isempty(obj.indEQ2n)
                x1 = U(obj.indEQ2n:obj.indEQ2n+2,1);
                x2 = U(obj.indEQ2n+3,1);
                x3 = U(obj.indEQ2n+4:obj.indEQ2n+(obj.M+3),1);
                
                R(obj.ResInd.indEQ2n,1) = (2/(3*obj.del))*x2*x2/(x1(3).^3) + (3/rho)*(x1(3)^2)*(obj.INT*(x3-obj.p0) + ...
                    obj.S0*obj.height*obj.width*(obj.D2(3,1:3)*(x1-obj.z0)) - ...
                    obj.E*In*(obj.D4(5,3:5)*(x1-obj.z0)))*obj.del;
                
                J21 = [0;0;-2*(x2.^2./(x1(3,:).^4))/obj.del + 6/rho*x1(3,:).*( obj.INT*(x3-obj.p0) + obj.S0*obj.height*obj.width*(obj.D2(3,2:3)*(x1(2:3,:)-obj.z0)) - ...
                    obj.E*In*(obj.D4(5,3:5)*(x1-obj.z0)) )*obj.del] + (3/rho*x1(3,:).^2).*(obj.S0*obj.height*obj.width*[0,obj.D2(3,2:3)]' - obj.E*In*obj.D4(5,3:5)')*obj.del;
                
                J22 = (4/3*(x2/obj.del)./(x1(3,:).^3));
                J23 = obj.INT'*(3*obj.del/rho*x1(3,:).^2);
                
                temp = [J21;J22;J23];
                J(obj.JacInd.indEQ2n,1) = temp(:);
            end
            
            
            if ~isempty(obj.indLT2npm.top)
                k = length(obj.indLT2npm.top);
                x1 = U(obj.indLT2npm.top:obj.indLT2npm.top+1,1);
                x2 = U(obj.indLT2npm.top+2,1);
                x3 = U(obj.indLT2npm.top+3:obj.indLT2npm.top+5,1);
                
                R(obj.ResInd.indLT2npmTop,1) = -x2/(3*x1(1).^3)*x3(1)/obj.del + ...
                    ((x1(1) + 4*obj.lambda)/(4*obj.mu))*(obj.gradx2(1,1:2)*(x1-obj.z0))*(x3(1)*(obj.gradx*(x3([1;3],1)-obj.p0))) + ...
                    ((x1(1)^2 + 6*obj.lambda*x1(1))/(12*obj.mu))*(x3(1)*(obj.LAP(1,3:end)*(x3-obj.p0)) + (obj.gradx*(x3([1;3],1)-obj.p0))^2 + (obj.grady(1,2)*(x3(2)-obj.p0))^2);
                
%                 R(obj.ResInd.indLT2npmTop,1) = -x2/(3*x1(2).^3)*x3(2)/obj.del + ...
%                     ((x1(2) + 4*obj.lambda)/(4*obj.mu))*(obj.gradx2(3,2:3)*(x1-obj.z0))*(x3(2)*obj.gradx*(x3(1:2,1)-obj.p0)) + ...
%                     ((x1(2)^2 + 6*obj.lambda*x1(2))/(12*obj.mu))*(x3(2)*obj.LAP(7,[1,3,4])*(x3-obj.p0) + (obj.gradx*(x3(1:2,1)-obj.p0))^2 + (obj.grady(1,2)*(x3(3)-obj.p0))^2);
                
                J1 = [(x2./x1(1,:).^4).*x3(1,:)/obj.del; zeros(1,k); -x3(1,:)./(3*x1(1,:).^3)/obj.del; -x2./(3*x1(1,:).^3)/obj.del; zeros(2,k)];
                
                J2 = [ obj.gradx2(1,1:2)'*(((x1(1,:)+4*obj.lambda)/(4*obj.mu)).*x3(1,:).*(obj.gradx*(x3([1,3],:) - obj.p0)));...
                    zeros(1,k);...
                    [obj.gradx(1),0,obj.gradx(2)]'*(((x1(1,:)+4*obj.lambda)/(4*obj.mu)).*(obj.gradx2(1,1:2)*(x1-obj.z0)).*x3(1,:))];
                J2(1,:) = J2(1,:) + (0.25/obj.mu)*(obj.gradx2(1,1:2)*(x1-obj.z0)).*x3(1,:).*(obj.gradx*(x3([1,3],:)-obj.p0));
                J2(4,:) = J2(4,:) + ((x1(1,:)+4*obj.lambda)/(4*obj.mu)).*(obj.gradx2(1,1:2)*(x1-obj.z0)).*(obj.gradx*(x3([1,3],:)-obj.p0));
                
                x1Adj = ((x1(1,:).^2 + 6*obj.lambda*x1(1,:))/(12*obj.mu));
                op = x3(1,:).*(obj.LAP(1,3:5)*(x3 - obj.p0)) + (obj.gradx*(x3([1,3],:) - obj.p0)).^2 + (obj.grady(1,2)*(x3(2,:) - obj.p0)).^2;
                op1 = obj.LAP(1,3:5)'*(x1Adj.*x3(1,:)) + 2*[obj.gradx(1),0,obj.gradx(2)]'*((obj.gradx*(x3([1,3],:) - obj.p0)).*x1Adj) + 2*[0,obj.grady(1,2),0]'*((obj.grady(1,2)*(x3(2,:) - obj.p0)).*x1Adj);
                J3 = [((x1(1,:) + 3*obj.lambda)/(6*obj.mu)).*op; zeros(2,k); op1];
                J3(4,:) = J3(4,:) + x1Adj.*(obj.LAP(1,3:5)*(x3 - obj.p0));
                
                temp = J1 + J2 + J3;
                J(obj.JacInd.indLT2npmTop,1) = temp(:);
            end
            
            if ~isempty(obj.indLT2npm.mid)
                k = length(obj.indLT2npm.mid);
                x1 = [U(obj.indLT2npm.mid,1)';U(obj.indLT2npm.mid+1,1)'];
                x2 = U(obj.indLT2npm.mid+2,1)';
                x3 = [U(obj.indLT2npm.mid+3,1)';U(obj.indLT2npm.mid+4,1)';U(obj.indLT2npm.mid+5,1)';U(obj.indLT2npm.mid+6,1)'];
                
                f = -x2./(3*x1(1,:).^3).*x3(2,:)/obj.del + ...
                    ((x1(1,:) + 4*obj.lambda)/(4*obj.mu)).*(obj.gradx2(1,1:2)*(x1-obj.z0)).*(x3(2,:).*(obj.gradx*(x3([2;4],:)-obj.p0))) + ...
                    ((x1(1,:).^2 + 6*obj.lambda*x1(1,:))/(12*obj.mu)).*(x3(2,:).*(obj.LAP(2,2:5)*(x3-obj.p0)) + (obj.gradx*(x3([2;4],:)-obj.p0)).^2 + (obj.grady(2,:)*(x3([1;3],:)-obj.p0)).^2);
                R(obj.ResInd.indLT2npmMid,1)=f';
                
                J1 = [(x2./x1(1,:).^4).*x3(2,:)/obj.del; zeros(1,k); -x3(2,:)./(3*x1(1,:).^3)/obj.del; zeros(1,k); -x2./(3*x1(1,:).^3)/obj.del; zeros(2,k)];
                
                J2 = [ obj.gradx2(1,1:2)'*(((x1(1,:)+4*obj.lambda)/(4*obj.mu)).*x3(2,:).*(obj.gradx*(x3([2,4],:) - obj.p0)));...
                    zeros(1,k);...
                    [0,obj.gradx(1),0,obj.gradx(2)]'*(((x1(1,:)+4*obj.lambda)/(4*obj.mu)).*(obj.gradx2(1,1:2)*(x1-obj.z0)).*x3(2,:))];
                J2(1,:) = J2(1,:) + (0.25/obj.mu)*(obj.gradx2(1,1:2)*(x1-obj.z0)).*x3(2,:).*(obj.gradx*(x3([2,4],:)-obj.p0));
                J2(5,:) = J2(5,:) + ((x1(1,:)+4*obj.lambda)/(4*obj.mu)).*(obj.gradx2(1,1:2)*(x1-obj.z0)).*(obj.gradx*(x3([2,4],:)-obj.p0));
                
                x1Adj = ((x1(1,:).^2 + 6*obj.lambda*x1(1,:))/(12*obj.mu));
                op = x3(2,:).*(obj.LAP(2,2:5)*(x3 - obj.p0)) + (obj.gradx*(x3([2,4],:) - obj.p0)).^2 + (obj.grady(2,:)*(x3([1,3],:) - obj.p0)).^2;
                op1 = obj.LAP(2,2:5)'*(x1Adj.*x3(2,:)) + 2*[0,obj.gradx(1),0,obj.gradx(2)]'*((obj.gradx*(x3([2,4],:) - obj.p0)).*x1Adj) + 2*[obj.grady(2,1),0,obj.grady(2,2),0]'*((obj.grady(2,:)*(x3([1,3],:) - obj.p0)).*x1Adj);
                J3 = [((x1(1,:) + 3*obj.lambda)/(6*obj.mu)).*op; zeros(2,k); op1];
                J3(5,:) = J3(5,:) + x1Adj.*(obj.LAP(2,2:5)*(x3 - obj.p0));
                
                temp = J1 + J2 + J3;
                J(obj.JacInd.indLT2npmMid,1) = temp(:);
            end
            
            if ~isempty(obj.indLT2npm.bot)
                k = length(obj.indLT2npm.bot);
                x1 = U(obj.indLT2npm.bot:obj.indLT2npm.bot+1,1);
                x2 = U(obj.indLT2npm.bot+2,1);
                x3 = U(obj.indLT2npm.bot+3:obj.indLT2npm.bot+5,1);
                
                R(obj.ResInd.indLT2npmBot,1) = -x2/(3*x1(1).^3)*x3(2)/obj.del + ...
                    ((x1(1) + 4*obj.lambda)/(4*obj.mu))*(obj.gradx2(1,1:2)*(x1-obj.z0))*(x3(2)*(obj.gradx*(x3([2;3],1)-obj.p0))) + ...
                    ((x1(1)^2 + 6*obj.lambda*x1(1))/(12*obj.mu))*(x3(2)*obj.LAP(3,[2;3;5])*(x3-obj.p0) + (obj.gradx*(x3([2;3],1)-obj.p0))^2 + (obj.grady(3,1)*(x3(1)-obj.p0))^2);
                
                
                J1 = [(x2./x1(1,:).^4).*x3(2,:)/obj.del; zeros(1,k); -x3(2,:)./(3*x1(1,:).^3)/obj.del; zeros(1,k); -x2./(3*x1(1,:).^3)/obj.del;zeros(1,k)];
                
                J2 = [ obj.gradx2(1,1:2)'*(((x1(1,:)+4*obj.lambda)/(4*obj.mu)).*x3(2,:).*(obj.gradx*(x3([2,3],:) - obj.p0)));...
                    zeros(1,k);...
                    [0,obj.gradx(1),obj.gradx(2)]'*(((x1(1,:)+4*obj.lambda)/(4*obj.mu)).*(obj.gradx2(1,1:2)*(x1-obj.z0)).*x3(2,:))];
                J2(1,:) = J2(1,:) + (0.25/obj.mu)*(obj.gradx2(1,1:2)*(x1-obj.z0)).*x3(2,:).*(obj.gradx*(x3([2,3],:)-obj.p0));
                J2(5,:) = J2(5,:) + ((x1(1,:)+4*obj.lambda)/(4*obj.mu)).*(obj.gradx2(1,1:2)*(x1-obj.z0)).*(obj.gradx*(x3([2,3],:)-obj.p0));
                
                x1Adj = ((x1(1,:).^2 + 6*obj.lambda*x1(1,:))/(12*obj.mu));
                op = x3(2,:).*(obj.LAP(3,[2,3,5])*(x3 - obj.p0)) + (obj.gradx*(x3([2,3],:) - obj.p0)).^2 + (obj.grady(3,1)*(x3(1,:) - obj.p0)).^2;
                op1 = obj.LAP(3,[2,3,5])'*(x1Adj.*x3(2,:)) + 2*[0,obj.gradx(1),obj.gradx(2)]'*((obj.gradx*(x3([2,3],:) - obj.p0)).*x1Adj) + 2*[obj.grady(3,1),0,0]'*((obj.grady(3,1)*(x3(1,:) - obj.p0)).*x1Adj);
                J3 = [((x1(1,:) + 3*obj.lambda)/(6*obj.mu)).*op; zeros(2,k); op1];
                J3(5,:) = J3(5,:) + x1Adj.*(obj.LAP(3,[2,3,5])*(x3 - obj.p0));
                
                temp = J1 + J2 + J3;
                J(obj.JacInd.indLT2npmBot,1) = temp(:);
                
            end
            
            
            if ~isempty(obj.indGTmnp2nsmp1.top)
                k = length(obj.indGTmnp2nsmp1.top);
                x1 = U(obj.indGTmnp2nsmp1.top:obj.indGTmnp2nsmp1.top+1,1);
                x2 = U(obj.indGTmnp2nsmp1.top+2,1);
                x3 = U(obj.indGTmnp2nsmp1.top+3:obj.indGTmnp2nsmp1.top+5,1);
                
                R(obj.ResInd.indGTmnp2nsmp1Top,1) = -x2/(3*x1(2).^3)*x3(2)/obj.del + ...
                    ((x1(2) + 4*obj.lambda)/(4*obj.mu))*(obj.gradx2(3,2:3)*(x1-obj.z0))*(x3(2)*obj.gradx*(x3(1:2,1)-obj.p0)) + ...
                    ((x1(2)^2 + 6*obj.lambda*x1(2))/(12*obj.mu))*(x3(2)*obj.LAP(7,[1,3,4])*(x3-obj.p0) + (obj.gradx*(x3(1:2,1)-obj.p0))^2 + (obj.grady(1,2)*(x3(3)-obj.p0))^2);
                
                J1 = [zeros(1,k);(x2./x1(2,:).^4).*x3(2,:)/obj.del; -x3(2,:)./(3*x1(2,:).^3)/obj.del; zeros(1,k); -x2./(3*x1(2,:).^3)/obj.del; zeros(1,k)];
                
                J2 = [ obj.gradx2(3,2:3)'*(((x1(2,:)+4*obj.lambda)/(4*obj.mu)).*x3(2,:).*(obj.gradx*(x3(1:2,:) - obj.p0)));...
                    zeros(1,k);...
                    [obj.gradx(1),obj.gradx(2),0]'*(((x1(2,:)+4*obj.lambda)/(4*obj.mu)).*(obj.gradx2(3,2:3)*(x1-obj.z0)).*x3(2,:))];
                J2(2,:) = J2(2,:) + (0.25/obj.mu)*(obj.gradx2(3,2:3)*(x1-obj.z0)).*x3(2,:).*(obj.gradx*(x3(1:2,:)-obj.p0));
                J2(5,:) = J2(5,:) + ((x1(2,:)+4*obj.lambda)/(4*obj.mu)).*(obj.gradx2(3,2:3)*(x1-obj.z0)).*(obj.gradx*(x3(1:2,:)-obj.p0));
                
                x1Adj = ((x1(2,:).^2 + 6*obj.lambda*x1(2,:))/(12*obj.mu));
                op = x3(2,:).*(obj.LAP(7,[1,3,4])*(x3 - obj.p0)) + (obj.gradx*(x3(1:2,:) - obj.p0)).^2 + (obj.grady(1,2)*(x3(3,:) - obj.p0)).^2;
                op1 = obj.LAP(7,[1,3,4])'*(x1Adj.*x3(2,:)) + 2*[obj.gradx(1),obj.gradx(2),0]'*((obj.gradx*(x3(1:2,:) - obj.p0)).*x1Adj) + 2*[0,0,obj.grady(2,2)]'*((obj.grady(2,2)*(x3(3,:) - obj.p0)).*x1Adj);
                J3 = [zeros(1,k); ((x1(2,:) + 3*obj.lambda)/(6*obj.mu)).*op; zeros(1,k); op1];
                J3(5,:) = J3(5,:) + x1Adj.*(obj.LAP(7,[1,3,4])*(x3 - obj.p0));
                
                temp = J1 + J2 + J3;
                J(obj.JacInd.indGTmnp2nsmp1Top,1) = temp(:);
                
            end
            
            if ~isempty(obj.indGTmnp2nsmp1.mid)
                k = length(obj.indGTmnp2nsmp1.mid);
                x1 = [U(obj.indGTmnp2nsmp1.mid,1)'; U(obj.indGTmnp2nsmp1.mid+1,1)'];
                x2 = U(obj.indGTmnp2nsmp1.mid+2,1)';
                x3 = [U(obj.indGTmnp2nsmp1.mid+3,1)';U(obj.indGTmnp2nsmp1.mid+4,1)';U(obj.indGTmnp2nsmp1.mid+5,1)';U(obj.indGTmnp2nsmp1.mid+6,1)'];
                
                f = -x2./(3*x1(2,:).^3).*x3(3,:)/obj.del + ...
                    ((x1(2,:) + 4*obj.lambda)/(4*obj.mu)).*(obj.gradx2(3,2:3)*(x1-obj.z0)).*((x3(3,:).*(obj.gradx*(x3([1,3],:)-obj.p0)))) + ...
                    ((x1(2,:).^2 + 6*obj.lambda*x1(2,:))/(12*obj.mu)).*(x3(3,:).*(obj.LAP(8,1:4)*(x3-obj.p0)) + (obj.gradx*(x3([1,3],:)-obj.p0)).^2 + (obj.grady(2,:)*(x3([2,4],:)-obj.p0)).^2);
                R(obj.ResInd.indGTmnp2nsmp1Mid,1)=f';
                
                J1 = [zeros(1,k);(x2./x1(2,:).^4).*x3(3,:)/obj.del; -x3(3,:)./(3*x1(2,:).^3)/obj.del; zeros(2,k); -x2./(3*x1(2,:).^3)/obj.del; zeros(1,k)];
                
                J2 = [ obj.gradx2(3,2:3)'*(((x1(2,:)+4*obj.lambda)/(4*obj.mu)).*x3(3,:).*(obj.gradx*(x3([1,3],:) - obj.p0)));...
                    zeros(1,k);...
                    [obj.gradx(1),0,obj.gradx(2),0]'*(((x1(2,:)+4*obj.lambda)/(4*obj.mu)).*(obj.gradx2(3,2:3)*(x1-obj.z0)).*x3(3,:))];
                J2(2,:) = J2(2,:) + (0.25/obj.mu)*(obj.gradx2(3,2:3)*(x1-obj.z0)).*x3(3,:).*(obj.gradx*(x3([1,3],:)-obj.p0));
                J2(6,:) = J2(6,:) + ((x1(2,:)+4*obj.lambda)/(4*obj.mu)).*(obj.gradx2(3,2:3)*(x1-obj.z0)).*(obj.gradx*(x3([1,3],:)-obj.p0));
                
                x1Adj = ((x1(2,:).^2 + 6*obj.lambda*x1(2,:))/(12*obj.mu));
                op = x3(3,:).*(obj.LAP(8,1:4)*(x3 - obj.p0)) + (obj.gradx*(x3([1,3],:) - obj.p0)).^2 + (obj.grady(2,:)*(x3([2,4],:) - obj.p0)).^2;
                op1 = obj.LAP(8,1:4)'*(x1Adj.*x3(3,:)) + 2*[obj.gradx(1),0,obj.gradx(2),0]'*((obj.gradx*(x3([1,3],:) - obj.p0)).*x1Adj) + 2*[0,obj.grady(2,1),0,obj.grady(2,2)]'*((obj.grady(2,:)*(x3([2,4],:) - obj.p0)).*x1Adj);
                J3 = [zeros(1,k); ((x1(2,:) + 3*obj.lambda)/(6*obj.mu)).*op; zeros(1,k); op1];
                J3(6,:) = J3(6,:) + x1Adj.*(obj.LAP(8,1:4)*(x3 - obj.p0));
                
                temp = J1 + J2 + J3;
                J(obj.JacInd.indGTmnp2nsmp1Mid,1) = temp(:);
            end
            
            if ~isempty(obj.indGTmnp2nsmp1.bot)
                k = length(obj.indGTmnp2nsmp1.bot);
                x1 = U(obj.indGTmnp2nsmp1.bot:obj.indGTmnp2nsmp1.bot+1,1);
                x2 = U(obj.indGTmnp2nsmp1.bot+2,1);
                x3 = U(obj.indGTmnp2nsmp1.bot+3:obj.indGTmnp2nsmp1.bot+5,1);
                
                R(obj.ResInd.indGTmnp2nsmp1Bot,1) = -x2/(3*x1(2).^3)*x3(3)/obj.del + ...
                    ((x1(2) + 4*obj.lambda)/(4*obj.mu))*obj.gradx2(3,2:3)*(x1-obj.z0)*(x3(3)*obj.gradx*(x3([1,3],:)-obj.p0)) + ...
                    ((x1(2)^2 + 6*obj.lambda*x1(2))/(12*obj.mu))*(x3(3)*obj.LAP(9,1:3)*(x3-obj.p0) + (obj.gradx*(x3([1,3],:)-obj.p0))^2 + (obj.grady(3,1)*(x3(2)-obj.p0))^2);
                
                J1 = [zeros(1,k);(x2./x1(2,:).^4).*x3(3,:)/obj.del; -x3(3,:)./(3*x1(2,:).^3)/obj.del; zeros(2,k); -x2./(3*x1(2,:).^3)/obj.del];
                
                J2 = [ obj.gradx2(3,2:3)'*(((x1(2,:)+4*obj.lambda)/(4*obj.mu)).*x3(3,:).*(obj.gradx*(x3([1,3],:) - obj.p0)));...
                    zeros(1,k);...
                    [obj.gradx(1),0,obj.gradx(2)]'*(((x1(2,:)+4*obj.lambda)/(4*obj.mu)).*(obj.gradx2(3,2:3)*(x1-obj.z0)).*x3(3,:))];
                J2(2,:) = J2(2,:) + (0.25/obj.mu)*(obj.gradx2(3,2:3)*(x1-obj.z0)).*x3(3,:).*(obj.gradx*(x3([1,3],:)-obj.p0));
                J2(6,:) = J2(6,:) + ((x1(2,:)+4*obj.lambda)/(4*obj.mu)).*(obj.gradx2(3,2:3)*(x1-obj.z0)).*(obj.gradx*(x3([1,3],:)-obj.p0));
                            
                x1Adj = ((x1(2,:).^2 + 6*obj.lambda*x1(2,:))/(12*obj.mu));
                op = x3(3,:).*(obj.LAP(9,1:3)*(x3 - obj.p0)) + (obj.gradx*(x3([1,3],:) - obj.p0)).^2 + (obj.grady(3,1)*(x3(2,:) - obj.p0)).^2;
                op1 = obj.LAP(9,1:3)'*(x1Adj.*x3(3,:)) + 2*[obj.gradx(1),0,obj.gradx(2)]'*((obj.gradx*(x3([1,3],:) - obj.p0)).*x1Adj) + 2*[0,obj.grady(3,1),0]'*((obj.grady(3,1)*(x3(2,:) - obj.p0)).*x1Adj);
                J3 = [zeros(1,k); ((x1(2,:) + 3*obj.lambda)/(6*obj.mu)).*op; zeros(1,k); op1];
                J3(6,:) = J3(6,:) + x1Adj.*(obj.LAP(9,1:3)*(x3 - obj.p0));
                
                temp = J1 + J2 + J3;
                J(obj.JacInd.indGTmnp2nsmp1Bot,1) = temp(:);
            end
            
            if ~isempty(obj.indGTmnElse.top)
                k = length(obj.indGTmnElse.top);
                x1 = [U(obj.indGTmnElse.top,1)';U(obj.indGTmnElse.top+1,1)';U(obj.indGTmnElse.top+2,1)'];
                x2 = U(obj.indGTmnElse.top+3,1)';
                x3 = [U(obj.indGTmnElse.top+4,1)';U(obj.indGTmnElse.top+5,1)';U(obj.indGTmnElse.top+6,1)';U(obj.indGTmnElse.top+7,1)'];
                
                f = -x2./(3*x1(2,:).^3).*x3(2,:)/obj.del + ...
                    ((x1(2,:) + 4*obj.lambda)/(4*obj.mu)).*(obj.gradx2(2,:)*(x1-obj.z0)).*(x3(2,:).*(obj.gradx*(x3([1,4],:)-obj.p0))) + ...
                    ((x1(2,:).^2 + 6*obj.lambda*x1(2,:))/(12*obj.mu)).*(x3(2,:).*(obj.LAP(4,[1,3:5])*(x3-obj.p0)) + (obj.gradx*(x3([1,4],:)-obj.p0)).^2 + (obj.grady(1,2)*(x3(3,:)-obj.p0)).^2);
                R(obj.ResInd.indGTmnElseTop,1)=f';
                
                J1 = [zeros(1,k);(x2./x1(2,:).^4).*x3(2,:)/obj.del;zeros(1,k); -x3(2,:)./(3*x1(2,:).^3)/obj.del; zeros(1,k); -x2./(3*x1(2,:).^3)/obj.del; zeros(2,k)];
                
                J2 = [ obj.gradx2(2,:)'*(((x1(2,:)+4*obj.lambda)/(4*obj.mu)).*x3(2,:).*(obj.gradx*(x3([1,4],:) - obj.p0)));...
                    zeros(1,k);...
                    [obj.gradx(1),zeros(1,2),obj.gradx(2)]'*(((x1(2,:)+4*obj.lambda)/(4*obj.mu)).*(obj.gradx2(2,:)*(x1-obj.z0)).*x3(2,:))];
                J2(2,:) = J2(2,:) + (0.25/obj.mu)*(obj.gradx2(2,:)*(x1-obj.z0)).*x3(2,:).*(obj.gradx*(x3([1,4],:)-obj.p0));
                J2(6,:) = J2(6,:) + ((x1(2,:)+4*obj.lambda)/(4*obj.mu)).*(obj.gradx2(2,:)*(x1-obj.z0)).*(obj.gradx*(x3([1,4],:)-obj.p0));
                
                x1Adj = ((x1(2,:).^2 + 6*obj.lambda*x1(2,:))/(12*obj.mu));
                op = x3(2,:).*(obj.LAP(4,[1,3:5])*(x3 - obj.p0)) + (obj.gradx*(x3([1,4],:) - obj.p0)).^2 + (obj.grady(2,2)*(x3(3,:) - obj.p0)).^2;
                op1 = obj.LAP(4,[1,3:5])'*(x1Adj.*x3(2,:)) + 2*[obj.gradx(1),zeros(1,2),obj.gradx(2)]'*((obj.gradx*(x3([1,4],:) - obj.p0)).*x1Adj) + 2*[0,0,obj.grady(2,2),0]'*((obj.grady(2,2)*(x3(3,:) - obj.p0)).*x1Adj);
                J3 = [zeros(1,k); ((x1(2,:) + 3*obj.lambda)/(6*obj.mu)).*op; zeros(2,k); op1];
                J3(6,:) = J3(6,:) + x1Adj.*(obj.LAP(4,[1,3:5])*(x3 - obj.p0));
                
                temp = J1 + J2 + J3;
                J(obj.JacInd.indGTmnElseTop,1) = temp(:);
            end
            
            if ~isempty(obj.indGTmnElse.mid)
                k = length(obj.indGTmnElse.mid);
                x1 = [U(obj.indGTmnElse.mid,1)';U(obj.indGTmnElse.mid+1,1)';U(obj.indGTmnElse.mid+2,1)'];
                x2 = U(obj.indGTmnElse.mid+3,1)';
                x3 = [U(obj.indGTmnElse.mid+4,1)';U(obj.indGTmnElse.mid+5,1)';U(obj.indGTmnElse.mid+6,1)';U(obj.indGTmnElse.mid+7,1)';U(obj.indGTmnElse.mid+8,1)'];
                
                f = -x2./(3*x1(2,:).^3).*x3(3,:)/obj.del + ...
                    ((x1(2,:) + 4*obj.lambda)/(4*obj.mu)).*(obj.gradx2(2,:)*(x1-obj.z0)).*(x3(3,:).*(obj.gradx*(x3([1,5],:)-obj.p0))) + ...
                    ((x1(2,:).^2 + 6*obj.lambda*x1(2,:))/(12*obj.mu)).*(x3(3,:).*(obj.LAP(5,:)*(x3-obj.p0)) + (obj.gradx*(x3([1,5],:)-obj.p0)).^2 + (obj.grady(2,:)*(x3([2,4],:)-obj.p0)).^2);
                R(obj.ResInd.indGTmnElseMid,1)=f';
                
                J1 = [zeros(1,k);(x2./x1(2,:).^4).*x3(3,:)/obj.del;zeros(1,k); -x3(3,:)./(3*x1(2,:).^3)/obj.del; zeros(2,k); -x2./(3*x1(2,:).^3)/obj.del; zeros(2,k)];
                
                J2 = [ obj.gradx2(2,:)'*(((x1(2,:)+4*obj.lambda)/(4*obj.mu)).*x3(3,:).*(obj.gradx*(x3([1,5],:) - obj.p0)));...
                    zeros(1,k);...
                    [obj.gradx(1),zeros(1,3),obj.gradx(2)]'*(((x1(2,:)+4*obj.lambda)/(4*obj.mu)).*(obj.gradx2(2,:)*(x1-obj.z0)).*x3(3,:))];
                J2(2,:) = J2(2,:) + (0.25/obj.mu)*(obj.gradx2(2,:)*(x1-obj.z0)).*x3(3,:).*(obj.gradx*(x3([1,5],:)-obj.p0));
                J2(7,:) = J2(7,:) + ((x1(2,:)+4*obj.lambda)/(4*obj.mu)).*(obj.gradx2(2,:)*(x1-obj.z0)).*(obj.gradx*(x3([1,5],:)-obj.p0));
                
                x1Adj = ((x1(2,:).^2 + 6*obj.lambda*x1(2,:))/(12*obj.mu));
                op = x3(3,:).*(obj.LAP(5,:)*(x3 - obj.p0)) + (obj.gradx*(x3([1,5],:) - obj.p0)).^2 + (obj.grady(2,:)*(x3([2,4],:) - obj.p0)).^2;
                op1 = obj.LAP(5,:)'*(x1Adj.*x3(3,:)) + 2*[obj.gradx(1),zeros(1,3),obj.gradx(2)]'*((obj.gradx*(x3([1,5],:) - obj.p0)).*x1Adj) + 2*[0,obj.grady(2,1),0,obj.grady(2,2),0]'*((obj.grady(2,:)*(x3([2,4],:) - obj.p0)).*x1Adj);
                J3 = [zeros(1,k); ((x1(2,:) + 3*obj.lambda)/(6*obj.mu)).*op; zeros(2,k); op1];
                J3(7,:) = J3(7,:) + x1Adj.*(obj.LAP(5,:)*(x3 - obj.p0));
                
                temp = J1 + J2 + J3;
                J(obj.JacInd.indGTmnElseMid,1) = temp(:);
            end
            
            if ~isempty(obj.indGTmnElse.bot)
                k = length(obj.indGTmnElse.bot);
                x1 = [U(obj.indGTmnElse.bot,1)';U(obj.indGTmnElse.bot+1,1)';U(obj.indGTmnElse.bot+2,1)'];
                x2 = U(obj.indGTmnElse.bot+3,1)';
                x3 = [U(obj.indGTmnElse.bot+4,1)';U(obj.indGTmnElse.bot+5,1)';U(obj.indGTmnElse.bot+6,1)';U(obj.indGTmnElse.bot+7,1)'];
                
                f = -x2./(3*x1(2,:).^3).*x3(3,:)/obj.del + ...
                    ((x1(2,:) + 4*obj.lambda)/(4*obj.mu)).*(obj.gradx2(2,:)*(x1-obj.z0)).*(x3(3,:).*(obj.gradx*(x3([1,4],:)-obj.p0))) + ...
                    ((x1(2,:).^2 + 6*obj.lambda*x1(2,:))/(12*obj.mu)).*(x3(3,:).*(obj.LAP(6,[1:3,5])*(x3-obj.p0)) + (obj.gradx*(x3([1,4],:)-obj.p0)).^2 + (obj.grady(3,1)*(x3(2,:)-obj.p0)).^2);
                R(obj.ResInd.indGTmnElseBot,1)=f';
                
                J1 = [zeros(1,k);(x2./x1(2,:).^4).*x3(3,:)/obj.del;zeros(1,k); -x3(3,:)./(3*x1(2,:).^3)/obj.del; zeros(2,k); -x2./(3*x1(2,:).^3)/obj.del; zeros(1,k)];
                
                J2 = [ obj.gradx2(2,:)'*(((x1(2,:)+4*obj.lambda)/(4*obj.mu)).*x3(3,:).*(obj.gradx*(x3([1,4],:) - obj.p0)));...
                    zeros(1,k);...
                    [obj.gradx(1),zeros(1,2),obj.gradx(2)]'*(((x1(2,:)+4*obj.lambda)/(4*obj.mu)).*(obj.gradx2(2,:)*(x1-obj.z0)).*x3(3,:))];
                J2(2,:) = J2(2,:) + (0.25/obj.mu)*(obj.gradx2(2,:)*(x1-obj.z0)).*x3(3,:).*(obj.gradx*(x3([1,4],:)-obj.p0));
                J2(7,:) = J2(7,:) + ((x1(2,:)+4*obj.lambda)/(4*obj.mu)).*(obj.gradx2(2,:)*(x1-obj.z0)).*(obj.gradx*(x3([1,4],:)-obj.p0));
                
                x1Adj = ((x1(2,:).^2 + 6*obj.lambda*x1(2,:))/(12*obj.mu));
                op = x3(3,:).*(obj.LAP(6,[1:3,5])*(x3 - obj.p0)) + (obj.gradx*(x3([1,4],:) - obj.p0)).^2 + (obj.grady(2,1)*(x3(2,:) - obj.p0)).^2;
                op1 = obj.LAP(6,[1:3,5])'*(x1Adj.*x3(3,:)) + 2*[obj.gradx(1),zeros(1,2),obj.gradx(2)]'*((obj.gradx*(x3([1,4],:) - obj.p0)).*x1Adj) + 2*[0,obj.grady(2,1),0,0]'*((obj.grady(2,1)*(x3(2,:) - obj.p0)).*x1Adj);
                J3 = [zeros(1,k); ((x1(2,:) + 3*obj.lambda)/(6*obj.mu)).*op; zeros(2,k); op1];
                J3(7,:) = J3(7,:) + x1Adj.*(obj.LAP(6,[1:3,5])*(x3 - obj.p0));
                
                temp = J1 + J2 + J3;
                J(obj.JacInd.indGTmnElseBot,1) = temp(:);
            end
            
            R = R  + obj.B*feval(obj.config.inFunc,t);
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
            %obj - MEMS object
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
            %obj     - MEMS object
            %nodenum - scalar specifying node numbers
            %
            %Outputs:
            %--------
            %ind     - vector specifying the corresponding indices in the
            %          full state vector
            %--------------------------------------------------------------
            
            ind = [];
            for i = 1:length(nodenum)
                if ~isempty(find(obj.mesh.dbc(:,1) == nodenum(i),1))
                    ind = [ind; nan];
                    continue;
                end
                
                if nodenum <= obj.N+2
                    ind = [ind; nodenum(i,1)-1; nodenum(i,1) - 1 + obj.N];
                else
                    temp = floor((nodenum(i,1) - obj.config.nNodes(1))/obj.config.nNodes(2));
                    ind = [ind; nodenum(i,1) + 2*obj.N - obj.config.nNodes(1) - obj.config.nNodes(2) - 2*(temp-1) - 1];
                end
            end
            
        end
        
        function  [JacStruct] = JacobianStructure(obj)
            %This function forms a boolean matrix whose true entries
            %correspond to nonzero entries in the jacobian
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - MEMS object
            %
            %Outputs:
            %--------
            %JacStruct - ndof x ndof boolean matrix whose true values
            %            correspond to nonzero entries in the jacobian
            %--------------------------------------------------------------
            
            tempSV1 = ones(obj.config.ndof,1);
            tempSV2 = rand(obj.config.ndof,1);
            tempSV3 = rand(obj.config.ndof,1);
            [~,J1] = ResJac(obj,tempSV1,rand);
            [~,J2] = ResJac(obj,tempSV2,rand);
            [~,J3] = ResJac(obj,tempSV1,rand);
            
            JacStruct = sparse(J1 | J2 | J3);
        end
                
        function  [probGNAT] = createCopy4GNAT(obj,gnat)
            %This function makes a copy of the current object and only
            %copies the source term and input matrix evaluated at indN
            %rows.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - MEMS object
            %indN    - scalar specifying node numbers
            %
            %Outputs:
            %--------
            %probGNAT - new "reduced" MEMS object
            %--------------------------------------------------------------
            
            probGNAT = MEMS([],obj);
            
            temp = probGNAT.config.param;
            temp2 = probGNAT.config.inFunc;
            temp3 = probGNAT.config.form;
            probGNAT.config = [];
            probGNAT.config.param = temp;
            probGNAT.config.inFunc = temp2;
            probGNAT.config.form = temp3;
            probGNAT.mesh   = [];
            probGNAT.ndim = [];
            probGNAT.Jstruct = [];
            
            probGNAT.INT = ones(1,obj.M)*obj.dy;
            probGNAT.D2 = full([obj.D2(1,1:3); obj.D2(2,1:3); obj.D2(end,end-2:end)]);
            probGNAT.D4 = full([obj.D4(1,1:5); obj.D4(2,1:5); obj.D4(3,1:5);...
                obj.D4(end-1,end-4:end);obj.D4(end,end-4:end)]);
            probGNAT.gradx2 = full([obj.gradx2(1,1:3); obj.gradx2(2,1:3); obj.gradx2(end,end-2:end)]);
            probGNAT.gradx  = [-1, 1]/(2*obj.dx);
            probGNAT.grady  = [0 1; -1 1; -1 0]/(2*obj.dy);
            m = obj.M; n = obj.N;
            probGNAT.LAP = full([0,0,obj.LAP(1,1:2),obj.LAP(1,m+1);...
                0,obj.LAP(2,1:3),obj.LAP(2,m+2);...
                0,obj.LAP(m,m-1:m),0,obj.LAP(m,2*m);...
                obj.LAP(m+1,1),0,obj.LAP(m+1,m+1:m+2),obj.LAP(m+1,2*m+1);...
                obj.LAP(m+2,2),obj.LAP(m+2,m+1:m+3),obj.LAP(m+2,2*m+2);...
                obj.LAP(2*m,m),obj.LAP(2*m,2*m-1:2*m),0,obj.LAP(2*m,3*m);...
                obj.LAP(m*(n-1)+1,m*(n-2)+1),0,obj.LAP(m*(n-1)+1,m*(n-1)+1:m*(n-1)+2),0;...
                obj.LAP(m*(n-1)+2,m*(n-2)+2),obj.LAP(m*(n-1)+2,m*(n-1)+1:m*(n-1)+3),0;...
                obj.LAP(m*n,m*(n-1)),obj.LAP(m*n,m*n-1:m*n),0,0]);
            
            
            %             probGNAT.LAP = (1/(obj.dy^2))*repmat([0 0 -2 1 0; 0 1 -2 1 0; 0 1 -2 0 0],3,1) + ...
            %                            (1/(obj.dx^2))*[repmat([0 0 -1 0 1],3,1); repmat([1 0 -2 0 1],3,1); repmat([1 0 1 0 0],3,1)];
            
            if length(obj.B) > 1
                probGNAT.B = obj.B(gnat.sampleInd,1);
            end
            
            if length(obj.G) > 1
                probGNAT.G = obj.G(gnat.sampleInd,1);
            end
            
            if length(obj.ic) > 1
                probGNAT.ic = obj.ic(unique(gnat.jrow),1);
            end
            
            probGNAT.lenRes = length(gnat.sampleInd);
            probGNAT.lenJac = length(gnat.jrow);
            
            probGNAT.indLTn = [];
            
            probGNAT.indEQnp1 = [];
            probGNAT.indEQnp2 = [];
            probGNAT.indEQ2ns1 = [];
            probGNAT.indEQ2n = [];
            probGNAT.indLT2nElse = [];
            
            probGNAT.indLT2npm.top = [];
            probGNAT.indLT2npm.mid = [];
            probGNAT.indLT2npm.bot = [];
            probGNAT.indGTmnp2nsmp1.top = [];
            probGNAT.indGTmnp2nsmp1.mid = [];
            probGNAT.indGTmnp2nsmp1.bot = [];
            probGNAT.indGTmnElse.top = [];
            probGNAT.indGTmnElse.mid = [];
            probGNAT.indGTmnElse.bot = [];
            
            probGNAT.ResInd.indLTn = [];
            probGNAT.ResInd.indEQnp1 = [];
            probGNAT.ResInd.indEQnp2 = [];
            probGNAT.ResInd.indEQ2ns1 = [];
            probGNAT.ResInd.indEQ2n = [];
            probGNAT.ResInd.indLT2nElse = [];
            probGNAT.ResInd.indLT2npmTop = [];
            probGNAT.ResInd.indLT2npmMid = [];
            probGNAT.ResInd.indLT2npmBot = [];
            probGNAT.ResInd.indGTmnp2nsmp1Top = [];
            probGNAT.ResInd.indGTmnp2nsmp1Mid = [];
            probGNAT.ResInd.indGTmnp2nsmp1Bot = [];
            probGNAT.ResInd.indGTmnElseTop = [];
            probGNAT.ResInd.indGTmnElseMid = [];
            probGNAT.ResInd.indGTmnElseBot = [];
            
            probGNAT.JacInd.indLTn = [];
            probGNAT.JacInd.indEQnp1 = [];
            probGNAT.JacInd.indEQnp2 = [];
            probGNAT.JacInd.indEQ2ns1 = [];
            probGNAT.JacInd.indEQ2n = [];
            probGNAT.JacInd.indLT2nElse = [];
            probGNAT.JacInd.indLT2npmTop = [];
            probGNAT.JacInd.indLT2npmMid = [];
            probGNAT.JacInd.indLT2npmBot = [];
            probGNAT.JacInd.indGTmnp2nsmp1Top = [];
            probGNAT.JacInd.indGTmnp2nsmp1Mid = [];
            probGNAT.JacInd.indGTmnp2nsmp1Bot = [];
            probGNAT.JacInd.indGTmnElseTop = [];
            probGNAT.JacInd.indGTmnElseMid = [];
            probGNAT.JacInd.indGTmnElseBot = [];
            
            cnt = 1;
            for i = 1:length(gnat.sampleInd)
                k = length(gnat.irstart(i):gnat.irstart(i+1)-1);
                if gnat.sampleInd(i) <= obj.N
                    probGNAT.indLTn = [probGNAT.indLTn; gnat.irstart(i)];
                    probGNAT.ResInd.indLTn = [probGNAT.ResInd.indLTn; i];
                    probGNAT.JacInd.indLTn = [probGNAT.JacInd.indLTn; (cnt:cnt+k-1)'];
                elseif gnat.sampleInd(i) == obj.N+1
                    probGNAT.indEQnp1 = [probGNAT.indEQnp1; gnat.irstart(i)];
                    probGNAT.ResInd.indEQnp1 = [probGNAT.ResInd.indEQnp1; i];
                    probGNAT.JacInd.indEQnp1 = [probGNAT.JacInd.indEQnp1; (cnt:cnt+k-1)'];
                elseif gnat.sampleInd(i) == obj.N+2
                    probGNAT.indEQnp2 = [probGNAT.indEQnp2; gnat.irstart(i)];
                    probGNAT.ResInd.indEQnp2 = [probGNAT.ResInd.indEQnp2; i];
                    probGNAT.JacInd.indEQnp2 = [probGNAT.JacInd.indEQnp2; (cnt:cnt+k-1)'];
                elseif gnat.sampleInd(i) == 2*obj.N-1
                    probGNAT.indEQ2ns1 = [probGNAT.indEQ2ns1; gnat.irstart(i)];
                    probGNAT.ResInd.indEQ2ns1 = [probGNAT.ResInd.indEQ2ns1; i];
                    probGNAT.JacInd.indEQ2ns1 = [probGNAT.JacInd.indEQ2ns1; (cnt:cnt+k-1)'];
                elseif gnat.sampleInd(i) == 2*obj.N
                    probGNAT.indEQ2n = [probGNAT.indEQ2n; gnat.irstart(i)];
                    probGNAT.ResInd.indEQ2n = [probGNAT.ResInd.indEQ2n; i];
                    probGNAT.JacInd.indEQ2n = [probGNAT.JacInd.indEQ2n; (cnt:cnt+k-1)'];
                elseif gnat.sampleInd(i) <= 2*obj.N
                    probGNAT.indLT2nElse = [probGNAT.indLT2nElse; gnat.irstart(i)];
                    probGNAT.ResInd.indLT2nElse = [probGNAT.ResInd.indLT2nElse; i];
                    probGNAT.JacInd.indLT2nElse = [probGNAT.JacInd.indLT2nElse; (cnt:cnt+k-1)'];
                elseif gnat.sampleInd(i) <= 2*obj.N + obj.M
                    if gnat.sampleInd(i) == 2*obj.N + 1
                        probGNAT.indLT2npm.top = [probGNAT.indLT2npm.top; gnat.irstart(i)];
                        probGNAT.ResInd.indLT2npmTop = [probGNAT.ResInd.indLT2npmTop; i];
                        probGNAT.JacInd.indLT2npmTop = [probGNAT.JacInd.indLT2npmTop; (cnt:cnt+k-1)'];
                    elseif gnat.sampleInd(i) == 2*obj.N + obj.M
                        probGNAT.indLT2npm.bot = [probGNAT.indLT2npm.bot; gnat.irstart(i)];
                        probGNAT.ResInd.indLT2npmBot = [probGNAT.ResInd.indLT2npmBot; i];
                        probGNAT.JacInd.indLT2npmBot = [probGNAT.JacInd.indLT2npmBot; (cnt:cnt+k-1)'];
                    else
                        probGNAT.indLT2npm.mid = [probGNAT.indLT2npm.mid; gnat.irstart(i)];
                        probGNAT.ResInd.indLT2npmMid = [probGNAT.ResInd.indLT2npmMid; i];
                        probGNAT.JacInd.indLT2npmMid = [probGNAT.JacInd.indLT2npmMid; (cnt:cnt+k-1)'];
                    end
                elseif gnat.sampleInd(i) >= obj.M*obj.N +2*obj.N - obj.M + 1
                    if gnat.sampleInd(i) == obj.M*obj.N +2*obj.N - obj.M + 1
                        probGNAT.indGTmnp2nsmp1.top = [probGNAT.indGTmnp2nsmp1.top; gnat.irstart(i)];
                        probGNAT.ResInd.indGTmnp2nsmp1Top = [probGNAT.ResInd.indGTmnp2nsmp1Top; i];
                        probGNAT.JacInd.indGTmnp2nsmp1Top = [probGNAT.JacInd.indGTmnp2nsmp1Top; (cnt:cnt+k-1)'];
                    elseif gnat.sampleInd(i) == obj.M*obj.N +2*obj.N
                        probGNAT.indGTmnp2nsmp1.bot = [probGNAT.indGTmnp2nsmp1.bot; gnat.irstart(i)];
                        probGNAT.ResInd.indGTmnp2nsmp1Bot = [probGNAT.ResInd.indGTmnp2nsmp1Bot; i];
                        probGNAT.JacInd.indGTmnp2nsmp1Bot = [probGNAT.JacInd.indGTmnp2nsmp1Bot; (cnt:cnt+k-1)'];
                    else
                        probGNAT.indGTmnp2nsmp1.mid = [probGNAT.indGTmnp2nsmp1.mid; gnat.irstart(i)];
                        probGNAT.ResInd.indGTmnp2nsmp1Mid = [probGNAT.ResInd.indGTmnp2nsmp1Mid; i];
                        probGNAT.JacInd.indGTmnp2nsmp1Mid = [probGNAT.JacInd.indGTmnp2nsmp1Mid; (cnt:cnt+k-1)'];
                    end
                else
                    if ~isempty(intersect(gnat.sampleInd(i),2*obj.N+1:obj.M:2*obj.N+obj.M*obj.N))
                        probGNAT.indGTmnElse.top = [probGNAT.indGTmnElse.top; gnat.irstart(i)];
                        probGNAT.ResInd.indGTmnElseTop = [probGNAT.ResInd.indGTmnElseTop; i];
                        probGNAT.JacInd.indGTmnElseTop = [probGNAT.JacInd.indGTmnElseTop; (cnt:cnt+k-1)'];
                    elseif ~isempty(intersect(gnat.sampleInd(i),2*obj.N+obj.M:obj.M:2*obj.N+obj.M*obj.N))
                        probGNAT.indGTmnElse.bot = [probGNAT.indGTmnElse.bot; gnat.irstart(i)];
                        probGNAT.ResInd.indGTmnElseBot = [probGNAT.ResInd.indGTmnElseBot; i];
                        probGNAT.JacInd.indGTmnElseBot = [probGNAT.JacInd.indGTmnElseBot; (cnt:cnt+k-1)'];
                    else
                        probGNAT.indGTmnElse.mid = [probGNAT.indGTmnElse.mid; gnat.irstart(i)];
                        probGNAT.ResInd.indGTmnElseMid = [probGNAT.ResInd.indGTmnElseMid; i];
                        probGNAT.JacInd.indGTmnElseMid = [probGNAT.JacInd.indGTmnElseMid; (cnt:cnt+k-1)'];
                    end
                end
                cnt = cnt + k;
            end
        end
        
        function  [probGNAT] = createCopy4GNATloc(obj,gnat)
            %This function makes a copy of the current object and only
            %copies the source term and input matrix evaluated at indN
            %rows.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - MEMS object
            %indN    - scalar specifying node numbers
            %
            %Outputs:
            %--------
            %probGNAT - new "reduced" MEMS object
            %--------------------------------------------------------------
            
            for kk = 1:gnat.nBases
                probGNAT(kk) = MEMS([],obj);
                
                temp = probGNAT(kk).config.param;
                temp2 = probGNAT(kk).config.inFunc;
                temp3 = probGNAT(kk).config.form;
                probGNAT(kk).config = [];
                probGNAT(kk).config.param = temp;
                probGNAT(kk).config.inFunc = temp2;
                probGNAT(kk).config.form = temp3;
                probGNAT(kk).mesh   = [];
                probGNAT(kk).ndim = [];
                probGNAT(kk).Jstruct = [];
                
                probGNAT(kk).INT = ones(1,obj.M)*obj.dy;
                probGNAT(kk).D2 = full([obj.D2(1,1:3); obj.D2(2,1:3); obj.D2(end,end-2:end)]);
                probGNAT(kk).D4 = full([obj.D4(1,1:5); obj.D4(2,1:5); obj.D4(3,1:5);...
                    obj.D4(end-1,end-4:end);obj.D4(end,end-4:end)]);
                probGNAT(kk).gradx2 = full([obj.gradx2(1,1:3); obj.gradx2(2,1:3); obj.gradx2(end,end-2:end)]);
                probGNAT(kk).gradx  = [-1, 1]/(2*obj.dx);
                probGNAT(kk).grady  = [0 1; -1 1; -1 0]/(2*obj.dy);
                m = obj.M; n = obj.N;
                probGNAT(kk).LAP = full([0,0,obj.LAP(1,1:2),obj.LAP(1,m+1);...
                    0,obj.LAP(2,1:3),obj.LAP(2,m+2);...
                    0,obj.LAP(m,m-1:m),0,obj.LAP(m,2*m);...
                    obj.LAP(m+1,1),0,obj.LAP(m+1,m+1:m+2),obj.LAP(m+1,2*m+1);...
                    obj.LAP(m+2,2),obj.LAP(m+2,m+1:m+3),obj.LAP(m+2,2*m+2);...
                    obj.LAP(2*m,m),obj.LAP(2*m,2*m-1:2*m),0,obj.LAP(2*m,3*m);...
                    obj.LAP(m*(n-1)+1,m*(n-2)+1),0,obj.LAP(m*(n-1)+1,m*(n-1)+1:m*(n-1)+2),0;...
                    obj.LAP(m*(n-1)+2,m*(n-2)+2),obj.LAP(m*(n-1)+2,m*(n-1)+1:m*(n-1)+3),0;...
                    obj.LAP(m*n,m*(n-1)),obj.LAP(m*n,m*n-1:m*n),0,0]);
                
                
                %             probGNAT.LAP = (1/(obj.dy^2))*repmat([0 0 -2 1 0; 0 1 -2 1 0; 0 1 -2 0 0],3,1) + ...
                %                            (1/(obj.dx^2))*[repmat([0 0 -1 0 1],3,1); repmat([1 0 -2 0 1],3,1); repmat([1 0 1 0 0],3,1)];
                
                if length(obj.B) > 1
                    probGNAT(kk).B = obj.B(gnat.sampleInd{kk,1},1);
                end
                
                if length(obj.G) > 1
                    probGNAT(kk).G = obj.G(gnat.sampleInd{kk,1},1);
                end
                
                if length(obj.ic) > 1
                    probGNAT(kk).ic = obj.ic(unique(gnat.jrow{kk,1}),1);
                end
                
                probGNAT(kk).lenRes = length(gnat.sampleInd{kk,1});
                probGNAT(kk).lenJac = length(gnat.jrow{kk,1});
                
                probGNAT(kk).indLTn = [];
                
                probGNAT(kk).indEQnp1 = [];
                probGNAT(kk).indEQnp2 = [];
                probGNAT(kk).indEQ2ns1 = [];
                probGNAT(kk).indEQ2n = [];
                probGNAT(kk).indLT2nElse = [];
                
                probGNAT(kk).indLT2npm.top = [];
                probGNAT(kk).indLT2npm.mid = [];
                probGNAT(kk).indLT2npm.bot = [];
                probGNAT(kk).indGTmnp2nsmp1.top = [];
                probGNAT(kk).indGTmnp2nsmp1.mid = [];
                probGNAT(kk).indGTmnp2nsmp1.bot = [];
                probGNAT(kk).indGTmnElse.top = [];
                probGNAT(kk).indGTmnElse.mid = [];
                probGNAT(kk).indGTmnElse.bot = [];
                
                probGNAT(kk).ResInd.indLTn = [];
                probGNAT(kk).ResInd.indEQnp1 = [];
                probGNAT(kk).ResInd.indEQnp2 = [];
                probGNAT(kk).ResInd.indEQ2ns1 = [];
                probGNAT(kk).ResInd.indEQ2n = [];
                probGNAT(kk).ResInd.indLT2nElse = [];
                probGNAT(kk).ResInd.indLT2npmTop = [];
                probGNAT(kk).ResInd.indLT2npmMid = [];
                probGNAT(kk).ResInd.indLT2npmBot = [];
                probGNAT(kk).ResInd.indGTmnp2nsmp1Top = [];
                probGNAT(kk).ResInd.indGTmnp2nsmp1Mid = [];
                probGNAT(kk).ResInd.indGTmnp2nsmp1Bot = [];
                probGNAT(kk).ResInd.indGTmnElseTop = [];
                probGNAT(kk).ResInd.indGTmnElseMid = [];
                probGNAT(kk).ResInd.indGTmnElseBot = [];
                
                probGNAT(kk).JacInd.indLTn = [];
                probGNAT(kk).JacInd.indEQnp1 = [];
                probGNAT(kk).JacInd.indEQnp2 = [];
                probGNAT(kk).JacInd.indEQ2ns1 = [];
                probGNAT(kk).JacInd.indEQ2n = [];
                probGNAT(kk).JacInd.indLT2nElse = [];
                probGNAT(kk).JacInd.indLT2npmTop = [];
                probGNAT(kk).JacInd.indLT2npmMid = [];
                probGNAT(kk).JacInd.indLT2npmBot = [];
                probGNAT(kk).JacInd.indGTmnp2nsmp1Top = [];
                probGNAT(kk).JacInd.indGTmnp2nsmp1Mid = [];
                probGNAT(kk).JacInd.indGTmnp2nsmp1Bot = [];
                probGNAT(kk).JacInd.indGTmnElseTop = [];
                probGNAT(kk).JacInd.indGTmnElseMid = [];
                probGNAT(kk).JacInd.indGTmnElseBot = [];
                
                cnt = 1;
                for i = 1:length(gnat.sampleInd{kk,1})
                    k = length(gnat.irstart{kk,1}(i):gnat.irstart{kk,1}(i+1)-1);
                    if gnat.sampleInd{kk,1}(i) <= obj.N
                        probGNAT(kk).indLTn = [probGNAT(kk).indLTn; gnat.irstart{kk,1}(i)];
                        probGNAT(kk).ResInd.indLTn = [probGNAT(kk).ResInd.indLTn; i];
                        probGNAT(kk).JacInd.indLTn = [probGNAT(kk).JacInd.indLTn; (cnt:cnt+k-1)'];
                    elseif gnat.sampleInd{kk,1}(i) == obj.N+1
                        probGNAT(kk).indEQnp1 = [probGNAT(kk).indEQnp1; gnat.irstart{kk,1}(i)];
                        probGNAT(kk).ResInd.indEQnp1 = [probGNAT(kk).ResInd.indEQnp1; i];
                        probGNAT(kk).JacInd.indEQnp1 = [probGNAT(kk).JacInd.indEQnp1; (cnt:cnt+k-1)'];
                    elseif gnat.sampleInd{kk,1}(i) == obj.N+2
                        probGNAT(kk).indEQnp2 = [probGNAT(kk).indEQnp2; gnat.irstart{kk,1}(i)];
                        probGNAT(kk).ResInd.indEQnp2 = [probGNAT(kk).ResInd.indEQnp2; i];
                        probGNAT(kk).JacInd.indEQnp2 = [probGNAT(kk).JacInd.indEQnp2; (cnt:cnt+k-1)'];
                    elseif gnat.sampleInd{kk,1}(i) == 2*obj.N-1
                        probGNAT(kk).indEQ2ns1 = [probGNAT(kk).indEQ2ns1; gnat.irstart{kk,1}(i)];
                        probGNAT(kk).ResInd.indEQ2ns1 = [probGNAT(kk).ResInd.indEQ2ns1; i];
                        probGNAT(kk).JacInd.indEQ2ns1 = [probGNAT(kk).JacInd.indEQ2ns1; (cnt:cnt+k-1)'];
                    elseif gnat.sampleInd{kk,1}(i) == 2*obj.N
                        probGNAT(kk).indEQ2n = [probGNAT(kk).indEQ2n; gnat.irstart{kk,1}(i)];
                        probGNAT(kk).ResInd.indEQ2n = [probGNAT(kk).ResInd.indEQ2n; i];
                        probGNAT(kk).JacInd.indEQ2n = [probGNAT(kk).JacInd.indEQ2n; (cnt:cnt+k-1)'];
                    elseif gnat.sampleInd{kk,1}(i) <= 2*obj.N
                        probGNAT(kk).indLT2nElse = [probGNAT(kk).indLT2nElse; gnat.irstart{kk,1}(i)];
                        probGNAT(kk).ResInd.indLT2nElse = [probGNAT(kk).ResInd.indLT2nElse; i];
                        probGNAT(kk).JacInd.indLT2nElse = [probGNAT(kk).JacInd.indLT2nElse; (cnt:cnt+k-1)'];
                    elseif gnat.sampleInd{kk,1}(i) <= 2*obj.N + obj.M
                        if gnat.sampleInd{kk,1}(i) == 2*obj.N + 1
                            probGNAT(kk).indLT2npm.top = [probGNAT(kk).indLT2npm.top; gnat.irstart{kk,1}(i)];
                            probGNAT(kk).ResInd.indLT2npmTop = [probGNAT(kk).ResInd.indLT2npmTop; i];
                            probGNAT(kk).JacInd.indLT2npmTop = [probGNAT(kk).JacInd.indLT2npmTop; (cnt:cnt+k-1)'];
                        elseif gnat.sampleInd{kk,1}(i) == 2*obj.N + obj.M
                            probGNAT(kk).indLT2npm.bot = [probGNAT(kk).indLT2npm.bot; gnat.irstart{kk,1}(i)];
                            probGNAT(kk).ResInd.indLT2npmBot = [probGNAT(kk).ResInd.indLT2npmBot; i];
                            probGNAT(kk).JacInd.indLT2npmBot = [probGNAT(kk).JacInd.indLT2npmBot; (cnt:cnt+k-1)'];
                        else
                            probGNAT(kk).indLT2npm.mid = [probGNAT(kk).indLT2npm.mid; gnat.irstart{kk,1}(i)];
                            probGNAT(kk).ResInd.indLT2npmMid = [probGNAT(kk).ResInd.indLT2npmMid; i];
                            probGNAT(kk).JacInd.indLT2npmMid = [probGNAT(kk).JacInd.indLT2npmMid; (cnt:cnt+k-1)'];
                        end
                    elseif gnat.sampleInd{kk,1}(i) >= obj.M*obj.N +2*obj.N - obj.M + 1
                        if gnat.sampleInd{kk,1}(i) == obj.M*obj.N +2*obj.N - obj.M + 1
                            probGNAT(kk).indGTmnp2nsmp1.top = [probGNAT(kk).indGTmnp2nsmp1.top; gnat.irstart{kk,1}(i)];
                            probGNAT(kk).ResInd.indGTmnp2nsmp1Top = [probGNAT(kk).ResInd.indGTmnp2nsmp1Top; i];
                            probGNAT(kk).JacInd.indGTmnp2nsmp1Top = [probGNAT(kk).JacInd.indGTmnp2nsmp1Top; (cnt:cnt+k-1)'];
                        elseif gnat.sampleInd{kk,1}(i) == obj.M*obj.N +2*obj.N
                            probGNAT(kk).indGTmnp2nsmp1.bot = [probGNAT(kk).indGTmnp2nsmp1.bot; gnat.irstart{kk,1}(i)];
                            probGNAT(kk).ResInd.indGTmnp2nsmp1Bot = [probGNAT(kk).ResInd.indGTmnp2nsmp1Bot; i];
                            probGNAT(kk).JacInd.indGTmnp2nsmp1Bot = [probGNAT(kk).JacInd.indGTmnp2nsmp1Bot; (cnt:cnt+k-1)'];
                        else
                            probGNAT(kk).indGTmnp2nsmp1.mid = [probGNAT(kk).indGTmnp2nsmp1.mid; gnat.irstart{kk,1}(i)];
                            probGNAT(kk).ResInd.indGTmnp2nsmp1Mid = [probGNAT(kk).ResInd.indGTmnp2nsmp1Mid; i];
                            probGNAT(kk).JacInd.indGTmnp2nsmp1Mid = [probGNAT(kk).JacInd.indGTmnp2nsmp1Mid; (cnt:cnt+k-1)'];
                        end
                    else
                        if ~isempty(intersect(gnat.sampleInd{kk,1}(i),2*obj.N+1:obj.M:2*obj.N+obj.M*obj.N))
                            probGNAT(kk).indGTmnElse.top = [probGNAT(kk).indGTmnElse.top; gnat.irstart{kk,1}(i)];
                            probGNAT(kk).ResInd.indGTmnElseTop = [probGNAT(kk).ResInd.indGTmnElseTop; i];
                            probGNAT(kk).JacInd.indGTmnElseTop = [probGNAT(kk).JacInd.indGTmnElseTop; (cnt:cnt+k-1)'];
                        elseif ~isempty(intersect(gnat.sampleInd{kk,1}(i),2*obj.N+obj.M:obj.M:2*obj.N+obj.M*obj.N))
                            probGNAT(kk).indGTmnElse.bot = [probGNAT(kk).indGTmnElse.bot; gnat.irstart{kk,1}(i)];
                            probGNAT(kk).ResInd.indGTmnElseBot = [probGNAT(kk).ResInd.indGTmnElseBot; i];
                            probGNAT(kk).JacInd.indGTmnElseBot = [probGNAT(kk).JacInd.indGTmnElseBot; (cnt:cnt+k-1)'];
                        else
                            probGNAT(kk).indGTmnElse.mid = [probGNAT(kk).indGTmnElse.mid; gnat.irstart{kk,1}(i)];
                            probGNAT(kk).ResInd.indGTmnElseMid = [probGNAT(kk).ResInd.indGTmnElseMid; i];
                            probGNAT(kk).JacInd.indGTmnElseMid = [probGNAT(kk).JacInd.indGTmnElseMid; (cnt:cnt+k-1)'];
                        end
                    end
                    cnt = cnt + k;
                end
            end
        end
        
        function   [axOUT] = problemPaperPlot(obj,modelobj,modelAuxobj,pstr,axIN)
            %This function generates the MEMS plot from Rewienski thesis.
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
            
            %Extract reconstructed state vectors (of the beam displacement)
            %with DBCs and time vector 
            %[~,~,svF] = obj.returnPosSV(modelobj,modelAuxobj,1,'sv',1);
            if mod(obj.N,2) == 0
                svF=modelobj.sv(obj.N/2,:);
            else
                svF=0.5*(modelobj.sv(floor(obj.N/2),:) + modelobj.sv(ceil(obj.N/2),:));
            end
            tvec = modelobj.time.T(1):modelobj.time.dt:modelobj.time.T(2);
            
            %Set up axes
            if isempty(axIN)
                ax = axes;
                title('Center Point Deflection Time History','interpreter','latex'); hold on;
                xlabel('Time [ms]','fontsize',14,'interpreter','latex');
                ylabel('Center Point Deflection $[\mu m]$','interpreter','latex');
            else
                ax = axIN;
            end
            
            plot(ax,tvec,svF,pstr{:});
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
            %              vector to extract (1 -> beam displacement, 2 ->
            %              d(u^3)/dt, 3 -> pressure).
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
            
            x1 = svF(1:obj.N,:); x2 = svF(obj.N+1:2*obj.N,:);
            x3 = svF(2*obj.N+1:end,:);
            
            n = size(svF,2);
            
            switch svid
                case 1
                    %Extract mesh spacing
                    X = obj.coord.x;
                    Y = [];
                    
                    %Include DBC in state
                    svDBC = [repmat(obj.z0,1,n);x1;repmat(obj.z0,1,n)];
                case 2
                    %Extract mesh spacing
                    X = obj.coord.x;
                    Y = [];
                    
                    %Include DBC in state
                    svDBC = [zeros(1,n);x2;zeros(1,n)];
                case 3
                    %Extract mesh spacing
                    [X,Y] = meshgrid(obj.coord.x,obj.coord.y);
                    
                    %Include DBC in state
                    x3 = reshape(x3,obj.M,obj.N,n);
                    x3DBCy = repmat(obj.p0,[1,obj.N+2,n]);
                    x3DBCx1 = x3(:,1,:);
                    x3DBCx2 = x3(:,end,:);
                    
                    x3 = [x3DBCx1, x3, x3DBCx2];
                    x3 = [x3DBCy; x3; x3DBCy];
                    
                    svDBC = x3;
                otherwise
                    error('svid must be 1, 2, or 3 for the MEMS problem');
            end
        end
        
        function   [] = problemAnimate(obj,modelobj,modelAuxobj,pstr,type,freq)
            %This function generates a MEMS animation.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj         - Problem object
            %modelobj    - FOM, ROM, GNAT, or TPWL object
            %modelAuxobj - ROM object if modelobj is of class GNAT
            %pstr        - cell array containing the plotting options
            %type        - scalar indicating which animation to use (1 ->
            %              beamdisp, 2 -> d(u^3)/dt, 3 -> pressure, 12 -> 1
            %              and 2, 13 -> 1 and 3, 23 -> 2 and 3.
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
            
%             %Extract reconstructed state vectors
%             svF = modelobj.reconstructFullState(modelAuxobj);
%             x1 = svF(1:obj.N,:); x2 = svF(obj.N+1:2*obj.N,:);
%             x3 = svF(2*obj.N+1:end,:);
            
            %Extract mesh spacing
            x = obj.coord.x;
            y = obj.coord.y;
            [X,Y] = meshgrid(x,y);
            
            %Establish time vector
            tvec = modelobj.time.T(1):modelobj.time.dt:modelobj.time.T(2);
            
            set(gcf,'position',[439   566   714   540]);
            switch type
                case 1 %Beam Displacement
                    ax = axes;
                    [~,~,x1] = obj.returnPosSV(modelobj,modelAuxobj,1);
                    %Determine plot limits
                    xmin = min(x); xmax = max(x);
                    ymin = min(min(x1)); ymax = max(max(x1));
                    
                    %Include DBC in plot
%                     x1 = [repmat(obj.z0,1,modelobj.time.nstep+1);x1;repmat(obj.z0,1,modelobj.time.nstep+1)];
                    
                    title(ax,'Beam Deflection as a Function of Time','fontsize',16,'fontname','georgia');
                    xlabel(ax,'Position','fontsize',14,'fontname','georgia');
                    ylabel(ax,'Beam Deflection','fontsize',14,'fontname','georgia');
                    set(ax,'nextplot','replacechildren','ylim',[ymin,ymax],'xlim',[xmin,xmax]);
                    for i = 1:freq:length(tvec)
                        plot(ax,x,x1(:,i),pstr{:}); pause(0.001);
                    end
                case 2 %d(u^3)/dt
                    ax = axes;
                    [~,~,x2] = obj.returnPosSV(modelobj,modelAuxobj,2);
                    %Determine plot limits
                    xmin = min(x); xmax = max(x);
                    ymin = min(min(x2)); ymax = max(max(x2));
                    
                    %Include DBC in plot
%                     x2 = [zeros(1,modelobj.time.nstep+1);x2;zeros(1,modelobj.time.nstep+1)];
                    
                    title(ax,'Time Derivative of Beam Deflection Cubed as a Function of Time','fontsize',16,'fontname','georgia');
                    xlabel(ax,'Position','fontsize',14,'fontname','georgia');
                    ylabel(ax,'Time Derivative of Beam Deflection Cubed','fontsize',14,'fontname','georgia');
                    set(ax,'nextplot','replacechildren','ylim',[ymin,ymax],'xlim',[xmin,xmax]);
                    for i = 1:freq:length(tvec)
                        plot(ax,x,x2(:,i),pstr{:}); pause(0.05);
                    end
                case 3 %pressure
                    ax = axes;
                    [~,~,x3] = obj.returnPosSV(modelobj,modelAuxobj,3);
                    %Determine plot limits
                    xmin = min(x); xmax = max(x);
                    ymin = min(y); ymax = max(y);
                    zmin = min(min(x3)); zmax = max(max(x3));
                    
%                     x3 = reshape(x3,obj.M,obj.N,modelobj.time.nstep+1);
%                     x3DBCy = repmat(obj.p0,[1,obj.N+2,modelobj.time.nstep+1]);
%                     x3DBCx1 = x3(:,1,:);
%                     x3DBCx2 = x3(:,end,:);
%                     
%                     x3 = [x3DBCx1, x3, x3DBCx2];
%                     x3 = [x3DBCy; x3; x3DBCy];

                    title(ax,'Pressure as a Function of Time','fontsize',16,'fontname','georgia');
                    xlabel(ax,'Position','fontsize',14,'fontname','georgia');
                    ylabel(ax,'Pressure (kPa)','fontsize',14,'fontname','georgia');
                    set(ax,'nextplot','replacechildren','ylim',[ymin,ymax],'xlim',[xmin,xmax],'zlim',[zmin,zmax],'clim',[zmin,zmax]);
                    for i = 1:freq:length(tvec)
                        surf(ax,X,Y,x3(:,:,i)); pause(0.001);
                    end
                case 12 %Beam Displacement and d(u^3)/dt
                    [~,~,x1] = obj.returnPosSV(modelobj,modelAuxobj,1);
                    [~,~,x2] = obj.returnPosSV(modelobj,modelAuxobj,2);
                    
                    %Determine plot limits
                    xmin1 = min(x); xmax1 = max(x);
                    ymin1 = min(min(x1)); ymax1 = max(max(x1));
                    
                    xmin2 = min(x); xmax2 = max(x);
                    ymin2 = min(min(x2)); ymax2 = max(max(x2));
                    
                    %Include DBC in plot
%                     x1 = [repmat(obj.z0,1,modelobj.time.nstep+1);x1;repmat(obj.z0,1,modelobj.time.nstep+1)];
%                     x2 =
%                     [zeros(1,modelobj.time.nstep+1);x2;zeros(1,modelobj.time.nstep+1)];
                    
                    ax(1) = subplot(2,1,1);
                    ax(2) = subplot(2,1,2);
                    
                    title(ax(1),'Beam Deflection as a Function of Time','fontsize',16,'fontname','georgia');
                    xlabel(ax(1),'Position','fontsize',14,'fontname','georgia');
                    ylabel(ax(1),'Beam Deflection','fontsize',14,'fontname','georgia');
                    set(ax(1),'nextplot','replacechildren','ylim',[ymin1,ymax1],'xlim',[xmin1,xmax1]);
                    
                    title(ax(2),'Time Derivative of Beam Deflection Cubed as a Function of Time','fontsize',16,'fontname','georgia');
                    xlabel(ax(2),'Position','fontsize',14,'fontname','georgia');
                    ylabel(ax(2),'Time Derivative of Beam Deflection Cubed','fontsize',14,'fontname','georgia');
                    set(ax(2),'nextplot','replacechildren','ylim',[ymin2,ymax2],'xlim',[xmin2,xmax2]);

                    for i = 1:freq:length(tvec)
                        plot(ax,x,x1(:,i),pstr{:});
                        plot(ax,x,x2(:,i),pstr{:}); pause(0.05);
                    end
                case 13 %Beam Displacement and pressure
                    [~,~,x1] = obj.returnPosSV(modelobj,modelAuxobj,1);
                    [~,~,x3] = obj.returnPosSV(modelobj,modelAuxobj,3);
                    %Determine plot limits
                    xmin1 = min(x); xmax1 = max(x);
                    ymin1 = min(min(x1)); ymax1 = max(max(x1));
                    
                    xmin3 = min(x); xmax3 = max(x);
                    ymin3 = min(y); ymax3 = max(y);
                    zmin3 = min(min(x3)); zmax3 = max(max(x3));
                    
                    %Include DBC in plot
%                     x1 = [repmat(obj.z0,1,modelobj.time.nstep+1);x1;repmat(obj.z0,1,modelobj.time.nstep+1)];
%                     x3 = reshape(x3,obj.M,obj.N,modelobj.time.nstep+1);
%                     x3DBCy = repmat(obj.p0,[1,obj.N+2,modelobj.time.nstep+1]);
%                     x3DBCx1 = x3(:,1,:);
%                     x3DBCx2 = x3(:,end,:);
%                     x3 = [x3DBCx1, x3, x3DBCx2];
%                     x3 = [x3DBCy; x3; x3DBCy];
                    
                    ax(1) = subplot(2,1,1);
                    ax(2) = subplot(2,1,2);
                    
                    title(ax(1),'Beam Deflection as a Function of Time','fontsize',16,'fontname','georgia');
                    xlabel(ax(1),'Position','fontsize',14,'fontname','georgia');
                    ylabel(ax(1),'Beam Deflection','fontsize',14,'fontname','georgia');
                    set(ax(1),'nextplot','replacechildren','ylim',[ymin1,ymax1],'xlim',[xmin1,xmax1]);
                    
                    title(ax(2),'Pressure as a Function of Time','fontsize',16,'fontname','georgia');
                    xlabel(ax(2),'Position','fontsize',14,'fontname','georgia');
                    ylabel(ax(2),'Pressure (kPa)','fontsize',14,'fontname','georgia');
                    set(ax(2),'nextplot','replacechildren','ylim',[ymin3,ymax3],'xlim',[xmin3,xmax3],'zlim',[zmin3,zmax3],'clim',[zmin3,zmax3]);

                    for i = 1:freq:length(tvec)
                        plot(ax,x,x1(:,i),pstr{:});
                        surf(ax,X,Y,x3(:,:,i),pstr{:}); pause(0.001);
                    end
                case 23 %d(u^3)/dt and pressure
                    [~,~,x2] = obj.returnPosSV(modelobj,modelAuxobj,2);
                    [~,~,x3] = obj.returnPosSV(modelobj,modelAuxobj,3);
                    
                    %Determine plot limits
                    xmin2 = min(x); xmax2 = max(x);
                    ymin2 = min(min(x2)); ymax2 = max(max(x2));
                    
                    xmin3 = min(x); xmax3 = max(x);
                    ymin3 = min(y); ymax3 = max(y);
                    zmin3 = min(min(x3)); zmax3 = max(max(x3));
                    
                    ax(1) = subplot(2,1,1);
                    ax(2) = subplot(2,1,2);
                    
                    title(ax(1),'Time Derivative of Beam Deflection Cubed as a Function of Time','fontsize',16,'fontname','georgia');
                    xlabel(ax(1),'Position','fontsize',14,'fontname','georgia');
                    ylabel(ax(1),'Time Derivative of Beam Deflection Cubed','fontsize',14,'fontname','georgia');
                    set(ax(1),'nextplot','replacechildren','ylim',[ymin2,ymax2],'xlim',[xmin2,xmax2]);
                    
                    title(ax(2),'Pressure as a Function of Time','fontsize',16,'fontname','georgia');
                    xlabel(ax(2),'Position','fontsize',14,'fontname','georgia');
                    ylabel(ax(2),'Pressure (kPa)','fontsize',14,'fontname','georgia');
                    set(ax(2),'nextplot','replacechildren','ylim',[ymin3,ymax3],'xlim',[xmin3,xmax3],'zlim',[zmin3,zmax3],'clim',[zmin3,zmax3]);
                    
                    for i = 1:freq:length(tvec)
                        plot(ax,x,x2(:,i),pstr{:}); pause(0.001);
                        surf(ax,X,Y,x3(:,:,i)); pause(0.001);
                    end
                case 123 %beam displacement and d(u^3)/dt and pressure
                    [~,~,x1] = obj.returnPosSV(modelobj,modelAuxobj,1);
                    [~,~,x2] = obj.returnPosSV(modelobj,modelAuxobj,2);
                    [~,~,x3] = obj.returnPosSV(modelobj,modelAuxobj,3);
                    
                    %Determine plot limits
                    xmin1 = min(x); xmax1 = max(x);
                    ymin1 = min(min(x1)); ymax1 = max(max(x1));
                    
                    xmin2 = min(x); xmax2 = max(x);
                    ymin2 = min(min(x2)); ymax2 = max(max(x2));
                    
                    xmin3 = min(x); xmax3 = max(x);
                    ymin3 = min(y); ymax3 = max(y);
                    zmin3 = min(min(x3)); zmax3 = max(max(x3));
                    
                    %Setup axes
                    ax(1) = subplot(3,1,1);
                    ax(2) = subplot(3,1,2);
                    ax(3) = subplot(3,1,3);
                    
                    %Setup lablels and axes properties
                    title(ax(1),'Beam Deflection as a Function of Time','fontsize',16,'fontname','georgia');
                    xlabel(ax(1),'Position','fontsize',14,'fontname','georgia');
                    ylabel(ax(1),'Beam Deflection','fontsize',14,'fontname','georgia');
                    set(ax(1),'nextplot','replacechildren','ylim',[ymin1,ymax1],'xlim',[xmin1,xmax1]);
                    
                    title(ax(2),'Time Derivative of Beam Deflection Cubed as a Function of Time','fontsize',16,'fontname','georgia');
                    xlabel(ax(2),'Position','fontsize',14,'fontname','georgia');
                    ylabel(ax(2),'Time Derivative of Beam Deflection Cubed','fontsize',14,'fontname','georgia');
                    set(ax(2),'nextplot','replacechildren','ylim',[ymin2,ymax2],'xlim',[xmin2,xmax2]);
                    
                    title(ax(3),'Pressure as a Function of Time','fontsize',16,'fontname','georgia');
                    xlabel(ax(3),'Position','fontsize',14,'fontname','georgia');
                    ylabel(ax(3),'Pressure (kPa)','fontsize',14,'fontname','georgia');
                    set(ax(3),'nextplot','replacechildren','ylim',[ymin3,ymax3],'xlim',[xmin3,xmax3],'zlim',[zmin3,zmax3],'clim',[zmin3,zmax3]);
                    
                    %Run animation
                    for i = 1:freq:length(tvec)
                        plot(ax,x,x1(:,i),pstr{:});
                        plot(ax,x,x2(:,i),pstr{:});
                        surf(ax,X,Y,x3(:,:,i)); pause(0.001);
                    end
            end
        end
        
        function  [varargout] = setPrevSV(varargin)
        end
    end
end
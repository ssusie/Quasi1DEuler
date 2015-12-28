classdef SteadyNozzle < handle
    
    properties (SetAccess = private, GetAccess = public)
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
        ndof;
        
        AreaParam = 'spline';
        A;
        dAdp;
        gamma;
        U0;
        UL;
        A_Bnd = [1;1];
        
        nSpPts;
        SpPts;
        
        pbar;
        Gamma;
        targetSoln = []; %This is going to be the target solution for parameter estimation
    end
    
    properties (Hidden=true)
        p;
        Atype = 1;
        
        ind1;
        ind2;
        indN;
        
        nsamp;
        dxGNAT;
        sampleInd;
        indOfUniqueJrow;
        reconstJhatInd;
        Jhat;
        
        dA_imH_dp;
        dA_ipH_dp;
    end
    
    methods
        function  [obj] = SteadyNozzle(cfgobj,oldobj)
            %This is the constructor for the 1D Steady Nozzle Eqn.  If
            %there is 1 input, the constructor will create the class
            %instance from the variables defined in this function.  If
            %there are 2 inputs, the first input will be unused and the
            %output will be a deep copy of the oldobj.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %cfgobj - CONFIG object
            %oldobj - steadyNozzle object whose properties will be copied
            %         into the new object.
            %
            %Outputs:
            %--------
            %obj - SteadyNozzle object
            %--------------------------------------------------------------
            
            if nargin == 0
                %Keep everything at defaults if no input arguments
                return;
            end
            
            % Construct a new object based on a deep copy of an old object
            % of this class by copying properties over.
            if nargin == 2 && isa(oldobj,'SteadyNozzle')
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
            obj.mesh.node = linspace(cfgobj.DLim(1,1),cfgobj.DLim(1,2),cfgobj.nNodes(1))';
            obj.mesh.conn = [1:cfgobj.nNodes(1)-1;2:cfgobj.nNodes(1)]';
            obj.mesh.dbc = [1,true;...
                length(obj.mesh.node),true];
            
            %Number of dofs is 2 less than number of nodes b/c of 2 DBCs
            obj.ndof = cfgobj.nNodes(1)-2;
            
            %Determine grid spacing
            obj.dx = diff(obj.mesh.node);
            
            %Number of unknowns per node
            obj.ndim = 1;
            
            %Determine grid spacing
            coord = obj.mesh.node(2:end-1,1); %want coordinate vector to only include nodes that remain after removing dbc
            
            %Extract parameters
            obj.U0 = obj.config.param{1};
            obj.UL = obj.config.param{2};
            obj.gamma = obj.config.param{3};
            obj.AreaParam = obj.config.param{4};
            p = obj.config.param{5};
            p = p(:);
            
            %Setup other properties
            xi = 0.5*(obj.mesh.node(2:end,1) + obj.mesh.node(1:end-1,1));
            
            switch obj.AreaParam
                case 'splineFixedEnds'
                    %Determine number of spline knots to parametrize shape
                    obj.nSpPts = length(p)+2;
                    p = [obj.A_Bnd(1);p;obj.A_Bnd(2)];
                    obj.SpPts = linspace(cfgobj.DLim(1,1),cfgobj.DLim(1,2),obj.nSpPts)';
                    obj.Atype = 1;
                    %%Setup Nozzle shape (cubic spline points)
                    [obj.A,obj.dAdp] = obj.NozzleAreaSpline(xi,p);
                case 'splineFreeEnds'
                    %Determine number of spline knots to parametrize shape
                    obj.nSpPts = length(p);
                    obj.SpPts = linspace(cfgobj.DLim(1,1),cfgobj.DLim(1,2),obj.nSpPts)';
                    obj.Atype = 2;
                    %%Setup Nozzle shape (cubic spline points)
                    [obj.A,obj.dAdp] = obj.NozzleAreaSpline(xi,p);
                case '1D'
                    obj.Atype = 3;
                    %%Setup Nozzle shape (cubic spline points)
                    [obj.A,obj.dAdp] = obj.NozzleArea1D(xi,p);
                case 'AIAA'
                    obj.Atype = 4;
                    [obj.A,obj.dAdp] = obj.NozzleAreaAIAA(xi,p);
            end

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
            
            %Determine output matrix from the output matrix function
            if isa(cfgobj.CFunc,'function_handle')
                obj.C = feval(cfgobj.CFunc);
            else
                obj.C = feval(cfgobj.CFunc);
            end
        end
                
        function [u,dudphi] = velocity(obj,phi)
            n = size(phi,1);
            u = (phi(2:end,1)-phi(1:end-1,1))./obj.dx;
            dudphi = spdiags([-1./obj.dx 1./obj.dx], 0:1,n-1,n); % checked by FD
        end
        
        function  [u_ipH,u_imH,u_im3H] = velocityGNAT(obj,phi)
            
            %du = [-ones(obj.nsamp,1),ones(obj.nsamp,1)];
            u_ipH = zeros(obj.nsamp,1);
            u_imH=zeros(obj.nsamp,1);
            u_im3H=zeros(obj.nsamp,1);
            if obj.ind1 && obj.ind2
                u_ipH(1)=(phi(2)-phi(1))/obj.dxGNAT;
                u_imH(1)=(phi(1)-obj.U0)/obj.dxGNAT;
                u_im3H(1)=0;
                
                u_ipH(2)=(phi(5)-phi(4))/obj.dxGNAT;
                u_imH(2)=(phi(4)-phi(3))/obj.dxGNAT;
                u_im3H(2)=(phi(3)-obj.U0)/obj.dxGNAT;
                
                if obj.indN
                    u_ipH(3:end-1) = (phi(9:4:end-3)-phi(8:4:end-3))./obj.dxGNAT;
                    u_imH(3:end-1) = (phi(8:4:end-3)-phi(7:4:end-3))./obj.dxGNAT;
                    u_im3H(3:end-1) = (phi(7:4:end-3)-phi(6:4:end-3))./obj.dxGNAT;
                    
                    u_ipH(end)=(obj.UL-phi(end))/obj.dxGNAT;
                    u_imH(end)=(phi(end)-phi(end-1))/obj.dxGNAT;
                    u_im3H(end)=(phi(end-1)-phi(end-2))/obj.dxGNAT;
                else
                    u_ipH(3:end) = (phi(9:4:end)-phi(8:4:end))./obj.dxGNAT;
                    u_imH(3:end) = (phi(8:4:end)-phi(7:4:end))./obj.dxGNAT;
                    u_im3H(3:end) = (phi(7:4:end)-phi(6:4:end))./obj.dxGNAT;
                end
            elseif obj.ind1
                u_ipH(1)=(phi(2)-phi(1))/obj.dxGNAT;
                u_imH(1)=(phi(1)-obj.U0)/obj.dxGNAT;
                u_im3H(1)=0;
                
                if obj.indN
                    u_ipH(2:end-1) = (phi(6:4:end-3)-phi(5:4:end-3))./obj.dxGNAT;
                    u_imH(2:end-1) = (phi(5:4:end-3)-phi(4:4:end-3))./obj.dxGNAT;
                    u_im3H(2:end-1) = (phi(4:4:end-3)-phi(3:4:end-3))./obj.dxGNAT;
                    
                    u_ipH(end)=(obj.UL-phi(end))/obj.dxGNAT;
                    u_imH(end)=(phi(end)-phi(end-1))/obj.dxGNAT;
                    u_im3H(end)=(phi(end-1)-phi(end-2))/obj.dxGNAT;
                else
                    u_ipH(2:end) = (phi(6:4:end)-phi(5:4:end))./obj.dxGNAT;
                    u_imH(2:end) = (phi(5:4:end)-phi(4:4:end))./obj.dxGNAT;
                    u_im3H(2:end) = (phi(4:4:end)-phi(3:4:end))./obj.dxGNAT;
                end
            elseif obj.ind2
                u_ipH(1)=(phi(3)-phi(2))/obj.dxGNAT;
                u_imH(1)=(phi(2)-phi(1))/obj.dxGNAT;
                u_im3H(1)=(phi(1)-obj.U0)/obj.dxGNAT;
                
                if obj.indN
                    u_ipH(2:end-1) = (phi(7:4:end-3)-phi(6:4:end-3))./obj.dxGNAT;
                    u_imH(2:end-1) = (phi(6:4:end-3)-phi(5:4:end-3))./obj.dxGNAT;
                    u_im3H(2:end-1) = (phi(5:4:end-3)-phi(4:4:end-3))./obj.dxGNAT;
                    
                    u_ipH(end)=(obj.UL-phi(end))/obj.dxGNAT;
                    u_imH(end)=(phi(end)-phi(end-1))/obj.dxGNAT;
                    u_im3H(end)=(phi(end-1)-phi(end-2))/obj.dxGNAT;
                else
                    u_ipH(2:end) = (phi(7:4:end)-phi(6:4:end))./obj.dxGNAT;
                    u_imH(2:end) = (phi(6:4:end)-phi(5:4:end))./obj.dxGNAT;
                    u_im3H(2:end) = (phi(5:4:end)-phi(4:4:end))./obj.dxGNAT;
                end
            end
        end
        
        function [rho,drhodu] = density(obj,u)
            n = size(u,1);
            rho = (1+(obj.gamma-1)/2*(1-u.^2)).^(1/(obj.gamma-1));
            drhodu = spdiags(-u.*rho.^(2-obj.gamma),0,n,n); % checked by FD
        end
        
        function [rho,drhodu] = densityGNAT(obj,u)
            n = size(u,1);
            rho = (1+(obj.gamma-1)/2*(1-u.^2)).^(1/(obj.gamma-1));
            drhodu = -u.*rho.^(2-obj.gamma); % checked by FD
        end
        
        function [rho_tilde,drho_tildedu,drho_tildedrho] = density_stab(obj,u,rho)
            rho_tilde = zeros(size(rho));
            n = size(rho_tilde,1);
            drho_tildedu = sparse(n,n);
            drho_tildedrho = sparse(n,n);
            
            [z,dzdrho,dzdu] = obj.bias(rho,u);
            
            rho_tilde(1,1) =  rho(1,1) - 1/(u(1,1))*z(1,1);  % subsonic inlet
            drho_tildedu(1,:) = -1/(u(1,1))*dzdu(1,:);
            drho_tildedu(1,1) = drho_tildedu(1,1) + 1/(u(1,1)^2)*z(1,1);
            drho_tildedrho(1,:) = -1/(u(1,1))*dzdrho(1,:);
            drho_tildedrho(1,1) = drho_tildedrho(1,1) + 1;
            
            rho_tilde(2:n,1) = rho(2:n,1) - 1./u(2:n,1).*(z(2:n,1) - z(1:n-1,1));
            drho_tildedu(2:n,2:n) = spdiags(1./(u(2:n,1).^2).*(z(2:n,1) - z(1:n-1,1)),0,n-1,n-1);
            drho_tildedu(2:n,:) =  drho_tildedu(2:n,:) + spdiags(- 1./u(2:n,1),0,n-1,n-1)*(dzdu(2:n,:) -dzdu(1:n-1,:)); % checked by FD
            
            drho_tildedrho(2:n,2:n) = speye(n-1);
            drho_tildedrho(2:n,:) = drho_tildedrho(2:n,:) + spdiags(- 1./u(2:n,1),0,n-1,n-1)*(dzdrho(2:n,:) -dzdrho(1:n-1,:)); % checked by FD
        end
        
        function [rho_tilde,Drho_tildeDu_ipH,Drho_tildeDu_imH] = density_stabGNAT(obj,u_imH,rho_imH,u_ipH,rho_ipH,drhodu_ipH,drhodu_imH)
%         function [rho_tilde,drho_tildedu,drho_tildedrho] = density_stabGNAT(obj,u_imH,rho_imH,u_ipH,rho_ipH,drhodu_ipH)

%             rho_tilde = zeros(size(rho));
%             n = size(rho_tilde,1);
%             drho_tildedu = sparse(n,n);
%             drho_tildedrho = sparse(n,n);
            
            [z_ipH,dzdrho_ipH,dzdu_ipH] = obj.biasGNAT(rho_ipH,u_ipH);
            [z_imH,dzdrho_imH,dzdu_imH] = obj.biasGNAT(rho_imH,u_imH);

            rho_tilde = rho_ipH - (1./u_ipH).*(z_ipH - z_imH);
            
            drho_tildedu = [(1./u_ipH).*dzdu_imH, (1./u_ipH.^2).*(z_ipH-z_imH)-(1./u_ipH).*dzdu_ipH];
            drho_tildedrho = [(1./u_ipH).*dzdrho_imH,1-(1./u_ipH).*dzdrho_ipH];
            
            Drho_tildeDu_ipH = drho_tildedu(:,2) + drho_tildedrho(:,2).*drhodu_ipH;
            Drho_tildeDu_imH = drho_tildedu(:,1) + drho_tildedrho(:,1).*drhodu_imH;
            
%             drho_tildedu(1,:) = -1/(u(1,1))*dzdu(1,:);
%             drho_tildedu(1,1) = drho_tildedu(1,1) + 1/(u(1,1)^2)*z(1,1);
%             drho_tildedrho(1,:) = -1/(u(1,1))*dzdrho(1,:);
%             drho_tildedrho(1,1) = drho_tildedrho(1,1) + 1;
%             
%             drho_tildedu(2:n,2:n) = spdiags(1./(u(2:n,1).^2).*(z(2:n,1) - z(1:n-1,1)),0,n-1,n-1);
%             drho_tildedu(2:n,:) =  drho_tildedu(2:n,:) + spdiags(- 1./u(2:n,1),0,n-1,n-1)*(dzdu(2:n,:) -dzdu(1:n-1,:)); % checked by FD
%             
%             drho_tildedrho(2:n,2:n) = speye(n-1);
%             drho_tildedrho(2:n,:) = drho_tildedrho(2:n,:) + spdiags(- 1./u(2:n,1),0,n-1,n-1)*(dzdrho(2:n,:) -dzdrho(1:n-1,:)); % checked by FD
        end
        
        function [rho_u_bar,drho_u_bardrho,drho_u_bardu] = bias(obj,rho,u)
            rho_star = 1;
            u_star = 1;
            n = size(u,1);
            [M,dMdu,dMdrho] = obj.mach(rho,u);
            % rho_u_bar = 0*(M<1) + (rho.*u - rho_star.*u_star).*(M>=1);
            % smooth version
            k = 10^(12);
            rho_u_bar = (rho.*u - rho_star.*u_star).*(0.5*(1+tanh(k*(M-1))));
            drho_u_bardu = spdiags(rho.*(0.5*(1+tanh(k*(M-1)))) +(rho.*u - rho_star.*u_star).*(0.5*k*(1-(tanh(k*(M-1))).^2).*diag(dMdu)),0,n,n); % checked by FD
            drho_u_bardrho = spdiags(u.*(0.5*(1+tanh(k*(M-1)))) +(rho.*u - rho_star.*u_star).*(0.5*k*(1-(tanh(k*(M-1))).^2).*diag(dMdrho)),0,n,n); % checked by FD
        end

        function [rho_u_bar,drho_u_bardrho,drho_u_bardu] = biasGNAT(obj,rho,u)
            rho_star = 1;
            u_star = 1;
            n = size(u,1);
            [M,dMdu,dMdrho] = obj.machGNAT(rho,u);
            % rho_u_bar = 0*(M<1) + (rho.*u - rho_star.*u_star).*(M>=1);
            % smooth version
            k = 10^(12);
            rho_u_bar = (rho.*u - rho_star.*u_star).*(0.5*(1+tanh(k*(M-1))));
            drho_u_bardu = rho.*(0.5*(1+tanh(k*(M-1)))) +(rho.*u - rho_star.*u_star).*(0.5*k*(1-(tanh(k*(M-1))).^2).*dMdu); % checked by FD
            drho_u_bardrho = u.*(0.5*(1+tanh(k*(M-1)))) +(rho.*u - rho_star.*u_star).*(0.5*k*(1-(tanh(k*(M-1))).^2).*dMdrho); % checked by FD
        end
        
        function [M,dMdu,dMdrho] = mach(obj,rho,u)
            n = size(u,1);
            %M = u./rho.^(obj.gamma-1);
            M = u./rho.^(0.5*(obj.gamma-1));
            dMdu = spdiags(1./rho.^(obj.gamma-1),0,n,n); % checked by FD
            %dMdu = spdiags(1./rho.^(0.5*(obj.gamma-1)),0,n,n); % checked by FD
            dMdrho = spdiags(0.5*(1-obj.gamma)*u.*rho.^(-0.5*(1+obj.gamma)),0,n,n); % checked by FD
            %dMdrho = spdiags((1-obj.gamma)*u.*rho.^(-obj.gamma),0,n,n); % checked by FD
        end
                
        function [M,dMdu,dMdrho] = machGNAT(obj,rho,u)
            n = size(u,1);
            %M = u./rho.^(obj.gamma-1);
            M = u./rho.^(0.5*(obj.gamma-1));
            dMdu = 1./rho.^(obj.gamma-1); % checked by FD
            %dMdu = spdiags(1./rho.^(0.5*(obj.gamma-1)),0,n,n); % checked by FD
            dMdrho = 0.5*(1-obj.gamma)*u.*rho.^(-0.5*(1+obj.gamma)); % checked by FD
            %dMdrho = spdiags((1-obj.gamma)*u.*rho.^(-obj.gamma),0,n,n); % checked by FD
        end

        function [rho,u,M] = getVariables(obj,phi)
            
            u = obj.velocity([obj.U0;phi;obj.UL]);
            rho = obj.density(u);
            M = obj.mach(rho,u);
            
        end
        
        function  [R,J] = ResJac(obj,phi,~)
            %This function calcuates the residual and jacobian of
            %SteadyNozzle
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj - SteadyNozzle object
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
            n = length(obj.A);
            [u,dudphi] = obj.velocity([obj.U0;phi;obj.UL]);
            [rho,drhodu] = obj.density(u);
            [rho_tilde,drho_tildedu,drho_tildedrho] = obj.density_stab(u,rho);
            R = obj.A(2:n,1).*rho_tilde(2:n,1).*u(2:n,1) - obj.A(1:n-1,1).*rho_tilde(1:n-1,1).*u(1:n-1,1);
            drho_tildedphi = drho_tildedu*dudphi + drho_tildedrho*drhodu*dudphi; % checked by FD
            J = spdiags(-obj.A(1:n-1,1).*rho_tilde(1:n-1,1),0,n-1,n-1)*dudphi(1:n-1,:) + spdiags(-obj.A(1:n-1,1).*u(1:n-1,1),0,n-1,n-1)*drho_tildedphi(1:n-1,:) ...
                + spdiags(obj.A(2:n,1).*rho_tilde(2:n,1),0,n-1,n-1)*dudphi(2:n,:) +spdiags(obj.A(2:n,1).*u(2:n,1),0,n-1,n-1)*drho_tildedphi(2:n,:); % checked by FD
            
            J = J(:,2:end-1);
        end
        
        function  [R,J] = ResJacGNAT(obj,phi,~)
            %This function calcuates the residual and jacobian of Burger's
            %Equation
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - SteadyNozzle object
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
            n=size(obj.A,1);
            [u_ipH,u_imH,u_im3H] = obj.velocityGNAT(phi);
            
            [rho_ipH,drhodu_ipH]=obj.densityGNAT(u_ipH);
            [rho_imH,drhodu_imH]=obj.densityGNAT(u_imH);
            [rho_im3H,drhodu_im3H]=obj.densityGNAT(u_im3H);
            
            [rhotilde_ipH,Drhotilde_ipH_Du_ipH,Drhotilde_ipH_Du_imH]=obj.density_stabGNAT(u_imH,rho_imH,u_ipH,rho_ipH,drhodu_ipH,drhodu_imH);
            [rhotilde_imH,Drhotilde_imH_Du_imH,Drhotilde_imH_Du_im3H]=obj.density_stabGNAT(u_im3H,rho_im3H,u_imH,rho_imH,drhodu_imH,drhodu_im3H);
%             [rho_tilde_ipH,drho_tildedu_ipH,drho_tildedrho_ipH]=obj.density_stabGNAT(u_imH,rho_imH,u_ipH,rho_ipH);
%             [rho_tilde_imH,drho_tildedu_imH,drho_tildedrho_imH]=obj.density_stabGNAT(u_im3H,rho_im3H,u_imH,rho_imH);
            
            R = obj.A(:,2).*rhotilde_ipH.*u_ipH - obj.A(:,1).*rhotilde_imH.*u_imH;
            
            J=[];
            if obj.ind1
                Drhotilde_ipH_dphi = Drhotilde_ipH_Du_ipH(1)*[-1/obj.dxGNAT;1/obj.dxGNAT] + Drhotilde_ipH_Du_imH(1)*[1/obj.dxGNAT;0];
                Drhotilde_imH_dphi = Drhotilde_imH_Du_imH(1)*[1/obj.dxGNAT;0];
                
                J = [J; obj.A(1,2)*(u_ipH(1)*Drhotilde_ipH_dphi + rhotilde_ipH(1)*[-1/obj.dxGNAT;1/obj.dxGNAT]) - ...
                        obj.A(1,1)*(u_imH(1)*Drhotilde_imH_dphi + rhotilde_imH(1)*[1/obj.dxGNAT;0])];
            end
            
            if obj.ind2
                if obj.ind1
                    ind = 2;
                else
                    ind = 1;
                end
                
                Drhotilde_ipH_dphi = Drhotilde_ipH_Du_ipH(ind)*[0;-1/obj.dxGNAT;1/obj.dxGNAT] + Drhotilde_ipH_Du_imH(ind)*[-1/obj.dxGNAT;1/obj.dxGNAT;0];
                Drhotilde_imH_dphi = Drhotilde_imH_Du_imH(ind)*[-1/obj.dxGNAT;1/obj.dxGNAT;0] + Drhotilde_imH_Du_im3H(ind)*[1/obj.dxGNAT;0;0];
                
                J = [J; obj.A(ind,2)*(u_ipH(ind)*Drhotilde_ipH_dphi + rhotilde_ipH(ind)*[0;-1/obj.dxGNAT;1/obj.dxGNAT]) - ...
                        obj.A(ind,1)*(u_imH(ind)*Drhotilde_imH_dphi + rhotilde_imH(ind)*[-1/obj.dxGNAT;1/obj.dxGNAT;0])];
            end
            
            if obj.ind1 && obj.ind2, startind=3; elseif obj.ind1||obj.ind2, startind=2; else startind=1; end
            if obj.indN, endind=obj.nsamp-1; else endind=obj.nsamp; end
            numint = endind-startind+1;
            
            du_ipH_dphi = [0,0,-1/obj.dxGNAT,1/obj.dxGNAT];
            du_imH_dphi = [0,-1/obj.dxGNAT,1/obj.dxGNAT,0];
            du_im3H_dphi = [-1/obj.dxGNAT,1/obj.dxGNAT,0,0];
            
            Drhotilde_ipH_dphi = Drhotilde_ipH_Du_ipH(startind:endind)*du_ipH_dphi + Drhotilde_ipH_Du_imH(startind:endind)*du_imH_dphi;
            Drhotilde_imH_dphi = Drhotilde_imH_Du_imH(startind:endind)*du_imH_dphi + Drhotilde_imH_Du_im3H(startind:endind)*du_im3H_dphi;
            
            J = [J; reshape((bsxfun(@times,(bsxfun(@times,Drhotilde_ipH_dphi,u_ipH(startind:endind)) + rhotilde_ipH(startind:endind)*du_ipH_dphi),obj.A(startind:endind,2)) - ...
                bsxfun(@times,(bsxfun(@times,Drhotilde_imH_dphi,u_imH(startind:endind)) + rhotilde_imH(startind:endind)*du_imH_dphi),obj.A(startind:endind,1)))',numint*4,1)];
            
            if obj.indN
                Drhotilde_ipH_dphi = Drhotilde_ipH_Du_ipH(end)*[0;0;-1/obj.dxGNAT] + Drhotilde_ipH_Du_imH(end)*[0;-1/obj.dxGNAT;1/obj.dxGNAT];
                Drhotilde_imH_dphi = Drhotilde_imH_Du_imH(end)*[0;-1/obj.dxGNAT;1/obj.dxGNAT] + Drhotilde_imH_Du_im3H(end)*[-1/obj.dxGNAT;1/obj.dxGNAT;0];
                
                J = [J; obj.A(end,2)*(u_ipH(end)*Drhotilde_ipH_dphi + rhotilde_ipH(end)*[0;0;-1/obj.dxGNAT]) - ...
                        obj.A(end,1)*(u_imH(end)*Drhotilde_imH_dphi + rhotilde_imH(end)*[0;-1/obj.dxGNAT;1/obj.dxGNAT])];
            end
        end
        
        function  [A,b] = linIneqROM(obj)
            
            n = obj.ndof;
            
            %FOR ROM, need A*obj.phi
            A = spdiags(repmat([-1 1],n,1),[0 1],n-1,n);
            b = obj.dx(2:end-1)*sqrt((1+obj.gamma)/(obj.gamma-1)) - A*obj.ic;
        end
        
        function  [A,b] = linIneqGNAT(obj)
            
            n = obj.ndof;
            
            %FOR GNAT, need A*obj.phi
            A = spdiags(repmat([-1 1],n,1),[0 1],n-1,n);
            b = obj.dxGNAT*sqrt((1+obj.gamma)/(obj.gamma-1))*ones(size(A,1),1);
        end
        
        function  [ind] = node2ind(obj,nodenum)
            %This function takes a node in the mesh and maps it to an index
            %in the full state vector
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - SteadyNozzle object
            %nodenum - scalar specifying node numbers
            %
            %Outputs:
            %--------
            %ind     - vector specifying the corresponding indices in the
            %          full state vector
            %--------------------------------------------------------------
            
            if nodenum == 1 || nodenum == length(obj.mesh.node)
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
            %[~,J] = ResJac(obj,rand(N,1),1);
            %obj.Jstruct = J ~= 0;
            
            temp = diag(ones(N,1),0) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1) + diag(ones(N-2,1),-2);
            
            %temp = diag(ones(N,1),0) + diag(ones(N-1,1),-1);
            obj.Jstruct = sparse(temp);
            
        end
        
        function  [probGNAT] = createCopy4GNAT(obj,gnat)
            %This function makes a copy of the current object and only
            %copies the source term and input matrix evaluated at indN
            %rows.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - SteadyNozzle object
            %gnat    - gnat object
            %
            %Outputs:
            %--------
            %probGNAT - new "reduced" OneDBurgers object
            %--------------------------------------------------------------
            
            probGNAT = SteadyNozzle([],obj);
            
            probGNAT.Atype = obj.Atype;
            temp = probGNAT.config.inFunc;
            temp2= probGNAT.config.form;
            probGNAT.config = [];
            probGNAT.config.inFunc = temp;
            probGNAT.config.form = temp2;
%             probGNAT.mesh   = [];
            probGNAT.ndim = [];
            probGNAT.Jstruct = [];
            probGNAT.dx = [];
            probGNAT.A = [];
            probGNAT.dAdp = [];
            
            probGNAT.dA_imH_dp = obj.dAdp(gnat.sampleInd,:);
            probGNAT.dA_ipH_dp = obj.dAdp(gnat.sampleInd+1,:);
            
            probGNAT.A(:,1) = obj.A(gnat.sampleInd,1);
            probGNAT.A(:,2) = obj.A(gnat.sampleInd+1,1);
            probGNAT.dAdp = obj.dAdp(gnat.sampleInd,1);
            
            [~,probGNAT.indOfUniqueJrow] = unique(gnat.jrow);
            probGNAT.sampleInd = gnat.sampleInd;
            probGNAT.reconstJhatInd = gnat.reconstJhatInd;
            probGNAT.Jhat = spalloc(gnat.nI,length(unique(gnat.jrow)),length(gnat.reconstJhatInd));
            
            probGNAT.nsamp=size(gnat.sampleInd,1);
            
            if ~isempty(intersect(gnat.sampleInd,1))
                probGNAT.ind1 = true;
            else
                probGNAT.ind1 = false;
            end
            
            if ~isempty(intersect(gnat.sampleInd,2))
                probGNAT.ind2 = true;
            else
                probGNAT.ind2 = false;
            end
            
            if ~isempty(intersect(gnat.sampleInd,obj.ndof))
                probGNAT.indN = true;
            else
                probGNAT.indN = false;
            end
     
            h = obj.mesh.node(2:end) - obj.mesh.node(1:end-1);
            probGNAT.dxGNAT = h(1);
                        
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
            %obj     - SteadyNozzle object
            %gnat    - gnat object
            %
            %Outputs:
            %--------
            %probGNAT - new "reduced" OneDBurgers object
            %--------------------------------------------------------------
            for kk = 1:gnat.nBases
                if isempty(gnat.sampleInd{kk,1})
                    continue;
                end
                probGNAT(kk) = SteadyNozzle([],obj);
                
                temp = probGNAT(kk).config.inFunc;
                temp2= probGNAT(kk).config.form;
                probGNAT(kk).config = [];
                probGNAT(kk).config.inFunc = temp;
                probGNAT(kk).config.form = temp2;
                probGNAT(kk).mesh   = [];
                probGNAT(kk).ndim = [];
                probGNAT(kk).Jstruct = [];
                probGNAT(kk).dx = [];
                probGNAT(kk).A = [];
                probGNAT(kk).dAdp = [];
                
                probGNAT(kk).dA_imH_dp = obj.dAdp(gnat.sampleInd{kk},:);
                probGNAT(kk).dA_ipH_dp = obj.dAdp(gnat.sampleInd{kk}+1,:);
                
                probGNAT(kk).A(:,1) = obj.A(gnat.sampleInd{kk},1);
                probGNAT(kk).A(:,2) = obj.A(gnat.sampleInd{kk}+1,1);
                
                probGNAT(kk).nsamp=size(gnat.sampleInd{kk},1);
                
                if ~isempty(intersect(gnat.sampleInd{kk},1))
                    probGNAT(kk).ind1 = true;
                else
                    probGNAT(kk).ind1 = false;
                end
                
                if ~isempty(intersect(gnat.sampleInd{kk},2))
                    probGNAT(kk).ind2 = true;
                else
                    probGNAT(kk).ind2 = false;
                end
                
                if ~isempty(intersect(gnat.sampleInd{kk},obj.ndof))
                    probGNAT(kk).indN = true;
                else
                    probGNAT(kk).indN = false;
                end
                
                h = obj.mesh.node(2:end) - obj.mesh.node(1:end-1);
                probGNAT(kk).dxGNAT = h(1);
                
                if length(obj.B) > 1
                    probGNAT(kk).B = obj.B(gnat.sampleInd{kk},1);
                end
                
                if length(obj.G) > 1
                    probGNAT(kk).G = obj.G(gnat.sampleInd{kk},1);
                end
                
                if length(obj.ic) > 1
                    probGNAT(kk).ic = obj.ic(unique(gnat.jrow{kk}),1);
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
        
        function  [CFL] = computeCFL(obj,phi)
            
            u=obj.velocity([obj.U0;phi;obj.UL]);
            CFL=min(obj.dx./u);
            
        end
        
        function  [] = setPrevSV(obj,~)
            return;
        end
        
        %%%%%%%%%%%%%%%Optimization functions%%%%%%%%%%%%%%%%%
        function  [R,dRdw,dRdp,d2RdwdpR,d2Rdw2R] = ResSens(obj,phi,~)
            n = length(obj.A);
            [u,dudphi] = obj.velocity([obj.U0;phi;obj.UL]);
            [rho,drhodu] = obj.density(u);
            [rho_tilde,drho_tildedu,drho_tildedrho] = obj.density_stab(u,rho);
            drho_tildedphi = drho_tildedu*dudphi + drho_tildedrho*drhodu*dudphi; % checked by FD
            R = obj.A(2:n,1).*rho_tilde(2:n,1).*u(2:n,1) - obj.A(1:n-1,1).*rho_tilde(1:n-1,1).*u(1:n-1,1);
            dRdw = spdiags(-obj.A(1:n-1,1).*rho_tilde(1:n-1,1),0,n-1,n-1)*dudphi(1:n-1,:) + spdiags(-obj.A(1:n-1,1).*u(1:n-1,1),0,n-1,n-1)*drho_tildedphi(1:n-1,:) ...
                + spdiags(obj.A(2:n,1).*rho_tilde(2:n,1),0,n-1,n-1)*dudphi(2:n,:) +spdiags(obj.A(2:n,1).*u(2:n,1),0,n-1,n-1)*drho_tildedphi(2:n,:); % checked by FD
            
            dRdw = dRdw(:,2:end-1);
            dRdA = spdiags(rho_tilde(2:n,1).*u(2:n,1),1,n-1,n) - spdiags(rho_tilde(1:n-1,1).*u(1:n-1,1),0,n-1,n); % checked by FD
            dRdp = dRdA*obj.dAdp;
            
            d2Rdw2R = [];
            d2RdwdpR = [];
            if (nargout == 4)
                h = 1e-6;
                d2resdAdphires = sparse(n,n+1);
                d2Rdw2R = zeros(n-1,n-1);
                for i=1:n-1
                    M = [-rho_tilde(i,1)*dudphi(i,:) - u(i,1)*drho_tildedphi(i,:), rho_tilde(i+1,1)*dudphi(i+1,:) + u(i+1,1)*drho_tildedphi(i+1,:)];
                    I = [i*ones(1,n+1),(i+1)*ones(1,n+1)];
                    J = [1:(n+1),1:(n+1)];
                    dresdphidA = sparse(I,J,M,n,n+1);
                    d2RdwdphiR = d2resdAdphires + dresdphidA*R(i,1);
                    
                    e = zeros(size(phi));
                    e(i) = 1;
                    [~,Jp1] = obj.ResJac(phi + h*e,[]);
                    [~,Jm1] = obj.ResJac(phi - h*e,[]);
                    
                    d2Rdw2R(:,i) = (Jp1' - Jm1')*(R/(2*h));
                end
                d2RdwdpR = d2RdwdphiR(:,2:end-1)'*obj.dAdp;
            end
            %else
            %    d2resdAdphires = 0;
            %end
        end
        
        function  [R,dRdw,dRdp,d2RdwdpR,d2Rdw2R] = ResSensGNAT(obj,phi,~)
            %This function calcuates the residual and jacobian of Burger's
            %Equation
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - SteadyNozzle object
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
            n=size(obj.A,1);
            [u_ipH,u_imH,u_im3H] = obj.velocityGNAT(phi);
            
            [rho_ipH,drhodu_ipH]=obj.densityGNAT(u_ipH);
            [rho_imH,drhodu_imH]=obj.densityGNAT(u_imH);
            [rho_im3H,drhodu_im3H]=obj.densityGNAT(u_im3H);
            
            [rhotilde_ipH,Drhotilde_ipH_Du_ipH,Drhotilde_ipH_Du_imH]=obj.density_stabGNAT(u_imH,rho_imH,u_ipH,rho_ipH,drhodu_ipH,drhodu_imH);
            [rhotilde_imH,Drhotilde_imH_Du_imH,Drhotilde_imH_Du_im3H]=obj.density_stabGNAT(u_im3H,rho_im3H,u_imH,rho_imH,drhodu_imH,drhodu_im3H);
            
            R = obj.A(:,2).*rhotilde_ipH.*u_ipH - obj.A(:,1).*rhotilde_imH.*u_imH;
            
            dRdw=[];
            if obj.ind1
                Drhotilde_ipH_dphi = Drhotilde_ipH_Du_ipH(1)*[-1/obj.dxGNAT;1/obj.dxGNAT] + Drhotilde_ipH_Du_imH(1)*[1/obj.dxGNAT;0];
                Drhotilde_imH_dphi = Drhotilde_imH_Du_imH(1)*[1/obj.dxGNAT;0];
                
                dRdw = [dRdw; obj.A(1,2)*(u_ipH(1)*Drhotilde_ipH_dphi + rhotilde_ipH(1)*[-1/obj.dxGNAT;1/obj.dxGNAT]) - ...
                        obj.A(1,1)*(u_imH(1)*Drhotilde_imH_dphi + rhotilde_imH(1)*[1/obj.dxGNAT;0])];
            end
            
            if obj.ind2
                if obj.ind1
                    ind = 2;
                else
                    ind = 1;
                end
                
                Drhotilde_ipH_dphi = Drhotilde_ipH_Du_ipH(ind)*[0;-1/obj.dxGNAT;1/obj.dxGNAT] + Drhotilde_ipH_Du_imH(ind)*[-1/obj.dxGNAT;1/obj.dxGNAT;0];
                Drhotilde_imH_dphi = Drhotilde_imH_Du_imH(ind)*[-1/obj.dxGNAT;1/obj.dxGNAT;0] + Drhotilde_imH_Du_im3H(ind)*[1/obj.dxGNAT;0;0];
                
                dRdw = [dRdw; obj.A(ind,2)*(u_ipH(ind)*Drhotilde_ipH_dphi + rhotilde_ipH(ind)*[0;-1/obj.dxGNAT;1/obj.dxGNAT]) - ...
                        obj.A(ind,1)*(u_imH(ind)*Drhotilde_imH_dphi + rhotilde_imH(ind)*[-1/obj.dxGNAT;1/obj.dxGNAT;0])];
            end
            
            if obj.ind1 && obj.ind2, startind=3; elseif obj.ind1||obj.ind2, startind=2; else startind=1; end
            if obj.indN, endind=obj.nsamp-1; else endind=obj.nsamp; end
            numint = endind-startind+1;
            
            du_ipH_dphi = [0,0,-1/obj.dxGNAT,1/obj.dxGNAT];
            du_imH_dphi = [0,-1/obj.dxGNAT,1/obj.dxGNAT,0];
            du_im3H_dphi = [-1/obj.dxGNAT,1/obj.dxGNAT,0,0];
            
            Drhotilde_ipH_dphi = Drhotilde_ipH_Du_ipH(startind:endind)*du_ipH_dphi + Drhotilde_ipH_Du_imH(startind:endind)*du_imH_dphi;
            Drhotilde_imH_dphi = Drhotilde_imH_Du_imH(startind:endind)*du_imH_dphi + Drhotilde_imH_Du_im3H(startind:endind)*du_im3H_dphi;
            
            dRdw = [dRdw; reshape((bsxfun(@times,(bsxfun(@times,Drhotilde_ipH_dphi,u_ipH(startind:endind)) + rhotilde_ipH(startind:endind)*du_ipH_dphi),obj.A(startind:endind,2)) - ...
                bsxfun(@times,(bsxfun(@times,Drhotilde_imH_dphi,u_imH(startind:endind)) + rhotilde_imH(startind:endind)*du_imH_dphi),obj.A(startind:endind,1)))',numint*4,1)];
            
            if obj.indN
                Drhotilde_ipH_dphi = Drhotilde_ipH_Du_ipH(end)*[0;0;-1/obj.dxGNAT] + Drhotilde_ipH_Du_imH(end)*[0;-1/obj.dxGNAT;1/obj.dxGNAT];
                Drhotilde_imH_dphi = Drhotilde_imH_Du_imH(end)*[0;-1/obj.dxGNAT;1/obj.dxGNAT] + Drhotilde_imH_Du_im3H(end)*[-1/obj.dxGNAT;1/obj.dxGNAT;0];
                
                dRdw = [dRdw; obj.A(end,2)*(u_ipH(end)*Drhotilde_ipH_dphi + rhotilde_ipH(end)*[0;0;-1/obj.dxGNAT]) - ...
                        obj.A(end,1)*(u_imH(end)*Drhotilde_imH_dphi + rhotilde_imH(end)*[0;-1/obj.dxGNAT;1/obj.dxGNAT])];
            end
            
            dRdp = bsxfun(@times,obj.dA_ipH_dp,rhotilde_ipH.*u_ipH) - bsxfun(@times,obj.dA_imH_dp,rhotilde_imH.*u_imH);
            
            d2Rdw2R = [];
            d2RdwdpR = [];
        end

        function  [f,dfdw,dfdp] = ParameterEstimation(obj,w,p)
            f = 0.5*sum((w-obj.targetSoln).^2);
            dfdw = (w - obj.targetSoln);
            dfdp = zeros(size(p,1),1);
        end

        function  [f,dfdw,dfdp] = ParameterEstimationROM(obj,w,p)
            f = 0.5*sum((w-obj.targetSoln).^2);
            dfdw = (w - obj.targetSoln);
            dfdp = zeros(size(p,1),1);
        end
        
        function  [f,dfdw,dfdp] = ParameterEstimationGNAT(obj,w,p)
            f = 0.5*sum((w(obj.indOfUniqueJrow)-obj.targetSoln).^2);
            dfdw = (w(obj.indOfUniqueJrow) - obj.targetSoln);
            dfdp = zeros(size(p,1),1);
        end

        function  [f,dfdw,dfdp] = OptSample(obj,w,p)
            
            [R,dRdw,dRdp] = obj.ResSens(w,[]);
            %[fo,dfodw,dfodp] = obj.ParameterEstimation(w,p);
            
            f = -0.5*(R'*R);% + obj.Gamma*fo;
            dfdw = -dRdw'*R;% + obj.Gamma*dfodw;
            dfdp = -dRdp'*R;% + obj.Gamma*dfodp;
            
%             for i = 1:size(obj.pbar,2)
%                 f = f - obj.Gamma*(p-obj.pbar(:,i))'*(p-obj.pbar(:,i));
%                 dfdp = dfdp - obj.Gamma*(p-obj.pbar(:,i));
%             end
        end

        function  [f,dfdw,dfdp] = OptSampleROM(obj,w,p)
            
            [R,dRdw,dRdp] = obj.ResSens(w,[]);
            [fo,dfodw,dfodp] = obj.ParameterEstimationROM(w,p);

            f = -0.5*(R'*R) + obj.Gamma*fo;
            dfdw = -dRdw'*R + obj.Gamma*dfodw;
            dfdp = -dRdp'*R + obj.Gamma*dfodp;
            
%             for i = 1:size(obj.pbar,2)
%                 f = f - obj.Gamma*(p-obj.pbar(:,i))'*(p-obj.pbar(:,i));
%                 dfdp = dfdp - obj.Gamma*(p-obj.pbar(:,i));
%             end
        end

        function  [f,dfdw,dfdp] = OptSampleGNAT(obj,w,p)
            
            [RHat,dRdwHat,dRdpHat] = obj.ResSensGNAT(w,[]);
            obj.Jhat(obj.reconstJhatInd) = dRdwHat;
            [fo,dfodw,dfodp] = obj.ParameterEstimationROM(w,p);
            
            f = -0.5*(RHat'*RHat) + obj.Gamma*fo;
            dfdw = -obj.Jhat'*RHat + obj.Gamma*dfodw;
            dfdp = -dRdpHat'*RHat + obj.Gamma*dfodp;
            
%             for i = 1:size(obj.pbar,2)
%                 f = f - obj.Gamma*(p-obj.pbar(:,i))'*(p-obj.pbar(:,i));
%                 dfdp = dfdp - obj.Gamma*(p-obj.pbar(:,i));
%             end
%             
        end
        
        function  [c,dcdw,dcdp] = ConstraintSet1(obj,w,p)
%             c = -1*ones(125,1);
%             dcdw = zeros(125,126);
%             dcdp = zeros(125,1);
%             return;
            %No large mach numbers
            [u,dudphi] = obj.velocity([obj.U0;w;obj.UL]);
            [rho,drhodu] = obj.density(u);
            [M,dMdu,dMdrho] = mach(obj,rho,u);
            
            Mmax = 1.2;
            %Don't mach at last cell center because it is implicitly
            %imposed through the BC on the right and the constraint on the
            %left
%             c = M(65) - Mmax;
%             
%             dcdw = dMdu(65,:)*dudphi(:,2:end-1) + dMdrho(65,:)*drhodu*dudphi(:,2:end-1);
%             dcdp = zeros(size(c,1),size(p,1));
            
% 
%             c = - M(50:80) + Mmax;
%             
%             dcdw = -dMdu(50:80,:)*dudphi(:,2:end-1) - dMdrho(50:80,:)*drhodu*dudphi(:,2:end-1);
%             dcdp = zeros(size(c,1),size(p,1));



            c = M(2:end-1) + Mmax;
            
            dcdw = -dMdu(2:end-1,:)*dudphi(:,2:end-1) - dMdrho(2:end-1,:)*drhodu*dudphi(:,2:end-1);
            dcdp = zeros(size(c,1),size(p,1));
        end
        
        function  [c,dcdw,dcdp] = ConstraintSet2(obj,w,p)
            %No strong shocks
            [u,dudphi] = obj.velocity([obj.U0;w;obj.UL]);
            [rho,drhodu] = obj.density(u);
            [M,dMdu,dMdrho] = mach(obj,rho,u);
            
            eps2 = 0.56^2;
            c = diff(M) - eps2;
            
            dcdw = (dMdu(2:end,:)*dudphi(:,2:end-1) + dMdrho(2:end,:)*drhodu*dudphi(:,2:end-1)) ...
                - (dMdu(1:end-1,:)*dudphi(:,2:end-1) + dMdrho(1:end-1,:)*drhodu*dudphi(:,2:end-1));
            dcdp = zeros(size(c,1),size(p,1));
        end
        
        function  [] = setTargetSolution(obj,wstar)
            obj.targetSoln = wstar;
        end
        
        function  [] = setPBar(obj,PBAR)
            obj.pbar = [obj.pbar,PBAR];
        end
        
        function  [] = setGamma(obj,GAMMA)
            obj.Gamma = GAMMA;
        end
        
        function  [] = updateParameters(obj,p)
            %Update the parameters
            if obj.Atype == 1
                %Nozzle Spline, Fixed Ends
                [obj.A, obj.dAdp] = obj.NozzleAreaSpline(0.5*(obj.mesh.node(2:end,1) + obj.mesh.node(1:end-1,1)),[obj.A_Bnd(1);p(:);obj.A_Bnd(2)]);
            elseif obj.Atype == 2
                %Nozzle Spline, Adjustable Ends
                [obj.A, obj.dAdp] = obj.NozzleAreaSpline(0.5*(obj.mesh.node(2:end,1) + obj.mesh.node(1:end-1,1)),p(:));
            elseif obj.Atype == 3
                [obj.A, obj.dAdp] = obj.NozzleArea1D(0.5*(obj.mesh.node(2:end,1) + obj.mesh.node(1:end-1,1)),p);
            elseif obj.Atype == 4
                [obj.A, obj.dAdp] = obj.NozzleAreaAIAA(0.5*(obj.mesh.node(2:end,1) + obj.mesh.node(1:end-1,1)),p);
            end
        end
        
        function  [] = updateParametersGNAT(obj,p)
            %Update the parameters
            if obj.Atype == 1
                %Nozzle Spline, Fixed Ends
                [obj.A, obj.dAdp] = obj.NozzleAreaSpline(0.5*(obj.mesh.node(2:end,1) + obj.mesh.node(1:end-1,1)),[obj.A_Bnd(1);p(:);obj.A_Bnd(2)]);
                
            elseif obj.Atype == 2
                %Nozzle Spline, Adjustable Ends
                [obj.A, obj.dAdp] = obj.NozzleAreaSpline(0.5*(obj.mesh.node(2:end,1) + obj.mesh.node(1:end-1,1)),p(:));
            else
                [obj.A, obj.dAdp] = obj.NozzleArea1D(0.5*(obj.mesh.node(2:end,1) + obj.mesh.node(1:end-1,1)),p);
            end
            tmp=obj.A; obj.A=[];
            obj.A(:,1) = tmp(obj.sampleInd,1);
            obj.A(:,2) = tmp(obj.sampleInd+1,1);
            obj.dAdp = obj.dAdp(obj.sampleInd,:);
        end
        
        function  [A,dAdp] = NozzleAreaAIAA(obj,xi,p)
            Ades=p(1);
            Ainlet=p(2);
            Aexit=p(3);
            Athroat=p(4);
            k=p(5);
            alpha=p(6);
            xstar=p(7);
            xb=p(8);
            xf=p(9);
            d=p(10);
            
            n=length(xi);
            
            A=zeros(n,1);
            dAdp = zeros(n,10);
            
            for i = 1:n
                if xi(i) <= xb
                    A1 = Ainlet;
                    A2 = Ades+k*(abs(xi(i)-xstar))^alpha;
                    if A1 <= A2
                        % First flat part
                        A(i) = A1;
                        dAdp(i,2) = 1;
                    else
                        % First parabolic part
                        A(i) = A2;
                        dAdp(i,1) = 1;
                        dAdp(i,5) = (abs(xi(i)-xstar))^alpha;
                        dAdp(i,6) = real(log(xi(i)-xstar))*k*(xi(i)-xstar)^alpha;
                        dAdp(i,7) = -alpha*k*(xi(i)-xstar)^(alpha-1);
                    end
                elseif xi(i) >= xf
                    A1 = Aexit;
                    A2 = Ades+k*(abs(xi(i)-xstar))^alpha;
                    if A1 <= A2
                        % Last flat part
                        A(i) = A1;
                        dAdp(i,3) = 1;
                    else
                        % Last parabolic part
                        A(i) = A2;
                        dAdp(i,1) = 1;
                        dAdp(i,5) = (abs(xi(i)-xstar))^alpha;
                        dAdp(i,6) = log(abs(xi(i)-xstar))*k*(xi(i)-xstar)^alpha;
                        dAdp(i,7) = -alpha*k*(xi(i)-xstar)^(alpha-1);
                    end
                else
                    
                    rthroat = sqrt(Athroat/pi);
                    rcirc = rthroat + .5*d;
                    
                    rxb = sqrt((Ades+k*(abs(xb-xstar))^alpha)/pi);
                    theta1b = atan((rxb-rcirc)/(xstar-xb));
                    theta2b = asin(.5*d/sqrt((rxb-rcirc)^2+(xstar-xb)^2));
                    mbc = -tan(theta1b+theta2b);
                    theta3b = atan(-1/mbc);
                    xc = xstar - .5*d*cos(theta3b);
                    
                    rxf = sqrt((Ades+k*(abs(xf-xstar))^alpha)/pi);
                    theta1f = atan((rxf-rcirc)/(xf-xstar));
                    theta2f = asin(.5*d/sqrt((rxf-rcirc)^2+(xf-xstar)^2));
                    mdf = tan(theta1f+theta2f);
                    theta3f = atan(1/mdf);
                    xd = xstar + .5*d*cos(theta3f);
                    if xi(i) < xc
                        % Descending line
                        r = mbc*(xi(i)-xb)+rxb;
                        hyp = sqrt((rxb-rcirc)^2+(xstar-xb)^2);
                        dAdr = 2*pi*r;
                        
                        drxbdAdes = .5*((Ades+k*(abs(xb-xstar))^alpha)/pi)^-.5/pi;
                        dhypdrxb = (rxb-rcirc)/hyp;
                        dTheta2bdhyp = (1-(.5*d/hyp)^2)^(-.5)*(-.5*d/hyp^2);
                        dTheta2bdrxb = dTheta2bdhyp*dhypdrxb;
                        dTheta2bdAdes = dTheta2bdrxb*drxbdAdes;
                        dTheta1bdrxb = ((1+((rxb-rcirc)/(xstar-xb))^2)*(xstar-xb))^-1;
                        dTheta1bdAdes = dTheta1bdrxb*drxbdAdes;
                        dmbcdTheta = -sec(theta1b+theta2b)^2;
                        dmbcdAdes = dmbcdTheta*(dTheta1bdAdes+dTheta2bdAdes);
                        drdmbc = xi(i)-xb;
                        drdAdes = drdmbc*dmbcdAdes+drxbdAdes;
                        dAdp(i,1) = dAdr*drdAdes;
                        
                        drcircdAthroat = 1/sqrt(4*pi*Athroat);
                        dTheta1bdrcirc = -dTheta1bdrxb;
                        dTheta1bdAthroat = dTheta1bdrcirc*drcircdAthroat;
                        dhypdrcirc = -dhypdrxb;
                        dTheta2bdAthroat = dTheta2bdhyp*dhypdrcirc*drcircdAthroat;
                        dmbcdAthroat = dmbcdTheta*(dTheta1bdAthroat + dTheta2bdAthroat);
                        drdAthroat = drdmbc*dmbcdAthroat;
                        dAdp(i,4) = dAdr*drdAthroat;
                        
                        drxbdk = drxbdAdes*(abs(xb-xstar))^alpha;
                        drdrxb = drdmbc*dmbcdTheta*(dTheta1bdrxb+dTheta2bdrxb)+1;
                        dAdp(i,5) = dAdr*drdrxb*drxbdk;
                        
                        drxbdalpha = drxbdAdes*(k*log(abs(xb-xstar))*(abs(xb-xstar))^alpha);
                        dAdp(i,6) = dAdr*drdrxb*drxbdalpha;
                        
                        drxbdxstar = .5*((Ades+k*(abs(xb-xstar))^alpha)/pi)^(-.5)*(alpha*k/pi*(abs(xb-xstar))^(alpha-1));
                        dTheta1bdxstar = (1+((rxb-rcirc)/(xstar-xb))^2)^-1*(drxbdxstar*(xstar-xb)^-1-rxb*(xstar-xb)^-2+rcirc*(xstar-xb)^-2);
                        dhypdxstar = ((rxb-rcirc)^2+(xstar-xb)^2)^-.5*((rxb-rcirc)*drxbdxstar+(xstar-xb));
                        dTheta2bdxstar = dTheta2bdhyp*dhypdxstar;
                        dmbcdxstar = dmbcdTheta*(dTheta1bdxstar+dTheta2bdxstar);
                        drdxstar = drdmbc*dmbcdxstar+drxbdxstar;
                        dAdp(i,7) = dAdr*drdxstar;
                        
                        drxbdxb = -drxbdxstar;
                        dhypdxb = ((rxb-rcirc)^2+(xstar-xb)^2)^-.5*((rxb-rcirc)*drxbdxb-(xstar-xb));
                        dTheta2bdxb = dTheta2bdhyp*dhypdxb;
                        dTheta1bdxb = (1+((rxb-rcirc)/(xstar-xb))^2)^-1*(drxbdxb*(xstar-xb)^-1+rxb*(xstar-xb)^-2-rcirc*(xstar-xb)^-2);
                        dmbcdxb = dmbcdTheta*(dTheta1bdxb+dTheta2bdxb);
                        drdxb = dmbcdxb*(xi(i)-xb)-mbc+drxbdxb;
                        dAdp(i,8) = dAdr*drdxb;
                        
                        dhypdd = ((rxb-rcirc)^2+(xstar-xb)^2)^-.5*-(rxb-rcirc)*.5;
                        dTheta2bdd = (1-(.5*d/hyp)^2)^(-.5)*(.5*hyp^-1-.5*d*hyp^-2*dhypdd);
                        dTheta1bdd = (1+((rxb-rcirc)/(xstar-xb))^2)^-1*(-(xstar-xb)^-1*.5);
                        dmbcdd = dmbcdTheta*(dTheta1bdd+dTheta2bdd);
                        drdd = drdmbc*dmbcdd;
                        dAdp(i,10) = dAdr*drdd;
                        
                        
                    elseif xi(i) > xd
                        % Ascending line
                        r = mdf*(xi(i)-xf)+rxf;
                        hyp = sqrt((rxf-rcirc)^2+(xf-xstar)^2);
                        dAdr = 2*pi*r;
                        
                        drxfdAdes = .5*((Ades+k*(abs(xf-xstar))^alpha)/pi)^-.5/pi;
                        dhypdrxf = (rxf-rcirc)/hyp;
                        dTheta2fdhyp = (1-(.5*d/hyp)^2)^(-.5)*(-.5*d/hyp^2);
                        dTheta2fdrxf = dTheta2fdhyp*dhypdrxf;
                        dTheta2fdAdes = dTheta2fdrxf*drxfdAdes;
                        dTheta1fdrxf = ((1+((rxf-rcirc)/(xstar-xf))^2)*(xf-xstar))^-1;
                        dTheta1fdAdes = dTheta1fdrxf*drxfdAdes;
                        dmdfdTheta = sec(theta1f+theta2f)^2;
                        dmdfdAdes = dmdfdTheta*(dTheta1fdAdes+dTheta2fdAdes);
                        drdmdf = xi(i)-xf;
                        drdAdes = drdmdf*dmdfdAdes+drxfdAdes;
                        dAdp(i,1) = dAdr*drdAdes;
                        
                        drcircdAthroat = 1/sqrt(4*pi*Athroat);
                        dTheta1fdrcirc = -dTheta1fdrxf;
                        dTheta1fdAthroat = dTheta1fdrcirc*drcircdAthroat;
                        dhypdrcirc = -dhypdrxf;
                        dTheta2fdAthroat = dTheta2fdhyp*dhypdrcirc*drcircdAthroat;
                        dmdfdAthroat = dmdfdTheta*(dTheta1fdAthroat + dTheta2fdAthroat);
                        drdAthroat = drdmdf*dmdfdAthroat;
                        dAdp(i,4) = dAdr*drdAthroat;
                        
                        drxfdk = drxfdAdes*(abs(xf-xstar))^alpha;
                        drdrxf = drdmdf*dmdfdTheta*(dTheta1fdrxf+dTheta2fdrxf)+1;
                        dAdp(i,5) = dAdr*drdrxf*drxfdk;
                        
                        drxfdalpha = drxfdAdes*(k*log(abs(xf-xstar))*(abs(xf-xstar))^alpha);
                        dAdp(i,6) = dAdr*drdrxf*drxfdalpha;
                        
                        drxfdxstar = -.5*((Ades+k*(abs(xf-xstar))^alpha)/pi)^(-.5)*(alpha*k/pi*(abs(xf-xstar))^(alpha-1));
                        dTheta1fdxstar = (1+((rxf-rcirc)/(xf-xstar))^2)^-1*(drxfdxstar*(xf-xstar)^-1+rxf*(xf-xstar)^-2-rcirc*(xf-xstar)^-2);
                        dhypdxstar = ((rxf-rcirc)^2+(xstar-xf)^2)^-.5*((rxf-rcirc)*drxfdxstar+(xstar-xf));
                        dTheta2fdxstar = dTheta2fdhyp*dhypdxstar;
                        dmdfdxstar = dmdfdTheta*(dTheta1fdxstar+dTheta2fdxstar);
                        drdxstar = drdmdf*dmdfdxstar+drxfdxstar;
                        dAdp(i,7) = dAdr*drdxstar;
                        
                        drxfdxf = -drxfdxstar;
                        dhypdxf = ((rxf-rcirc)^2+(xstar-xf)^2)^-.5*((rxf-rcirc)*drxfdxf-(xstar-xf));
                        dTheta2fdxf = dTheta2fdhyp*dhypdxf;
                        dTheta1fdxf = -(1+((rxf-rcirc)/(xstar-xf))^2)^-1*(drxfdxf*(xstar-xf)^-1+rxf*(xstar-xf)^-2-rcirc*(xstar-xf)^-2);
                        dmdfdxf = dmdfdTheta*(dTheta1fdxf+dTheta2fdxf);
                        drdxf = dmdfdxf*(xi(i)-xf)-mdf+drxfdxf;
                        dAdp(i,9) = dAdr*drdxf;
                        
                        dhypdd = ((rxf-rcirc)^2+(xstar-xf)^2)^-.5*-(rxf-rcirc)*.5;
                        dTheta2fdd = (1-(.5*d/hyp)^2)^(-.5)*(.5*hyp^-1-.5*d*hyp^-2*dhypdd);
                        dTheta1fdd = -(1+((rxf-rcirc)/(xstar-xf))^2)^-1*(-(xstar-xf)^-1*.5);
                        dmdfdd = dmdfdTheta*(dTheta1fdd+dTheta2fdd);
                        drdd = drdmdf*dmdfdd;
                        dAdp(i,10) = dAdr*drdd;
                        
                    else
                        % Circular part
                        r = rcirc-sqrt((.5*d)^2-(xi(i)-xstar)^2);
                        dAdr = 2*pi*r;
                        dAdp(i,4) = dAdr/sqrt(4*pi*Athroat);
                        dAdp(i,7) = dAdr*-((d/2)^2-(xi(i)-xstar)^2)^(-1/2)*(xi(i)-xstar);
                        dAdp(i,10) = dAdr*(.5-.25*d*((d/2)^2-(xi(i)-xstar)^2)^-.5);
                    end
                    A(i) = pi*r^2;
                end
            end
            
        end
        
        function  [A,dAdp] = NozzleAreaSpline(obj,xi,p)
            %obj.SpPts = [0;1;3;4];
            
            Dx = diff(obj.SpPts);
            
            M = diag([0;Dx(2:end)],1) + diag([Dx(1:end-1); 0],-1) + diag([0;2*(Dx(1:end-1)+Dx(2:end));0],0);
            M(1,1) = Dx(2);
            M(end,end) = Dx(end-1);
            M(1,2) = -(Dx(1) + Dx(2));
            M(end,end-1) = -(Dx(end) + Dx(end-1));
            M(1,3) = Dx(1);
            M(end,end-2) = Dx(end);
            
            f_sig = zeros(obj.nSpPts,1);
            f_sig(2:end-1,1)  = 3*((p(3:end) - p(2:end-1))./Dx(2:end) - (p(2:end-1) - p(1:end-2))./Dx(1:end-1));
            
            f_dsig = zeros(obj.nSpPts, obj.nSpPts);
            f_dsig(2:end-1,1:end-2) = diag(3./Dx(1:end-1));
            f_dsig(2:end-1,2:end-1) = f_dsig(2:end-1,2:end-1) + diag(-3*(1./Dx(1:end-1) + 1./Dx(2:end)));
            f_dsig(2:end-1,3:end) = f_dsig(2:end-1,3:end) + diag(3./Dx(2:end));
            
            sig = M\f_sig;
            dsig = M\f_dsig;
            
            a = (sig(2:end) - sig(1:end-1))./(3*Dx);
            b = sig(1:end-1);
            c = -(sig(2:end) + 2*sig(1:end-1)).*(Dx/3) + (p(2:end) - p(1:end-1))./Dx;
            d = p(1:end-1);
            
            dadp = repmat(1./(3*Dx),1,obj.nSpPts).*(dsig(2:end,:) - dsig(1:end-1,:));
            dbdp = dsig(1:end-1,:);
            
            dcdp = zeros(obj.nSpPts-1,obj.nSpPts);
            dcdp(:,1:end-1) = -diag(1./Dx);
            dcdp(:,2:end) = dcdp(:,2:end) + diag(1./Dx);
            dcdp = dcdp + repmat(-Dx/3,1,obj.nSpPts).*(dsig(2:end,:) + 2*dsig(1:end-1,:));
            
            dddp = [eye(obj.nSpPts-1),zeros(obj.nSpPts-1,1)];
            
            n = length(xi);
            A = zeros(n,1);
            dAdp = zeros(n,obj.nSpPts);
            for i = 1:n
                tmp = find(obj.SpPts < xi(i));
                ind = tmp(end);
                
                %A(i) = a(ind)*(xi(i)-obj.SpPts(ind))^3 + b(ind)*(xi(i)-obj.SpPts(ind))^2 + c(ind)*(xi(i)-obj.SpPts(ind)) + d(ind);
                %if i == 1 || i == n
                %  continue;
                %end
                %Want the sensitivity of A with respect to the inlet and
                %outlet area to be zero because I want these fixed
                A(i) = a(ind)*(xi(i)-obj.SpPts(ind))^3 + b(ind)*(xi(i)-obj.SpPts(ind))^2 + c(ind)*(xi(i)-obj.SpPts(ind)) + d(ind);
                dAdp(i,:) = dadp(ind,:)*(xi(i)-obj.SpPts(ind))^3 + dbdp(ind,:)*(xi(i)-obj.SpPts(ind))^2 + dcdp(ind,:)*(xi(i)-obj.SpPts(ind)) + dddp(ind,:);
            end
            
            if obj.Atype == 1
                dAdp = dAdp(:,2:end-1);
            end
        end

        function  [A,dAdp] = NozzleArea1D(obj,xi,p)
            
            A = 0.6*(xi-1).^2 + p;
            dAdp = ones(size(A));
            
        end
    end
end
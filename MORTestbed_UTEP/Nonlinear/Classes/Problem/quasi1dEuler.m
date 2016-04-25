classdef quasi1dEuler < handle
    
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
        
        sigma0=10;
        beta;
        nondim;
        
        nVol;
        A_Bnd;
        Atype = 1;
                
        S;
        SVol;
        dSdp;
        p;
        gamma;
        M0;
        T0;
        ShockLocRatio;
%         P0;
        pExitIncrease;
        pExit;
        rhoExit;
        uExit;
        R;
        Cv;
        c0;
        Pt;
        Tt;
        astar2;
        
        xVol;
        
        nSpPts;
        SpPts;
        
        Gamma;
        targetSoln = []; %This is going to be the target solution for parameter estimation
        
        pbar; 
        
        ind1;
        indN;
        nVolSamp;
        nVolMask;
        
        iarray;
        iarrayMask;
        iarrayFace;
        iarrayFaceI;
        iarraySFull;
    end
    
%     properties (Hidden=true)
%         pbar; 
%         
%         ind1;
%         indN;
%         nVolSamp;
%         nVolMask;
%         
%         iarray;
%         iarrayMask;
%         iarrayFace;
%         iarrayFaceI;
%         iarraySFull;
%     end
    
    methods
        %General Functions
        function  [obj] = quasi1dEuler(cfgobj,oldobj)
            %This is the constructor for the 1D Euler Nozzle Eqn.  If
            %there is 1 input, the constructor will create the class
            %instance from the variables defined in this function.  If
            %there are 2 inputs, the first input will be unused and the
            %output will be a deep copy of the oldobj.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %cfgobj - CONFIG object
            %oldobj - quasiEuler1D object whose properties will be copied
            %         into the new object.
            %
            %Outputs:
            %--------
            %obj - quasiEuler1D object
            %--------------------------------------------------------------
            
            if nargin == 0
                %Keep everything at defaults if no input arguments
                return;
            end
            
            % Construct a new object based on a deep copy of an old object
            % of this class by copying properties over.
            if nargin == 2 && isa(oldobj,'quasi1dEuler')
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
            obj.mesh.gridpt = linspace(cfgobj.DLim(1,1),cfgobj.DLim(1,2),cfgobj.nNodes(1))';
            obj.mesh.node = 0.5*(obj.mesh.gridpt(1:end-1) + obj.mesh.gridpt(2:end)); %Location of cell-centers
            obj.mesh.conn = [1:cfgobj.nNodes(1)-1;2:cfgobj.nNodes(1)]';
            obj.mesh.dbc = [];
                        
            %Number of unknowns per node
            obj.ndim = 3;
            
            %Number of control volumes
            obj.nVol = cfgobj.nNodes(1)-1;
            
            %Number of dofs is the number of grid points
            obj.ndof = obj.ndim*obj.nVol;
                    
            %Extract parameters
            obj.gamma = obj.config.param{1};
            obj.R = obj.config.param{2};
            obj.Pt = obj.config.param{3};
            obj.Tt = obj.config.param{4};
            obj.pExitIncrease = obj.config.param{5};
            obj.p = obj.config.param{6}; obj.p = obj.p(:);
            AreaParam = obj.config.param{7};
            obj.M0 = obj.config.param{8};
            obj.ShockLocRatio = obj.config.param{9};
            if length(obj.config.param)==8
                obj.A_Bnd=obj.config.param{8};
            end
            
            % Compute S(x) and SVol(x)
            switch AreaParam
                case 'splineFixedEnds'
                    %Determine number of spline knots to parametrize shape
                    obj.nSpPts = length(obj.p)+2;
                    obj.p = [obj.A_Bnd(1);obj.p;obj.A_Bnd(2)];
                    obj.SpPts = linspace(cfgobj.DLim(1,1),cfgobj.DLim(1,2),obj.nSpPts)';
                    obj.A_Bnd=obj.A_Bnd/obj.nondim.L;
                    obj.Atype = 1;
                case 'splineFreeEnds'
                    %Determine number of spline knots to parametrize shape
                    obj.nSpPts = length(obj.p);
                    obj.SpPts = linspace(cfgobj.DLim(1,1),cfgobj.DLim(1,2),obj.nSpPts)';
                    obj.Atype = 2;
                case '1D'
                    obj.Atype = 3;
                case 'AIAA'
                    obj.Atype = 4;
            end
            [obj.S,obj.SVol,obj.xVol,obj.dSdp] = obj.computeNozzleArea();
            plot(obj.SVol)
            factor = 1;
            obj.S = factor*obj.S;
            obj.SVol = factor*obj.SVol;
            obj.xVol = factor*obj.xVol;
            obj.dSdp = factor*obj.dSdp;
            Ainlet = obj.SVol(1);
            Aexit  = obj.SVol(end);
            Athroat = min(min(obj.S),min(obj.SVol));
            
            %Extract information from altfile
            ShockLoc = ceil(obj.ShockLocRatio*obj.ndof/3);
            [~, P, rho, u, T] = feval(cfgobj.altFile,obj.SVol,Athroat,obj.M0,...
                                      obj.gamma,obj.R,obj.Tt,obj.Pt,obj.nVol,ShockLoc,obj.pExitIncrease);
            obj.pExit = P(end);
            obj.rhoExit = rho(end);
            obj.uExit = u(end);
            
            %Nondimensionalization values
            rho0 = rho(1);
            u0   = u(1);

            %Parameters to nondimensionalize with
            obj.nondim.u   = u0;
            obj.nondim.rho = rho0;
            obj.nondim.p   = rho0*u0*u0;
            obj.nondim.L   = cfgobj.DLim(1,2) - cfgobj.DLim(1,1);
            obj.nondim.t = obj.nondim.L/obj.nondim.u;
            
            obj.ic = [rho(1)/obj.nondim.rho;u(1)/obj.nondim.u;P(1)/obj.nondim.p;...
                      reshape(obj.primitiveToConservative(rho(2:end-1)/obj.nondim.rho,u(2:end-1)/obj.nondim.u,P(2:end-1)/obj.nondim.p)',obj.ndim*(obj.nVol-2),1);...
                      rho(end)/obj.nondim.rho;u(end)/obj.nondim.u;P(end)/obj.nondim.p];            
                  
            obj.Cv = (obj.R/(obj.gamma-1))/(obj.nondim.u^2);
            obj.astar2 = 2*obj.gamma*(obj.gamma-1)/(obj.gamma+1)*obj.Cv*obj.Tt;
            obj.Pt = obj.Pt/obj.nondim.p;      
            
            obj.dx = diff(obj.mesh.gridpt)';
            obj.dx=(1/obj.nondim.L)*obj.dx;
            
            obj.S=(1/(obj.nondim.L^2))*obj.S;
            obj.SVol=(1/obj.nondim.L^2)*obj.SVol;
            obj.xVol=(1/obj.nondim.L)*obj.xVol;
            obj.dSdp=(1/obj.nondim.L)*obj.dSdp;
 
            obj.beta = obj.gamma-1;
            
            %Setup Jacobian Structure
            obj.JacobianStructure();
        end
                
        function U = primitiveToConservative(obj,rho,u,P)
            
            U(:,1) = rho;
            U(:,2) = rho.*u;
            U(:,3) = P/(obj.gamma-1) + rho.*u.^2/2;
            
        end
        
        function  [rho,u,P,c,dcdU] = conservativeToPrimitive(obj,U)
            %This function extracts the primitive variables from the 3 x N
            %array of conservative variables.
            %--------------------------------------------------------------
            %Input:
            %------
            %U - array of conservative variables (first row is rho, second
            %    row is rho*u, third row is e).
            %
            %Output:
            %-------
            %rho  - 1 x N array of density
            %u    - 1 x N array of velocity
            %P    - 1 x N array of pressure
            %c    - 1 x N array of sound speed
            %dcdU - 3 x N array of sound speed derivatives; dcdU(:,i) =
            %       [dc_drho(i); dc_drhou(i); dc_de(i)]
            %--------------------------------------------------------------
            
            rho = U(1,:);
            u = U(2,:)./U(1,:);
            P = (obj.gamma-1)*(U(3,:) - 0.5*rho.*u.*u);
            c=sqrt(obj.gamma*P./rho);
            dcdU = [0.5*obj.gamma./(c.*rho).*(0.5*(obj.gamma-1)*u.*u - P./rho);...
                    -0.5*obj.gamma*(obj.gamma-1)*u./(rho.*c);...
                    0.5*obj.gamma*(obj.gamma-1)./(rho.*c)];
            %M = u./ (sqrt(obj.gamma*P./rho));
        end
        
        function  [R,J,Cp,Cpp,dCp,dCpp] = ResJac(obj,U,~)
            %This function computes the residual and jacobian of the
            %quasi-1D Euler equations.  Also, terms needed for the time
            %stepping residual and jacobian are computed.  Governing
            %equations are:
            %
            %Inlet: 
            %    [ P - Pt*(1-(gamma-1)/(gamma+1)*u^2/astar2)^(gamma/(gamma-1));...
            %      rho - P/((gamma-1)*Cv*T(u));...
            %      dPdt - rho*c*dudt+(u-c)*(dPdx - rho*c*dudx) - gamma*P*(u/S)*dSdx] = 0
            %which can be written:
            % Cpp*dVdt + [ P - Pt*(1-(gamma-1)/(gamma+1)*u^2/astar2)^(gamma/(gamma-1));...
            %              rho - P/((gamma-1)*Cv*T(u));...
            %              (u-c)*(dPdx - rho*c*dudx) -  gamma*P*(u/S)*dSdx] = 0
            %where
            %  Cpp = [0,0,0;0,0,0;0,-rho*c,1]
            %
            %Governing:
            %    dUdt + (1/S)*d(FS)dx = Q
            %
            %Outlet:
            %    Cp*dVdt = -Lambda*Cp*dVdx + Q'  
            %
            % where V = vector of primitive variables, U = vector of
            % conservative variables.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj - quasiEuler1D object
            %U   - state vector: U(1:3) = primitive variables at inlet,
            %      U(end-2:end) = primitive variables at outlet, and
            %      U(4:end-3) = conservative variables on domain interior.
            %
            %Outputs:
            %--------
            %R    - nonlinear residual (without time derivative terms)
            %J    - nonlinear jacobian
            %Cp   - matrix multiplying outlet time derivative term
            %Cpp  - matrix multiplying inlet time derivative term
            %dCp  - jacobian of Cp
            %dCpp - jacobian of Cpp
            %--------------------------------------------------------------
           
            %Reshape U from a stacked vector ([rho_1;rho_1*u_1;e_1;...]) of
            %size 3nVol x 1 to a matrix of size 3 x nVol where U(1,:) =
            %rho, U(2,:) = rho.*u, and U(3,:) = e.  Except first three and
            %last three values are primitive variables while others are
            %conservative.
            U = reshape(U,3,obj.nVol);
            %Extract primitive variables.  Take into account that the first
            %three and last three variables are PRIMITIVE VARIABLES while
            %all other are conservative variables.
            [rho,u,P,c] = obj.conservativeToPrimitive(U(:,2:end-1));
            rho = [U(1,1),rho,U(1,end)]; %Density
            u   = [U(2,1),u,U(2,end)]; %Velocity
            P   = [U(3,1),P,U(3,end)]; %Pressure
            c   = [sqrt(obj.gamma*P(1)/rho(1)),c,sqrt(obj.gamma*P(end)/rho(end))]; %Speed of sound
            e   = [P(1)/(obj.gamma-1)+rho(1)*u(1)^2/2,U(3,2:end-1),P(end)/(obj.gamma-1)+rho(end)*u(end)^2/2]; %Energy
            %Derivative of sound speed w.r.t. primitive variables.
            %dc_prim = [dcdrho_1, ..., dcdrho_nVol; ...
            %           dcdu_1  , ..., dcdu_nVol; ...
            %           dcdP_1  , ..., dcdP_nVol]
            dc_prim = [-0.5*obj.gamma*P./(c.*rho.*rho);zeros(1,obj.nVol);0.5*obj.gamma./(c.*rho)]';
            %Derivative of sound speed w.r.t. conservative variables.
            %dc_cons = [dcdrho_1 , ..., dcdrho_nVol; ...
            %           dcdrhou_1, ..., dcdrhou_nVol; ...
            %           dcde_1  , ..., dcde_nVol]
            dc_cons = [0.5*obj.gamma./(c.*rho).*(0.5*(obj.gamma-1)*u.*u - P./rho);...
                      -0.5*obj.gamma*(obj.gamma-1)*u./(rho.*c);...
                       0.5*obj.gamma*(obj.gamma-1)./(rho.*c)]';
            
            %Get residual and Jacobian contributions from inlet, governing
            %equations (interior), and outlet.  Also compute matrices
            %multiplying time derivatives.
            if nargout == 1
                R1=obj.fullImplicitInletBC(rho(1:2),u(1:2),P(1:2),c(1:2),dc_prim(1:2,:));
                R2=obj.governEqn(rho,u,P,c,e,dc_cons); %R2 is a 3 x nVol-2
                R3=obj.fullImplicitOutletBC(rho(end-1:end),u(end-1:end),P(end-1:end),c(end-1:end),dc_prim(end-1:end,:));

                R=[R1;R2(:);R3];
                
                J=[]; Cp=[]; dCp=[]; Cpp=[]; dCpp=[];
            else
                [R1,J1,Cpp,dCpp]=obj.fullImplicitInletBC(rho(1:2),u(1:2),P(1:2),c(1:2),dc_prim(1:2,:));
                [R2,J2]=obj.governEqn(rho,u,P,c,e,dc_cons); %R2 is a 3 x nVol-2
                [R3,J3,Cp,dCp]=obj.fullImplicitOutletBC(rho(end-1:end),u(end-1:end),P(end-1:end),c(end-1:end),dc_prim(end-1:end,:));
                
                %Form full residual and Jacobian.
                R=[R1;R2(:);R3];
                J=[J1,zeros(3,3*(obj.nVol-2)); J2; zeros(3,3*(obj.nVol-2)), J3];
            end
        end
        
        function  [roeF,droeF] = roeFlux(obj,rho,u,P,c,e,dc)
            %This function computes the Roe flux and Jacobian at each
            %inter-cell face (since there are nVol cells, there are nVol-1
            %interfaces). 
            %Roe Flux at the interface between cell i and i+1:
            %F_{i+1/2} = 0.5*(F_{i} + F_{i+1}) - 0.5*|\hat{A}_{i+1/2}|*(U_{i+1}-U_{i})
            %for i = 1, ..., nVol-1 (need for each INTERIOR face)
            %This requires F_{i}, U_{i} for i = 1, ..., nVol
            %--------------------------------------------------------------
            %Input:
            %------
            %obj - quasi1dEuler object
            %rho - 1xnVol array of densities in each volume
            %u   - 1xnVol array of velocities in each volume
            %P   - 1xnVol array of pressures in each volume
            %c   - 1xnVol array of sound speeds in each volume
            %e   - 1xnVol array of energies in each volume
            %dc  - 3xnVol array of derivatives of sound speed w.r.t.
            %      CONSERVATIVE VARIABLES in each volume.
            %      dc = [dcdrho_1 , ..., dcdrho_nVol; ...
            %            dcdrhou_1, ..., dcdrhou_nVol; ...
            %            dcde_1  , ..., dcde_nVol]
            %
            %Output:
            %-------
            %roeF  - 3x(nVol-1) array of Roe flux at each interface.  The
            %        first dimension respresents the component of the flux
            %        (rho, rho*u, e) and the second dimension represents
            %        the interface number.
            %droeF - 3x6x(nVol-1) array of Roe flux derivatives at each
            %        interface.  The first dimension represents the
            %        component of the flux (rho, rho*u, e), the second
            %        dimension represents the derivative taken w.r.t. and
            %        the third dimension represents the interface number.
            %        Since F_{i+1/2} only depends on U_{i} and U_{i+1}, the
            %        Jacobian of F_{i+1/2} only depends on the state at
            %        volumes i and i+1.  Therefore (in the second dimension
            %        of drhoF), instead of storing the derivative w.r.t.
            %        3*nVol variables (where most are zero), we store the
            %        derivative w.r.t. only 6 variables (rho_i, rhou_i,
            %        e_i, rho_{i+1}, rhou_{i+1}, e_{i+1}) 
            %--------------------------------------------------------------

            %Compute nonlinear term at each cell volume
            F=obj.eulerNonlin(rho,u,e,P);
            %Compute the Roe adjusted values at each interface as well as
            %their derivatives w.r.t. state.  See rhoTerms for description
            %of each variable.
            [rhoH,uH,cH,drhoH,duH,dcH] = obj.roeTerms(rho,u,P,e);
     
            %Initialize roeF and droeF
            roeF = 0.5*(F(:,1:end-1)+F(:,2:end));
            droeF = zeros(3,6,obj.nVol-1);
            for i = 1:obj.nVol-1
                %Form conservative state variables at volume i and i+1
                U_ip1=[rho(i+1);rho(i+1)*u(i+1);e(i+1)];
                U_i  =[rho(i);rho(i)*u(i);e(i)];
                dU = U_ip1-U_i;
                
                %Compute terms for \hat{A} and its derivative w.r.t. state
                [S,Si,CA,CAi,Lam,dS,dSi,dCA,dCAi] = obj.eulerJac(rhoH(i),uH(i),cH(i));
                
                %Apply entropy correction to Lam
                [Lam(1,1),Lam(2,2),Lam(3,3),dlam_u,dlam_upc,dlam_umc]=obj.entropyCorrect(rho(i:i+1),u(i:i+1),c(i:i+1),dc(i:i+1,:),uH(i),cH(i),duH(i,:,:),dcH(i,:,:));
                sgn=sign(diag(Lam));
                
                %Derivative contributions from absolute value
                dlam_u=sgn(1)*dlam_u;
                dlam_upc=sgn(2)*dlam_upc;
                dlam_umc=sgn(3)*dlam_umc;
                
                %Form \hat{A}
                Ahat=(Si*CAi*abs(Lam)*CA*S);
                
                %Compute the derivative of \hat{A} w.r.t. state time dU (chain rule).
                %Data will be stored according to:
                %dAhatTimesDU(:,1,1) = dAhat/drho_{i}*dU
                %dAhatTimesDU(:,1,2) = dAhat/drho_{i+1}*dU
                %The middle index corresponds to the derivative that is
                %taken: 1 = rho, 2 = rho*u, 3 = e
                dAhatTimesDU = zeros(3,3,2);
                for j = 1:2 %j = 1 -> d/d()_i, j = 2 -> d/d()_{i+1}
                    %dAhat/drho * dU
                    dAhatTimesDU(:,1,j) = (dSi(:,:,1)*drhoH(i,j)+dSi(:,:,2)*duH(i,j,1))*(CAi*(abs(Lam)*(CA*(S*dU)))) + ...
                        Si*((dCAi(:,:,1)*drhoH(i,j)+dCAi(:,:,2)*dcH(i,j,1))*(abs(Lam)*(CA*(S*dU)))) + ...
                        Si*(CAi*(diag([dlam_u(j,1),dlam_upc(j,1),dlam_umc(j,1)])*(CA*(S*dU)))) + ...
                        Si*(CAi*(abs(Lam)*((dCA(:,:,1)*drhoH(i,j)+dCA(:,:,2)*dcH(i,j,1))*(S*dU)))) + ...
                        Si*(CAi*(abs(Lam)*(CA*((dS(:,:,1)*drhoH(i,j)+dS(:,:,2)*duH(i,j,1))*dU))));
                    %dAhat/drhou * dU
                    dAhatTimesDU(:,2,j) = (dSi(:,:,2)*duH(i,j,2))*(CAi*(abs(Lam)*(CA*(S*dU)))) + ...
                        Si*((dCAi(:,:,2)*dcH(i,j,2))*(abs(Lam)*(CA*(S*dU)))) + ...
                        Si*(CAi*(diag([dlam_u(j,2),dlam_upc(j,2),dlam_umc(j,2)])*(CA*(S*dU)))) + ...
                        Si*(CAi*(abs(Lam)*((dCA(:,:,2)*dcH(i,j,2))*(S*dU)))) + ...
                        Si*(CAi*(abs(Lam)*(CA*((dS(:,:,2)*duH(i,j,2))*dU))));
                    %dAhat/de * dU
                    dAhatTimesDU(:,3,j) = Si*((dCAi(:,:,2)*dcH(i,j,3))*(abs(Lam)*(CA*(S*dU)))) + ...
                        Si*(CAi*(diag([dlam_u(j,3),dlam_upc(j,3),dlam_umc(j,3)])*(CA*(S*dU)))) + ...
                        Si*(CAi*(abs(Lam)*((dCA(:,:,2)*dcH(i,j,3))*(S*dU))));
                end
                
                %Compute Roe flux between volume i and i+1
                roeF(:,i) = roeF(:,i) - 0.5*Ahat*dU;
                
                %Compute Jacobian of Roe flux between volumne i and i+1
                [S,Si,CA,CAi,Lam] = obj.eulerJac(rho(i),u(i),c(i));
                Jblk1 = 0.5*(Si*CAi*Lam*CA*S - dAhatTimesDU(:,:,1) + Ahat);
                
                [S,Si,CA,CAi,Lam] = obj.eulerJac(rho(i+1),u(i+1),c(i+1));
                Jblk2 = 0.5*(Si*CAi*Lam*CA*S - dAhatTimesDU(:,:,2) - Ahat);
                               
                droeF(:,:,i)=[Jblk1,Jblk2];
            end
            
        end
        
        function  [lam_u,lam_upc,lam_umc,dlam_u,dlam_upc,dlam_umc] = entropyCorrect(obj,rho,u,c,dc,uH,cH,duH,dcH)
            %This function applies entropy correction to the flow
            %eigenvalues: u, u+c, and u-c.  Let lambda represent one of
            %these eigenvalues and \hat{lam} be the corresponding
            %eigenvalue evaluated at "hatted" variables (uH, cH).  Then,
            %define eps = sigma0*max(0,\hat{lam}-lam_i,lam_{i+1}-\hat{lam})
            %and if |\hat{lam}| < eps apply entropy correction as 
            %\hat{lam} <- 0.5*(\hat{lam}^2/eps+eps).
            %--------------------------------------------------------------
            %Input:
            %------
            %obj - quasi1dEuler object
            %rho - 1x2 array of densities at volumes i and i+1
            %u   - 1x2 array of velocities at volumes i and i+1
            %P   - 1x2 array of pressures at pressures i and i+1
            %c   - 1x2 array of sounds speeds at volumes i and i+1
            %dc  - 3x2 array of derivatives of sound speed w.r.t.
            %      CONSERVATIVE VARIABLES in at volumes i and i+1.
            %      dc = [dcdrho_i , dcdrho_{i+1}; ...
            %            dcdrhou_i, dcdrhou_{i+1}; ...
            %            dcde_i   , dcde_{i+1}]
            %uH  - Roe averaged velocities at interface between volumes i
            %      and i+1 
            %cH  - Roe averaged sound speed at interface between volumes i
            %      and i+1 
            %duH - 1x2x2 derivative of Roe averaged velocity at interface
            %      between volume i and i+1 w.r.t. conservative variables
            %      at volumes i and i+1.  The second dimension represents
            %      the volume that derivatives are taken w.r.t. (1 ->
            %      volume i, 2 -> volume i+1) and the third dimension
            %      indicates the conservative variable the derivative is
            %      taken with respect to.
            %dcH - same explanation as duH except with sound speeds
            %
            %Output:
            %-------
            %lam_u    - entropy corrected u eigenvalue at interface between 
            %           volume i and i+1
            %lam_upc  - entropy corrected u+c eigenvalue at interface
            %           between  volume i and i+1
            %lam_umc  - entropy corrected u-c eigenvalue at interface
            %           between volume i and i+1
            %dlam_u   - 1x6 derivative of entropy corrected u eigenvalue at
            %           interface between volume i and i+1 w.r.t.
            %           conservative variables at volume i and i+1,
            %           respectively. [dudrho_i,dudrhou_i,dude_i,dudrho_{i+1},dudrhou_{i+1},dude_{i+1}]
            %dlam_upc - same explanation as dlam_u for eigenvalue u+c
            %dlam_umc - same explanation as dlam_u for eigenvalue u-c
            %--------------------------------------------------------------

            %First eigenvalue (u)
            lamH=uH;
            lam = u;
            lam_u=lamH;
            %Determine eps and track which value was selected (for purposes
            %on computing derivatives of the max function numerically).
            [eps,ind]=max([0,lamH-lam(1),lam(2)-lamH],[],2);
            eps=obj.sigma0*eps;
            %Apply entropy correction if necessary
            entrop=(abs(lamH)<eps);
            if entrop
                lam_u=0.5*(lamH.^2./eps+eps);                
            end
            %lam_u = lamH.*entrop+0.5*(lamH.^2./eps+eps).*(abs(lamH)<eps);
            
            %deps = 2x3 because this function only works for one volume at
            %a time and the six columns are to hold the derivatives w.r.t.
            %rho_i, rhou_i, e_i, rho_{i+1}, rhou_{i+1}, e_{i+1}
            dlam_u = [duH(1,1,1),duH(1,1,2),duH(1,1,3);...
                      duH(1,2,1),duH(1,2,2),duH(1,2,3)];
%             dlam_u = [squeeze(duH(1,1,:))';squeeze(duH(1,2,:))'];
            if ind == 1
                deps=zeros(2,3);
            elseif ind == 2
                deps = obj.sigma0*(dlam_u - [-u(1)/rho(1),1/rho(1),0;zeros(1,3)]);
            elseif ind == 3
                deps = obj.sigma0*([zeros(1,3);-u(2)/rho(2),1/rho(2),0] - dlam_u);
            end
            %Apply entropy correction to derivative if necessary
            if entrop
                dlam_u = 0.5*(2*(lamH/eps)*dlam_u - (lamH/eps)^2*deps + deps);
            end
            
            %Second eigenvalue (u + c) - see u eigenvalue for additional
            %comments
            lamH=uH+cH;
            lam = u+c;
            lam_upc=lamH;
            [eps,ind]=max([0,lamH-lam(1),lam(2)-lamH],[],2);
            eps=obj.sigma0*eps;
            entrop=(abs(lamH)<eps);
            if entrop
                lam_upc = 0.5*(lamH.^2./eps+eps);
            end
            
            dlam_upc = [duH(1,1,1)+dcH(1,1,1),duH(1,1,2)+dcH(1,1,2),duH(1,1,3)+dcH(1,1,3);...
                        duH(1,2,1)+dcH(1,2,1),duH(1,2,2)+dcH(1,2,2),duH(1,2,3)+dcH(1,2,3)];
%             dlam_upc = [squeeze(duH(1,1,:)+dcH(1,1,:))';squeeze(duH(1,2,:)+dcH(1,2,:))'];
            if ind == 1
                deps=zeros(2,3);
            elseif ind == 2
                deps = obj.sigma0*(dlam_upc - [-u(1)/rho(1)+dc(1,1),1/rho(1)+dc(1,2),dc(1,3);zeros(1,3)]);
            elseif ind == 3
                deps = obj.sigma0*([zeros(1,3);-u(2)/rho(2)+dc(2,1),1/rho(2)+dc(2,2),dc(2,3)] - dlam_upc);
            end
            
            if entrop
                dlam_upc = 0.5*(2*(lamH/eps)*dlam_upc - (lamH/eps)^2*deps + deps);
            end
            
            %Third eigenvalue (u - c) - see u eigenvalue above for
            %additional comments
            lamH=uH-cH;
            lam = u-c;
            lam_umc=lamH;
            [eps,ind]=max([0,lamH-lam(1),lam(2)-lamH],[],2);
            eps=obj.sigma0*eps;
            entrop=(abs(lamH)<eps);
            if entrop
                lam_umc = 0.5*(lamH.^2./eps+eps);
            end
            
            dlam_umc = [duH(1,1,1)-dcH(1,1,1),duH(1,1,2)-dcH(1,1,2),duH(1,1,3)-dcH(1,1,3);...
                        duH(1,2,1)-dcH(1,2,1),duH(1,2,2)-dcH(1,2,2),duH(1,2,3)-dcH(1,2,3)];
%             dlam_umc = [squeeze(duH(1,1,:)-dcH(1,1,:))';squeeze(duH(1,2,:)-dcH(1,2,:))'];
            if ind == 1
                deps=zeros(2,3);
            elseif ind == 2
                deps = obj.sigma0*(dlam_umc - [-u(1)/rho(1)-dc(1,1),1/rho(1)-dc(1,2),-dc(1,3);zeros(1,3)]);
            elseif ind == 3
                deps = obj.sigma0*([zeros(1,3);-u(2)/rho(2)-dc(2,1),1/rho(2)-dc(2,2),-dc(2,3)] - dlam_umc);
            end
            
            if entrop
                dlam_umc = 0.5*(2*(lamH/eps)*dlam_umc - (lamH/eps)^2*deps + deps);
            end
        end
        
        function  [rhoH,uH,cH,drhoH,duH,dcH] = roeTerms(obj,rho,u,P,e)
            %This function returns the Roe averaged density, velocity, and
            %sound speed at the faces, and their derivatives.
            %--------------------------------------------------------------
            %Input:
            %------
            %obj          - quasi1dEuler object
            %rho, u, P, e - 1 x nVol array of density, velocity, pressure,
            %               and energy, respectively
            %
            %Output:
            %-------
            %rhoH  - 1 x nVol-1 array of Roe adjusted densities at cell
            %        interfaces.  rhoH(i) = \hat{\rho}{i+1/2} = see notes
            %        for expression in terms of \rho_{i} and \rho_{i+1}
            %uH    - 1 x nVol-1 array of Roe adjusted velocity at cell
            %        interfaces.  uH(i) = \hat{u}{i+1/2} = see notes
            %        for expression in terms of u_{i} and \u_{i+1}
            %cH    - 1 x nVol-1 array of Roe adjusted sound speed at cell
            %        interfaces.  c(i) = \hat{c}{i+1/2} = see notes
            %        for expression in terms of c_{i} and c_{i+1}
            %drhoH - nVol-1 x 2 x 2 array of derivatives of Roe adjusted
            %        densities at cell interface.  The first index
            %        corresponds to the cell interface number, the second
            %        index corresponds to whether the derivative is with
            %        respect to the quantity at the cell left (=1) of the
            %        interface or right (=2) of it, and the third index
            %        corresponds to the variable with which the derivative
            %        is taken with respect to.  There are only 2 pages in
            %        the third dimension b/c drhoHde = 0.  drhoH(i,:,1) = 
            %        [d(\rho_{i+1/2})/d(\rho_i),
            %        d(\rho_{i+1/2})/d(\rho_{i+1})].   Changing the third
            %        dimension to a 2 takes all derivatives w.r.t. rho*u
            %duH   - nVol-1 x 2 x 3 array of derivatives of Roe adjusted
            %        velocities at cell interface. Explanation similar to
            %        drhoH case.  
            %dcH   - nVol-1 x 2 x 3 array of derivatives of Roe adjusted
            %        sound speeds at cell interface.Explanation similar to
            %        drhoH case.  
            %--------------------------------------------------------------
            
            %To easily extend this to GNAT, pre-compute an index array that
            %will allow you to write:
            %rhoH=sqrt(rho(indA1)).*sqrt(rho(indA2));
            
            %See notes for straightforward (but tedious) derivative
            %computations.
            dPdrho = 0.5*(obj.gamma-1)*u.*u;
            dPdrhou= -(obj.gamma-1)*u;
            
            rhoH=sqrt(rho(1:end-1)).*sqrt(rho(2:end));
            rhoGAvg=sqrt(rho(1:end-1))+sqrt(rho(2:end));
            uH=(sqrt(rho(1:end-1)).*u(1:end-1) + sqrt(rho(2:end)).*u(2:end))./rhoGAvg;
            hH=((e(1:end-1)+P(1:end-1))./sqrt(rho(1:end-1)) + (e(2:end)+P(2:end))./sqrt(rho(2:end)))./rhoGAvg;
            cH=sqrt((obj.gamma-1)*(hH-0.5*uH.^2));
            
            drhoH = zeros(obj.nVol-1,2); %only 1 pages bc derivatives wrt rho*u and e are zero
            drhoH(:,1) = 0.5*(1./sqrt(rho(1:end-1))).*sqrt(rho(2:end));
            drhoH(:,2) = 0.5*(1./sqrt(rho(2:end))).*sqrt(rho(1:end-1));
            
            irhoSum1 = 1./(rhoH+rho(1:end-1));
            irhoSum2 = 1./(rho(2:end)+rhoH);
            
            duH = zeros(obj.nVol-1,2,3); %only 2 pages bc independent of e
            duH(:,1,1)=-0.5*irhoSum1.*(uH+u(1:end-1));
            duH(:,2,1)=-0.5*irhoSum2.*(u(2:end)+uH);
            duH(:,1,2)=irhoSum1;
            duH(:,2,2)=irhoSum2;
            
            dhH = zeros(obj.nVol-1,2,3);
            dhH(:,1,1) = irhoSum1.*(-0.5*hH - 0.5*(e(1:end-1)+P(1:end-1))./rho(1:end-1) + dPdrho(1:end-1));
            dhH(:,2,1) = irhoSum2.*(-0.5*hH - 0.5*(e(2:end)+P(2:end))./rho(2:end) + dPdrho(2:end));
            dhH(:,1,2) = dPdrhou(1:end-1).*irhoSum1;
            dhH(:,2,2) = dPdrhou(2:end).*irhoSum2;
            dhH(:,1,3) = obj.gamma.*irhoSum1;
            dhH(:,2,3) = obj.gamma.*irhoSum2;
            
            dcH = zeros(obj.nVol-1,2,3);
            dcH(:,1,1) = (obj.gamma-1)*(0.5./cH).*(dhH(:,1,1)' - uH.*duH(:,1,1)');
            dcH(:,2,1) = (obj.gamma-1)*(0.5./cH).*(dhH(:,2,1)' - uH.*duH(:,2,1)');
            
            dcH(:,1,2) = (obj.gamma-1)*(0.5./cH).*(dhH(:,1,2)' - uH.*duH(:,1,2)');
            dcH(:,2,2) = (obj.gamma-1)*(0.5./cH).*(dhH(:,2,2)' - uH.*duH(:,2,2)');
            
            dcH(:,1,3) = (obj.gamma-1)*(0.5./cH).*dhH(:,1,3)';
            dcH(:,2,3) = (obj.gamma-1)*(0.5./cH).*dhH(:,2,3)';
            
        end 
        
        function  [F] = eulerNonlin(obj,rho,u,e,P)
            %This function computes the nonlinear Euler flux term.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %rho, u, e, P - 1 x N array of densities, velocities, energies,
            %               and pressures, respectively.
            %
            %Outputs:
            %--------
            %F - 3 x N array of euler flux terms.  Each column corresponds
            %    to the flux at rho(i), u(i), e(i), P(i)
            %--------------------------------------------------------------
            
                F = [rho.*u; rho.*u.*u+P; (e+P).*u];
                
        end
        
        function  [S,Si,CA,CAi,Lam,dS,dSi,dCA,dCAi] = eulerJac(obj,rho,u,c)
            %This function computes the nonlinear Euler flux Jacobian term.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %rho, u, e, P - scalar density, velocity, energy, and pressure
            %               (could be at cell-center or cell interface)
            %
            %Outputs:
            %--------
            %S,Si,CA,CAi,Lam - 3 x 3 matrices defining the Euler flux
            %                  Jacobian corresponding to the density,
            %                  velocity, energy, and pressure given.
            %dS,dSi,dCA,dCAi,dLam - 3 x 3 x 3 matrices defining the
            %                       derivative of each matrix above w.r.t.
            %                       rho, rho*u, and e.  I.E. dS(i,j,1) =
            %                       d(S_{i,j})/d(rho)
            %--------------------------------------------------------------
 
            
            %Each matrix has size 3x3 and the derivatives have size 3x3x3.
            %rho, u, c assumed scalar.
            
            %Compute S, Si, CA, CAi, Lam according to their definitions,
            %see notes.
            alpha=0.5*u.^2;
                        
            S=[1,0,0;-u/rho,1/rho,0;alpha*obj.beta,-u*obj.beta,obj.beta];
            Si=[1,0,0;u,rho,0;alpha,rho*u,(1/obj.beta)];
            CA=[1,0,-1/(c*c);0,rho*c,1;0,-rho*c,1];
            CAi=[1,0.5/(c*c),0.5/(c*c);0,0.5/(rho*c),-0.5/(rho*c);0,0.5,0.5];
            Lam = [u,0,0;0,u+c,0;0,0,u-c];
            
            %If derivatives of the matricies requested, compute them, see
            %notes.
            dS=zeros(3,3,2);
            dSi=zeros(3,3,2);
            dCA=zeros(3,3,2);
            dCAi=zeros(3,3,2);
            if nargout>5
                %dS% First page = dSdrho, Second page = dSdu
                dS(2,1:2,1)=[u/(rho*rho),-1/(rho*rho)];
                dS(2:3,1:2,2)=[-1/rho,0;u*obj.beta,-obj.beta];
                
                %dSi$ First page = dSidrho, Second page = dSidu
                dSi(2:3,2,1)=[1;u];
                dSi(2:3,1:2,2)=[1,0;u,rho];
                
                %dCA% First page = dCAdrho, Second page = dCAdc
                dCA(2:3,2,1)=[c;-c];
                dCA(:,2:3,2)=[0,2/(c*c*c);rho,0;-rho,0];
                
                %dCAi% First page = dCAidrho, Second page = dCAidc
                dCAi(2,2:3,1) = [-0.5/(rho*rho*c),0.5/(rho*rho*c)];
                dCAi(1:2,2:3,2) = [-1/(c*c*c),-1/(c*c*c);-0.5/(rho*c*c),0.5/(rho*c*c)];
            end
        end
        
        function  [Q,dQ,dQdS] = forceTerm(obj,u,P)
            %This function computes the forcing term in the quasi-1D Euler
            %equations. 
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %P - 1 x nVol array of pressures.
            %
            %Outputs:
            %--------
            %Q - 3 x nVol array of euler source terms.  Each column
            %    corresponds to the source at P(i)
            %--------------------------------------------------------------
           
            %Compute derivatives of pressure w.r.t. conservative variables.
            dPdrho = 0.5*(obj.gamma-1)*u.*u;
            dPdrhou= -(obj.gamma-1)*u;
            dPde=obj.gamma-1;
            
            %Compute Q and its derivative w.r.t. conservative variables.
            %See notes.
            Q = [zeros(1,obj.nVol); (P./obj.SVol).*(obj.S(2:end)-obj.S(1:end-1))./obj.dx; zeros(1,obj.nVol)];
            dQ = zeros(3,3,obj.nVol); %First dimension is component (since Q is 3x1 vector), 
                                      %second dimension is the index the derivative is taken w.r.t., 
                                      %the third is the volume
            dQ(2,1,:)=dPdrho.*(obj.S(2:end)-obj.S(1:end-1))./(obj.SVol.*obj.dx);
            dQ(2,2,:)=dPdrhou.*(obj.S(2:end)-obj.S(1:end-1))./(obj.SVol.*obj.dx);
            dQ(2,3,:)=dPde.*(obj.S(2:end)-obj.S(1:end-1))./(obj.SVol.*obj.dx);
            
            dQdS=[];
            if nargout==3
                dQdS=zeros(1,3,obj.nVol); %First dimension is component (since Q is 3x1 vector) -> since the first and third row of Q are zero, we only store the 2nd row, 
                                          %second dimension is the index the derivative is taken w.r.t. (1 -> dS_{i-1/2}, 2 -> dS_{i}, 3 -> dS_{i+1/2}), 
                                          %the third is the volume
                dQdS(1,1,:) = -P./(obj.SVol.*obj.dx);
                dQdS(1,2,:) = -P.*(obj.S(2:end)-obj.S(1:end-1))./(obj.SVol.*obj.SVol.*obj.dx);
                dQdS(1,3,:) =  P./(obj.SVol.*obj.dx);
            end
        end
        
        function  [R1,J1,Cpp,dCpp,dR1dp] = fullImplicitInletBC(obj,rho,u,P,c,dc)
            %This function computes the residual and Jacobian contributions
            %of the fully implicit inlet boundary condition:
            %    [ P - Pt*(1-(gamma-1)/(gamma+1)*u^2/astar2)^(gamma/(gamma-1));...
            %      rho - P/((gamma-1)*Cv*T(u));...
            %      dPdt - rho*c*dudt+(u-c)*(dPdx - rho*c*dudx) - gamma*P*(u/S)*dSdx] = 0
            %which can be written:
            % Cpp*dVdt + [ P - Pt*(1-(gamma-1)/(gamma+1)*u^2/astar2)^(gamma/(gamma-1));...
            %              rho - P/((gamma-1)*Cv*T(u));...
            %              (u-c)*(dPdx - rho*c*dudx) -  gamma*P*(u/S)*dSdx] = 0
            %where
            %  Cpp = [0,0,0;0,0,0;0,-rho*c,1]
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %rho - 1x2 array of densities at volumes 1 and 2
            %u   - 1x2 array of velocities at volumes 1 and 2
            %P   - 1x2 array of pressures at volumes 1 and 2
            %c   - 1x2 array of sound speeds at volumes 1 and 2
            %dc  - 3x2 array of derivatives of sound speeds at volumes 1
            %      and 2 w.r.t. primitive variables:
            %      dc = [ dcdrho_1, dcdrho_2;...
            %             dcdu_1  , dcdu_2  ;...
            %             dcdP_1  , dcdP_2  ]
            %
            %Outputs:
            %--------
            %R1    - nonlinear inlet residual (without time derivative terms)
            %J1    - nonlinear inlet jacobian
            %Cpp  - matrix multiplying inlet time derivative term
            %dCpp - jacobian of Cpp
            %--------------------------------------------------------------
            
            %Compute kappa and its derivative for readability later
            kappa  = 1 - ((obj.gamma-1)/(obj.gamma+1))*u(1)^2/obj.astar2;
            dkappa = -2*((obj.gamma-1)/(obj.gamma+1))*u(1)/obj.astar2;
            
            %Compute temperature and its derivative.
            T    = obj.Tt*kappa;
            dTdu = obj.Tt*dkappa;
            
            %Compute the cell center to cell center distance (between cell
            %1 and 2)
            dxBar=0.5*(obj.dx(2)+obj.dx(1));
            
            %Compute residual from definition (notes pg. 7)
            R1=zeros(3,1);
            R1(1) = rho(1) - P(1)/((obj.gamma-1)*obj.Cv*T);
            R1(2) = P(1) - obj.Pt*kappa^(obj.gamma/(obj.gamma-1));
            R1(3) = (u(1)-c(1))*(P(2)-P(1) - rho(1)*c(1)*(u(2)-u(1)))/dxBar + ... 
                obj.gamma*P(1)*u(1)*(obj.S(2)-obj.S(1))/(obj.SVol(1)*obj.dx(1));
      
            if nargout == 1
                J1=[]; Cpp=[]; dCpp=[]; dR1dp=[]; 
                return;
            end
            
            %Compute Cpp and derivative from its definition.  See notes
            %pg.7
            Cpp=[0,0,0;0,0,0;0,-rho(1)*c(1),1];
            dCpp=zeros(3,3,3);
            dCpp(3,2,1)= -(c(1)+rho(1)*dc(1,1));
            dCpp(3,2,2)= -rho(1)*dc(1,2);
            dCpp(3,2,3)= -rho(1)*dc(1,3);
            
            %Fill in Jacobian contributions (notes pg. 7)
            J1=zeros(3,6);
            J1(1,1)=1;
            J1(1,2) = P(1)/((obj.gamma-1)*obj.Cv*T*T)*dTdu;
            J1(1,3) = -1/((obj.gamma-1)*obj.Cv*T);
            J1(2,2) = -obj.Pt*(obj.gamma/(obj.gamma-1))*kappa^(1/(obj.gamma-1))*dkappa;
            J1(2,3)=1;
            J1(3,1)=-(u(1)-c(1))*(c(1)+rho(1)*dc(1,1))*(u(2)-u(1))/dxBar - dc(1,1)*(P(2)-P(1)-rho(1)*c(1)*(u(2)-u(1)))/dxBar;
            J1(3,2)=(u(1)-c(1))*(-rho(1)*dc(1,2)*(u(2)-u(1))+rho(1)*c(1))/dxBar + (1-dc(1,2))*(P(2)-P(1)-rho(1)*c(1)*(u(2)-u(1)))/dxBar + obj.gamma*P(1)*(obj.S(2)-obj.S(1))/(obj.dx(1)*obj.SVol(1));
            J1(3,3)=-(u(1)-c(1))*(1+rho(1)*dc(1,3)*(u(2)-u(1)))/dxBar + (-dc(1,3))*(P(2)-P(1)-rho(1)*c(1)*(u(2)-u(1)))/dxBar + obj.gamma*u(1)*(obj.S(2)-obj.S(1))/(obj.dx(1)*obj.SVol(1));
            
            %Derivatives wrt conservative variables (notes pg. 8)
            J1(3,4)=(u(1)-c(1))*(0.5*(obj.gamma-1)*u(2)^2 + rho(1)*c(1)*u(2)/rho(2))/dxBar;
            J1(3,5)=(u(1)-c(1))*(-(obj.gamma-1)*u(2) - rho(1)*c(1)/rho(2))/dxBar;
            J1(3,6)=(u(1)-c(1))*(obj.gamma-1)/dxBar;
        
            dR1dp=[];
            if nargout == 5
                dR1dp = [zeros(2,3);...
                    -obj.gamma*P(1)*u(1)/(obj.SVol(1)*obj.dx(1)),...
                    -obj.gamma*P(1)*u(1)*(obj.S(2)-obj.S(1))/(obj.SVol(1)^2*obj.dx(1)),...
                    obj.gamma*P(1)*u(1)/(obj.SVol(1)*obj.dx(1))]*obj.dSdp(1:3,:);
            end
        end
        
        function  [R2,J2,dR2dp] = governEqn(obj,rho,u,P,c,e,dc)
            %This function computes the residual and Jacobian contributions
            %of the governing Quasi-1D Euler equations:
            %    dUdt + (1/S)*d(FS)dx = Q,  
            % where U = vector of conservative variables
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %rho - 1xnVol array of densities at each volume
            %u   - 1xnVol array of velocities at each volume
            %P   - 1xnVol array of pressures at each volume
            %c   - 1xnVol array of sound speeds at each volume
            %dc  - 3xnVol array of derivatives of sound speeds at each
            %      volume w.r.t. conservative variables:
            %      dc = [ dcdrho_1 , ..., dcdrho_nVol ;...
            %             dcdrhou_1, ..., dcdrhou_nVol;...
            %             dcde_1   , ..., dcde_nVol   ]
            %
            %Outputs:
            %--------
            %R2    - nonlinear residual (without time derivative terms)
            %J2    - nonlinear jacobian
            %--------------------------------------------------------------
            
            %Compute Roe Flux (MacCormack notes)
            [roeF,droeF] = obj.roeFlux(rho,u,P,c,e,dc); %F_{i+1/2} for i=1,...,nVol-1
            
            %ADJUST DERIVATIVES OF 1st and LAST VOLUME TO TAKE INTO ACCOUNT
            %THAT FIRST THREE AND LAST THREE VARIABLES ARE PRIMITIVE
            %VARIABLES (see notes page 17)
            %Compute the change of variables Jacobian between conservative
            %and primitive variables at the inflow.
            dUdV=[1,0,0;u(1),rho(1),0;0.5*u(1)*u(1),rho(1)*u(1),1/(obj.gamma-1)];
            %Use change of variables to ensure derivatives w.r.t. first
            %three variables are w.r.t. PRIMITIVE
            droeF(:,1:3,1)=droeF(:,1:3,1)*dUdV;
            %Compute the change of variables Jacobian between conservative
            %and primitive variables at the outflow.
            dUdV=[1,0,0;u(end),rho(end),0;0.5*u(end)*u(end),rho(end)*u(end),1/(obj.gamma-1)];
            %Use change of variables to ensure derivatives w.r.t. first
            %three variables are w.r.t. PRIMITIVE
            droeF(:,4:6,end)=droeF(:,4:6,end)*dUdV;
            
            %Get forcing term at each volume
            if nargout == 3
                [Q,dQ,dQdS] = obj.forceTerm(u,P);
            else
                [Q,dQ] = obj.forceTerm(u,P);
            end
            
            
           right = bsxfun(@rdivide,bsxfun(@times,roeF(:,2:end),obj.S(3:end-1)),(obj.SVol(2:end-1).*obj.dx(2:end-1)));
	   left=  bsxfun(@rdivide,-bsxfun(@times,roeF(:,1:end-1),obj.S(2:end-2)),(obj.SVol(2:end-1).*obj.dx(2:end-1)));

           
 
            %Compute residual
            R2 = bsxfun(@rdivide,(bsxfun(@times,roeF(:,2:end),obj.S(3:end-1)) ...
                               -  bsxfun(@times,roeF(:,1:end-1),obj.S(2:end-2))),...
                (obj.SVol(2:end-1).*obj.dx(2:end-1))) - Q(:,2:end-1);

% KTC: -  bsxfun(@times,roeF(:,1:end-1),obj.S(2:end-2): from first cell (divide by (obj.SVol(2:end-1).*obj.dx(2:end-1)))
%- bsxfun(@times,roeF(:,2:end),obj.S(3:end-1)): from last cell (divide by (obj.SVol(2:end-1).*obj.dx(2:end-1)))            
%- Q(:,2:end-1): from every cell
%- Also add U^n+1 - U^n from every cell

            if nargout == 1
                J2=[]; dR2dp=[];
                return;
            end
            
            %Allocate space for Jacobian (may be faster to directly fill
            %sparse arrays and use J2=sparse(i,j,s,m,n))
            J2=spalloc(3*(obj.nVol-2),3*obj.nVol,3*9*(obj.nVol-2));
            
            %Fill entries of Jacobian according to:
            %J2(3*(i-1)+1:3*i,3*(i-1)+1:3*i)     = dR2_{i}/dU_{i-1}
            %J2(3*(i-1)+1:3*i,3*i+1:3*(i+1))     = dR2_{i}/dU_{i}
            %J2(3*(i-1)+1:3*i,3*(i+1)+1:3*(i+2)) =  dR2_{i}/dU_{i+1}
            for k = 2:obj.nVol-1
                i=k-1;
                J2(3*(i-1)+1:3*i,3*(i-1)+1:3*(i+2)) = [-(1./(obj.SVol(k)*obj.dx(k)))*obj.S(k)*droeF(:,1:3,k-1),...
                                                       (1./(obj.SVol(k)*obj.dx(k)))*(obj.S(k+1)*droeF(:,1:3,k) - obj.S(k)*droeF(:,4:6,k-1))-dQ(:,:,k),...
                                                       (1./(obj.SVol(k)*obj.dx(k)))*obj.S(k+1)*droeF(:,4:6,k)];
            end
            
            dR2dp=[];
            if nargout == 3
                n=length(obj.p);
                dR2dp=zeros(3*obj.nVol-6,n);
                
                for k = 2:obj.nVol-1
                    i=k-1;
                    
                    ind = ((3*(i-1)+1-(i-1)):3*i-(i-1))+2;
                    dR2dp(3*(i-1)+1:3*i,:) = ([-(1/(obj.SVol(k)*obj.dx(k)))*roeF(:,k-1),-(1/(obj.SVol(k)^2*obj.dx(k)))*(roeF(:,k)*obj.S(k+1)-roeF(:,k-1)*obj.S(k)),(1/(obj.SVol(k)*obj.dx(k)))*roeF(:,k)] ...
                                   - [zeros(1,3);dQdS(1,:,k);zeros(1,3)])*obj.dSdp(ind,:);
                end
            end
        end
        
        function  [R3,J3,Cp,dCp,dR3dp] = fullImplicitOutletBC(obj,rho,u,P,c,dc)
            %This function computes the residual and Jacobian contributions
            %of the fully implicit outlet boundary condition:
            %    Cp*dVdt = -Lambda*Cp*dVdx + Q'  
            %
            % where V = vector of primitive variables
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %rho - 1x2 array of densities at volumes nVol-1 and nVol
            %u   - 1x2 array of velocities at volumes nVol-1 and nVol
            %P   - 1x2 array of pressures at volumes nVol-1 and nVol
            %c   - 1x2 array of sound speeds at volumes nVol-1 and nVol
            %dc  - 3x2 array of derivatives of sound speeds at volumes
            %      nVol-1 and nVol w.r.t. primitive variables:
            %      dc = [ dcdrho_1, dcdrho_2;...
            %             dcdu_1  , dcdu_2  ;...
            %             dcdP_1  , dcdP_2  ]
            %
            %Outputs:
            %--------
            %R3    - nonlinear outlet residual (without time derivative terms)
            %J3    - nonlinear outlet jacobian
            %Cp  - matrix multiplying outlet time derivative term
            %dCp - jacobian of Cp
            %--------------------------------------------------------------

            
            %All input variables are 2x1 vectors, with the 1,1 entry
            %corresponding to the value in cell nVol-1 ad the 2,1 entry
            %corresponding to the value in cell nVol
            
            %MIGHT NEED TO SMOOTH g for optimization!
            g = (u(2) >= c(2));
            %Compute Lambda, Cp, Qp from definition (MacCormack notes: Ch.13,
            %pg.10) along with their derivatives (see notes pgs. 8-9)
            Lam=[u(2),0,0;0,u(2)+c(2),0;0,0,(u(2)-c(2))*g];
            dLam(:,:,1)=diag([0,dc(2,1),-g*dc(2,1)]);
            dLam(:,:,2)=diag([1,1+dc(2,2),g*(1-dc(2,2))]);
            dLam(:,:,3)=diag([0,dc(2,3),-g*dc(2,3)]);
            
            Cp = [1,0,-1/(c(2)*c(2));0,rho(2)*c(2),1;0,-rho(2)*c(2)*g,1];
                       
            q = -obj.gamma*P(2)*u(2)/obj.SVol(end)*(obj.S(end)-obj.S(end-1))/obj.dx(end);
            dq = [0, -obj.gamma*P(2)/obj.SVol(end)*(obj.S(end)-obj.S(end-1))/obj.dx(end), -obj.gamma*u(2)/obj.SVol(end)*(obj.S(end)-obj.S(end-1))/obj.dx(end)];
                
            Qp = [0;q;g*q];
            dQp = [0,0,0; dq; g*dq];
            
            %Compute cell center to cell center distance between cells
            %nVol-1 and nVol
            dxBar=0.5*(obj.dx(end)+obj.dx(end-1));
            %Form primitive variables and precompute various quantities.
            V = [rho;u;P];
            dVdx = (V(:,2)-V(:,1))/dxBar;
            CpdVdx = Cp*dVdx;
            
            dVdU1=[1,0,0;-u(1)/rho(1),1/rho(1),0;0.5*(obj.gamma-1)*u(1)*u(1),-(obj.gamma-1)*u(1),obj.gamma-1];
            
            %Compute outlet residual and Jacobian (notes pg. 8-9)
            R3 = Lam*CpdVdx - Qp;
            
            if nargout == 1
                J3=[]; Cp=[]; dCp=[]; dR3dp=[];
                return;
            end
            
            dCp(:,:,1)= [0,0,2/(c(2)^3)*dc(2,1);0,c(2)+rho(2)*dc(2,1),0;0,-g*(c(2)+rho(2)*dc(2,1)),0];
            dCp(:,:,2)= [0,0,2/(c(2)^3)*dc(2,2);0,rho(2)*dc(2,2),0;0,-rho(2)*g*dc(2,2),0];
            dCp(:,:,3)= [0,0,2/(c(2)^3)*dc(2,3);0,rho(2)*dc(2,3),0;0,-rho(2)*g*dc(2,3),0];

            J3 = [dLam(:,:,1)*CpdVdx + Lam*dCp(:,:,1)*dVdx, ...
                  dLam(:,:,2)*CpdVdx + Lam*dCp(:,:,2)*dVdx, ...
                  dLam(:,:,3)*CpdVdx + Lam*dCp(:,:,3)*dVdx];
            J3 = J3 + Lam*Cp/dxBar - dQp;
            J3 = [-Lam*Cp*dVdU1/dxBar,J3];
            
            dR3dp=[];
            if nargout == 5
                dR3dp = -[zeros(1,3);...
                         obj.gamma*P(2)*u(2)/(obj.SVol(end)*obj.dx(end)), obj.gamma*P(2)*u(2)*(obj.S(end)-obj.S(end-1))/(obj.SVol(end)^2*obj.dx(end)) ,-obj.gamma*P(2)*u(2)/(obj.SVol(end)*obj.dx(end));...
                         g*obj.gamma*P(2)*u(2)/(obj.SVol(end)*obj.dx(end)), g*obj.gamma*P(2)*u(2)*(obj.S(end)-obj.S(end-1))/(obj.SVol(end)^2*obj.dx(end)) ,-g*obj.gamma*P(2)*u(2)/(obj.SVol(end)*obj.dx(end))];
                dR3dp = dR3dp*obj.dSdp(end-2:end,:);
            end
        end
        
        function  [cfl] = computeCFL(obj,U)
            
            %Reshape U from a stacked vector ([rho_1;rho_1*u_1;e_1;...]) of
            %size 3nVol x 1 to a matrix of size 3 x nVol where U(1,:) =
            %rho, U(2,:) = rho.*u, and U(3,:) = e
            U = reshape(U,3,obj.nVol);
            %Extract primitive variables
            [~,u,~,c] = obj.conservativeToPrimitive(U(:,2:end-1));
            c=[sqrt(obj.gamma*U(3,1)/U(1,1)),c,sqrt(obj.gamma*U(3,end)/U(1,end))];
            u   = [U(2,1),u,U(2,end)];
            cfl = min(obj.dx./(u+c));
            
        end
                
        function  [] = setPrevSV(obj,~)
            return;
        end
        
        function  [] = setProperty(obj,prop,val)
            obj.(prop)=val;
        end

        %GNAT Functions
        function  [] = JacobianStructure(obj)
            %This function computes a logical matrix of the Jacobian
            %indicating the nonzero entries.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - quasi1dEuler object
            %
            %Outputs:
            %--------
            %No outputs.
            %--------------------------------------------------------------
            
            obj.Jstruct = spalloc(3*obj.nVol,3*obj.nVol,3*3*(2*2 + 3*(obj.nVol-2)));
            
            %First three equations (depend on first two volumes = first 6
            %variables)
            obj.Jstruct(1:3,1:6)=1;
            %All other equations
            for k = 2:obj.nVol-1
                i=k-1;
                obj.Jstruct(3*(k-1)+1:3*k,3*(i-1)+1:3*(i+2))=1;
            end
            %Last three equations (depend on last two volumes = last 6
            %variables)
            obj.Jstruct(end-2:end,end-5:end)=1;
            
        end
        
        function  [ind] = node2ind(obj,nodenum)
            %This function takes a node in the mesh and maps it to an index
            %in the full state vector
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - quasi1dEuler object
            %nodenum - scalar specifying node numbers
            %
            %Outputs:
            %--------
            %ind     - vector specifying the corresponding indices in the
            %          full state vector
            %--------------------------------------------------------------
            
            ind = 3*(nodenum-1)+1:3*nodenum;
        end
        
        function  [R,J,Cp,Cpp,dCp,dCpp] = ResJacGNAT(obj,U,~)
            %This function computes the residual and jacobian of the
            %quasi-1D Euler equations on the mask.  Also, terms needed for
            %the time stepping residual and jacobian are computed.
            %Governing equations are:
            %
            %Inlet (only if inlet included in mask): 
            %    [ P - Pt*(1-(gamma-1)/(gamma+1)*u^2/astar2)^(gamma/(gamma-1));...
            %      rho - P/((gamma-1)*Cv*T(u));...
            %      dPdt - rho*c*dudt+(u-c)*(dPdx - rho*c*dudx) - gamma*P*(u/S)*dSdx] = 0
            %which can be written:
            % Cpp*dVdt + [ P - Pt*(1-(gamma-1)/(gamma+1)*u^2/astar2)^(gamma/(gamma-1));...
            %              rho - P/((gamma-1)*Cv*T(u));...
            %              (u-c)*(dPdx - rho*c*dudx) -  gamma*P*(u/S)*dSdx] = 0
            %where
            %  Cpp = [0,0,0;0,0,0;0,-rho*c,1]
            %
            %Governing:
            %    dUdt + (1/S)*d(FS)dx = Q
            %
            %Outlet (only if outlet included in mask):
            %    Cp*dVdt = -Lambda*Cp*dVdx + Q'  
            %
            % where V = vector of primitive variables, U = vector of
            % conservative variables.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj - quasiEuler1D object
            %U   - state vector: U(1:3) = primitive variables at inlet,
            %      U(end-2:end) = primitive variables at outlet, and
            %      U(4:end-3) = conservative variables on domain interior.
            %
            %Outputs:
            %--------
            %R    - nonlinear residual (without time derivative terms)
            %J    - nonlinear jacobian
            %Cp   - matrix multiplying outlet time derivative term
            %Cpp  - matrix multiplying inlet time derivative term
            %dCp  - jacobian of Cp
            %dCpp - jacobian of Cpp
            %--------------------------------------------------------------
           
            %Reshape U from a stacked vector ([rho_1;rho_1*u_1;e_1;...]) of
            %size 3nVol x 1 to a matrix of size 3 x nVol where U(1,:) =
            %rho, U(2,:) = rho.*u, and U(3,:) = e.  Except first three and
            %last three values are primitive variables while others are
            %conservative.
% 	keyboard
            U = reshape(U,3,obj.nVolSamp);
            %Extract primitive variables.  Take into account that the first
            %three and last three variables are PRIMITIVE VARIABLES while
            %all other are conservative variables.
            startInd = 1 + obj.ind1;
            endInd   = obj.nVolSamp - obj.indN;
            [rho,u,P,c] = obj.conservativeToPrimitive(U(:,startInd:endInd));
            e=U(3,startInd:endInd);
            if obj.ind1
                rho = [U(1,1),rho]; %Density
                u   = [U(2,1),u]; %Velocity
                P   = [U(3,1),P]; %Pressure
                c   = [sqrt(obj.gamma*P(1)/rho(1)),c]; %Speed of sound
                e   = [P(1)/(obj.gamma-1)+rho(1)*u(1)^2/2,e]; %Energy
            end
            if obj.indN
                rho = [rho,U(1,end)]; %Density
                u   = [u,U(2,end)]; %Velocity
                P   = [P,U(3,end)]; %Pressure
                c   = [c,sqrt(obj.gamma*P(end)/rho(end))]; %Speed of sound
                e   = [e,P(end)/(obj.gamma-1)+rho(end)*u(end)^2/2]; %Energy
            end
            %Derivative of sound speed w.r.t. primitive variables.
            %dc_prim = [dcdrho_1, ..., dcdrho_nVol; ...
            %           dcdu_1  , ..., dcdu_nVol; ...
            %           dcdP_1  , ..., dcdP_nVol]
            dc_prim = [-0.5*obj.gamma*P./(c.*rho.*rho);zeros(1,obj.nVolSamp);0.5*obj.gamma./(c.*rho)]';
            %Derivative of sound speed w.r.t. conservative variables.
            %dc_cons = [dcdrho_1 , ..., dcdrho_nVol; ...
            %           dcdrhou_1, ..., dcdrhou_nVol; ...
            %           dcde_1  , ..., dcde_nVol]
            dc_cons = [0.5*obj.gamma./(c.*rho).*(0.5*(obj.gamma-1)*u.*u - P./rho);...
                      -0.5*obj.gamma*(obj.gamma-1)*u./(rho.*c);...
                       0.5*obj.gamma*(obj.gamma-1)./(rho.*c)]';
            
            %Get residual and Jacobian contributions from inlet, governing
            %equations (interior), and outlet.  Also compute matrices
            %multiplying time derivatives.
            R1=[]; J1=[]; R3=[]; J3=[];
            Cp=[]; Cpp=[]; dCp=[]; dCpp=[];
            if obj.ind1
                [R1,J1,Cpp,dCpp]=obj.fullImplicitInletBC(rho(1:2),u(1:2),P(1:2),c(1:2),dc_prim(1:2,:));
            end
            [R2,J2]=obj.governEqnGNAT(rho,u,P,c,e,dc_cons); %R2 is a 3 x nVol-2
            if obj.indN
                [R3,J3,Cp,dCp]=obj.fullImplicitOutletBC(rho(end-1:end),u(end-1:end),P(end-1:end),c(end-1:end),dc_prim(end-1:end,:));
            end
            
            %Form residual and Jacobian.
            R=[R1;R2(:);R3];
            J1=J1'; J3=J3';
            %keyboard
            J=[J1(:);J2(:);J3(:)];
         end
        
        function  [roeF,droeF] = roeFluxGNAT(obj,rho,u,P,c,e,dc)
            %This function computes the Roe flux and Jacobian at each
            %inter-cell face in the mask (since there are nVolMask cells,
            %there are nVolMask-1  interfaces). 
            %Roe Flux at the interface between cell i and i+1:
            %F_{i+1/2} = 0.5*(F_{i} + F_{i+1}) - 0.5*|\hat{A}_{i+1/2}|*(U_{i+1}-U_{i})
            %for i = 1, ..., nVol-1 (need for each INTERIOR face)
            %This requires F_{i}, U_{i} for i = 1, ..., nVol
            %--------------------------------------------------------------
            %Input:
            %------
            %obj - quasi1dEuler object
            %rho - 1xnVolSamp array of densities in each volume
            %u   - 1xnVolSamp array of velocities in each volume
            %P   - 1xnVolSamp array of pressures in each volume
            %c   - 1xnVolSamp array of sound speeds in each volume
            %e   - 1xnVolSamp array of energies in each volume
            %dc  - 3xnVolSamp array of derivatives of sound speed w.r.t.
            %      CONSERVATIVE VARIABLES in each volume.
            %      dc = [dcdrho_1 , ..., dcdrho_nVolSamp; ...
            %            dcdrhou_1, ..., dcdrhou_nVolSamp; ...
            %            dcde_1  , ..., dcde_nVolSamp]
            %
            %Output:
            %-------
            %roeF  - 3x(nVolSamp-1) array of Roe flux at each interface.  The
            %        first dimension respresents the component of the flux
            %        (rho, rho*u, e) and the second dimension represents
            %        the interface number.
            %droeF - 3x6x(nVolSamp-1) array of Roe flux derivatives at each
            %        interface.  The first dimension represents the
            %        component of the flux (rho, rho*u, e), the second
            %        dimension represents the derivative taken w.r.t. and
            %        the third dimension represents the interface number.
            %        Since F_{i+1/2} only depends on U_{i} and U_{i+1}, the
            %        Jacobian of F_{i+1/2} only depends on the state at
            %        volumes i and i+1.  Therefore (in the second dimension
            %        of drhoF), instead of storing the derivative w.r.t.
            %        3*nVol variables (where most are zero), we store the
            %        derivative w.r.t. only 6 variables (rho_i, rhou_i,
            %        e_i, rho_{i+1}, rhou_{i+1}, e_{i+1}) 
            %--------------------------------------------------------------

            %Compute nonlinear term at each cell volume in mask
            F=obj.eulerNonlin(rho,u,e,P);
            %Compute the Roe adjusted values at each interface as well as
            %their derivatives w.r.t. state.  See rhoTerms for description
            %of each variable.
            [rhoH,uH,cH,drhoH,duH,dcH] = obj.roeTermsGNAT(rho,u,P,e);
     
            %Initialize roeF and droeF
            roeF = 0.5*(F(:,obj.iarray)+F(:,obj.iarray+1));
            droeF = zeros(3,6,obj.nVolSamp-1);
            
            cnt=0;
            for i = obj.iarray(:)'
                cnt=cnt+1;
                %Form conservative state variables at volume i and i+1
                U_ip1=[rho(i+1);rho(i+1)*u(i+1);e(i+1)];
                U_i  =[rho(i);rho(i)*u(i);e(i)];
                dU = U_ip1-U_i;
                
                %Compute terms for \hat{A} and its derivative w.r.t. state
                [S,Si,CA,CAi,Lam,dS,dSi,dCA,dCAi] = obj.eulerJac(rhoH(cnt),uH(cnt),cH(cnt));
                
                %Apply entropy correction to Lam
                [Lam(1,1),Lam(2,2),Lam(3,3),dlam_u,dlam_upc,dlam_umc]=obj.entropyCorrect(rho(i:i+1),u(i:i+1),c(i:i+1),dc(i:i+1,:),uH(cnt),cH(cnt),duH(cnt,:,:),dcH(cnt,:,:));
                sgn=sign(diag(Lam));
                
                %Derivative contributions from absolute value
                dlam_u=sgn(1)*dlam_u;
                dlam_upc=sgn(2)*dlam_upc;
                dlam_umc=sgn(3)*dlam_umc;
                
                %Form \hat{A}
                Ahat=(Si*CAi*abs(Lam)*CA*S);
                
                %Compute the derivative of \hat{A} w.r.t. state time dU (chain rule).
                %Data will be stored according to:
                %dAhatTimesDU(:,1,1) = dAhat/drho_{i}*dU
                %dAhatTimesDU(:,1,2) = dAhat/drho_{i+1}*dU
                %The middle index corresponds to the derivative that is
                %taken: 1 = rho, 2 = rho*u, 3 = e
                dAhatTimesDU = zeros(3,3,2);
                for j = 1:2 %j = 1 -> d/d()_i, j = 2 -> d/d()_{i+1}
                    %dAhat/drho * dU
                    dAhatTimesDU(:,1,j) = (dSi(:,:,1)*drhoH(cnt,j)+dSi(:,:,2)*duH(cnt,j,1))*(CAi*(abs(Lam)*(CA*(S*dU)))) + ...
                        Si*((dCAi(:,:,1)*drhoH(cnt,j)+dCAi(:,:,2)*dcH(cnt,j,1))*(abs(Lam)*(CA*(S*dU)))) + ...
                        Si*(CAi*(diag([dlam_u(j,1),dlam_upc(j,1),dlam_umc(j,1)])*(CA*(S*dU)))) + ...
                        Si*(CAi*(abs(Lam)*((dCA(:,:,1)*drhoH(cnt,j)+dCA(:,:,2)*dcH(cnt,j,1))*(S*dU)))) + ...
                        Si*(CAi*(abs(Lam)*(CA*((dS(:,:,1)*drhoH(cnt,j)+dS(:,:,2)*duH(cnt,j,1))*dU))));
                    %dAhat/drhou * dU
                    dAhatTimesDU(:,2,j) = (dSi(:,:,2)*duH(cnt,j,2))*(CAi*(abs(Lam)*(CA*(S*dU)))) + ...
                        Si*((dCAi(:,:,2)*dcH(cnt,j,2))*(abs(Lam)*(CA*(S*dU)))) + ...
                        Si*(CAi*(diag([dlam_u(j,2),dlam_upc(j,2),dlam_umc(j,2)])*(CA*(S*dU)))) + ...
                        Si*(CAi*(abs(Lam)*((dCA(:,:,2)*dcH(cnt,j,2))*(S*dU)))) + ...
                        Si*(CAi*(abs(Lam)*(CA*((dS(:,:,2)*duH(cnt,j,2))*dU))));
                    %dAhat/de * dU
                    dAhatTimesDU(:,3,j) = Si*((dCAi(:,:,2)*dcH(cnt,j,3))*(abs(Lam)*(CA*(S*dU)))) + ...
                        Si*(CAi*(diag([dlam_u(j,3),dlam_upc(j,3),dlam_umc(j,3)])*(CA*(S*dU)))) + ...
                        Si*(CAi*(abs(Lam)*((dCA(:,:,2)*dcH(cnt,j,3))*(S*dU))));
                end
                
                %Compute Roe flux between volume i and i+1
                roeF(:,cnt) = roeF(:,cnt) - 0.5*Ahat*(U_ip1-U_i);
                
                %Compute Jacobian of Roe flux between volumne i and i+1
                [S,Si,CA,CAi,Lam] = obj.eulerJac(rho(i),u(i),c(i));
                Jblk1 = 0.5*(Si*CAi*Lam*CA*S - dAhatTimesDU(:,:,1) + Ahat);
                
                [S,Si,CA,CAi,Lam] = obj.eulerJac(rho(i+1),u(i+1),c(i+1));
                Jblk2 = 0.5*(Si*CAi*Lam*CA*S - dAhatTimesDU(:,:,2) - Ahat);
                               
                droeF(:,:,cnt)=[Jblk1,Jblk2];
            end
        end
        
        function  [rhoH,uH,cH,drhoH,duH,dcH] = roeTermsGNAT(obj,rho,u,P,e)
            %This function returns the Roe averaged density, velocity, and
            %sound speed at the faces of mask, and their derivatives.
            %Note: even though all output arrays have some dimension
            %nVolSamp-1 this is not exactly true since there will be extra
            %entries whenever there is a jump in the elements stored (i.e.
            %if volumes 1, 2, 69, 70, 72 are in the mask, there are two
            %jumps: from 2 -> 69 and 70 -> 72).
            %--------------------------------------------------------------
            %Input:
            %------
            %obj          - quasi1dEuler object
            %rho, u, P, e - 1 x nVolSamp array of density, velocity,
            %               pressure, and energy, respectively
            %
            %Output:
            %-------
            %rhoH  - 1 x nVolSamp-1 array of Roe adjusted densities at cell
            %        interfaces.  rhoH(i) = \hat{\rho}{i+1/2} = see notes
            %        for expression in terms of \rho_{i} and \rho_{i+1}
            %uH    - 1 x nVolSamp-1 array of Roe adjusted velocity at cell
            %        interfaces.  uH(i) = \hat{u}{i+1/2} = see notes
            %        for expression in terms of u_{i} and \u_{i+1}
            %cH    - 1 x nVolSamp-1 array of Roe adjusted sound speed at cell
            %        interfaces.  c(i) = \hat{c}{i+1/2} = see notes
            %        for expression in terms of c_{i} and c_{i+1}
            %drhoH - nVolSamp-1 x 2 x 2 array of derivatives of Roe adjusted
            %        densities at cell interface.  The first index
            %        corresponds to the cell interface number, the second
            %        index corresponds to whether the derivative is with
            %        respect to the quantity at the cell left (=1) of the
            %        interface or right (=2) of it, and the third index
            %        corresponds to the variable with which the derivative
            %        is taken with respect to.  There are only 2 pages in
            %        the third dimension b/c drhoHde = 0.  drhoH(i,:,1) = 
            %        [d(\rho_{i+1/2})/d(\rho_i),
            %        d(\rho_{i+1/2})/d(\rho_{i+1})].   Changing the third
            %        dimension to a 2 takes all derivatives w.r.t. rho*u
            %duH   - nVolSamp-1 x 2 x 3 array of derivatives of Roe adjusted
            %        velocities at cell interface. Explanation similar to
            %        drhoH case.  
            %dcH   - nVolSamp-1 x 2 x 3 array of derivatives of Roe adjusted
            %        sound speeds at cell interface.Explanation similar to
            %        drhoH case.  
            %--------------------------------------------------------------
            
            %To easily extend this to GNAT, pre-compute an index array that
            %will allow you to write:
            %rhoH=sqrt(rho(indA1)).*sqrt(rho(indA2));
            
            %See notes for straightforward (but tedious) derivative
            %computations.
            dPdrho = 0.5*(obj.gamma-1)*u.*u;
            dPdrhou= -(obj.gamma-1)*u;
            
            rhoH=sqrt(rho(obj.iarray)).*sqrt(rho(obj.iarray+1));
            rhoGAvg=sqrt(rho(obj.iarray))+sqrt(rho(obj.iarray+1));
            uH=(sqrt(rho(obj.iarray)).*u(obj.iarray) + sqrt(rho(obj.iarray+1)).*u(obj.iarray+1))./rhoGAvg;
            hH=((e(obj.iarray)+P(obj.iarray))./sqrt(rho(obj.iarray)) + (e(obj.iarray+1)+P(obj.iarray+1))./sqrt(rho(obj.iarray+1)))./rhoGAvg;
            cH=sqrt((obj.gamma-1)*(hH-0.5*uH.^2));
                     
            drhoH = zeros(length(rhoH),2); %only 1 pages bc derivatives wrt rho*u and e are zero
            drhoH(:,1) = 0.5*(1./sqrt(rho(obj.iarray))).*sqrt(rho(obj.iarray+1));
            drhoH(:,2) = 0.5*(1./sqrt(rho(obj.iarray+1))).*sqrt(rho(obj.iarray));
            
            irhoSum1 = 1./(rhoH+rho(obj.iarray));
            irhoSum2 = 1./(rho(obj.iarray+1)+rhoH);
            
            duH = zeros(length(uH),2,3); %only 2 pages bc independent of e
            duH(:,1,1)=-0.5*irhoSum1.*(uH+u(obj.iarray));
            duH(:,2,1)=-0.5*irhoSum2.*(u(obj.iarray+1)+uH);
            duH(:,1,2)=irhoSum1;
            duH(:,2,2)=irhoSum2;
            
            dhH = zeros(length(hH),2,3);
            dhH(:,1,1) = irhoSum1.*(-0.5*hH - 0.5*(e(obj.iarray)+P(obj.iarray))./rho(obj.iarray) + dPdrho(obj.iarray));
            dhH(:,2,1) = irhoSum2.*(-0.5*hH - 0.5*(e(obj.iarray+1)+P(obj.iarray+1))./rho(obj.iarray+1) + dPdrho(obj.iarray+1));
            dhH(:,1,2) = dPdrhou(obj.iarray).*irhoSum1;
            dhH(:,2,2) = dPdrhou(obj.iarray+1).*irhoSum2;
            dhH(:,1,3) = obj.gamma.*irhoSum1;
            dhH(:,2,3) = obj.gamma.*irhoSum2;
            
            dcH = zeros(length(cH),2,3);
            dcH(:,1,1) = (obj.gamma-1)*(0.5./cH).*(dhH(:,1,1)' - uH.*duH(:,1,1)');
            dcH(:,2,1) = (obj.gamma-1)*(0.5./cH).*(dhH(:,2,1)' - uH.*duH(:,2,1)');
            
            dcH(:,1,2) = (obj.gamma-1)*(0.5./cH).*(dhH(:,1,2)' - uH.*duH(:,1,2)');
            dcH(:,2,2) = (obj.gamma-1)*(0.5./cH).*(dhH(:,2,2)' - uH.*duH(:,2,2)');
            
            dcH(:,1,3) = (obj.gamma-1)*(0.5./cH).*dhH(:,1,3)';
            dcH(:,2,3) = (obj.gamma-1)*(0.5./cH).*dhH(:,2,3)';
        end
        
        function  [Q,dQ,dQdS] = forceTermGNAT(obj,u,P)
            %This function computes the forcing term in the quasi-1D Euler
            %equations.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %P - 1 x nVolMask array of pressures.
            %
            %Outputs:
            %--------
            %Q  - 3 x nVolMask array of euler source terms.  Each column
            %     corresponds to the source at P(i)
            %dQ - 1x3xnVolMask array of derivative of euler source
            %     terms. The first index is the component (only the second
            %     componenet of dQ is retained for space efficiency), the
            %     second index indicates the derivative, the third index is
            %     the volume.
            %dQdS - 1x3xnVolMask array of derivative of euler source
            %     terms w.r.t S. The first index is the component (only the
            %     second componenet of dQdS is retained for space
            %     efficiency), the second index indicates the derivative,
            %     the third index is the volume.
            %--------------------------------------------------------------
            
            %Compute derivatives of pressure w.r.t. conservative variables
            %on volumes in mask.
            dPdrho = 0.5*(obj.gamma-1)*u(obj.iarrayMask).*u(obj.iarrayMask);
            dPdrhou= -(obj.gamma-1)*u(obj.iarrayMask);
            dPde=obj.gamma-1;
            
            %Compute Q and its derivative w.r.t. conservative variables.
            %See notes.
            Q = [zeros(1,obj.nVolMask); (P(obj.iarrayMask)./obj.SVol).*(obj.S(obj.iarrayFace+1)-obj.S(obj.iarrayFace))./obj.dx; zeros(1,obj.nVolMask)];
            dQ = zeros(1,3,obj.nVolMask); %First dimension is component (since Q is 3x1 vector),
            %second dimension is the index the derivative is taken w.r.t.,
            %the third is the volume
            dQ(1,1,:)=dPdrho.*(obj.S(obj.iarrayFace+1)-obj.S(obj.iarrayFace))./(obj.SVol.*obj.dx);
            dQ(1,2,:)=dPdrhou.*(obj.S(obj.iarrayFace+1)-obj.S(obj.iarrayFace))./(obj.SVol.*obj.dx);
            dQ(1,3,:)=dPde.*(obj.S(obj.iarrayFace+1)-obj.S(obj.iarrayFace))./(obj.SVol.*obj.dx);
            
            dQdS=[];
            if nargout==3
                dQdS=zeros(1,3,obj.nVolMask); %First dimension is component (since Q is 3x1 vector) -> since the first and third row of Q are zero, we only store the 2nd row,
                %second dimension is the index the derivative is taken w.r.t. (1 -> dS_{i-1/2}, 2 -> dS_{i}, 3 -> dS_{i+1/2}),
                %the third is the volume
                dQdS(1,1,:) = -P(obj.iarrayMask)./(obj.SVol.*obj.dx);
                dQdS(1,2,:) = -P(obj.iarrayMask).*(obj.S(obj.iarrayFace+1)-obj.S(obj.iarrayFace))./(obj.SVol.*obj.SVol.*obj.dx);
                dQdS(1,3,:) = P(obj.iarrayMask)./(obj.SVol.*obj.dx);
            end
        end
        
        function  [R2,J2,dR2dp] = governEqnGNAT(obj,rho,u,P,c,e,dc)
            %This function computes the residual and Jacobian contributions
            %of the governing Quasi-1D Euler equations:
            %    dUdt + (1/S)*d(FS)dx = Q,
            % where U = vector of conservative variables
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %rho - 1xnVol array of densities at each volume
            %u   - 1xnVol array of velocities at each volume
            %P   - 1xnVol array of pressures at each volume
            %c   - 1xnVol array of sound speeds at each volume
            %dc  - 3xnVol array of derivatives of sound speeds at each
            %      volume w.r.t. conservative variables:
            %      dc = [ dcdrho_1 , ..., dcdrho_nVol ;...
            %             dcdrhou_1, ..., dcdrhou_nVol;...
            %             dcde_1   , ..., dcde_nVol   ]
            %
            %Outputs:
            %--------
            %R2    - nonlinear residual (without time derivative terms)
            %J2    - nonlinear jacobian
            %--------------------------------------------------------------
            
            %Compute Roe Flux (MacCormack notes)
            [roeF,droeF] = obj.roeFluxGNAT(rho,u,P,c,e,dc); %F_{i+1/2} for i=1,...,nVol-1
            
            %ADJUST DERIVATIVES OF 1st and LAST VOLUME TO TAKE INTO ACCOUNT
            %THAT FIRST THREE AND LAST THREE VARIABLES ARE PRIMITIVE
            %VARIABLES (see notes page 17)
            if obj.ind1
                %Compute the change of variables Jacobian between conservative
                %and primitive variables at the inflow.
                dUdV=[1,0,0;u(1),rho(1),0;0.5*u(1)*u(1),rho(1)*u(1),1/(obj.gamma-1)];
                %Use change of variables to ensure derivatives w.r.t. first
                %three variables are w.r.t. PRIMITIVE
                droeF(:,1:3,1)=droeF(:,1:3,1)*dUdV;
            end
            if obj.indN
                %Compute the change of variables Jacobian between conservative
                %and primitive variables at the outflow.
                dUdV=[1,0,0;u(end),rho(end),0;0.5*u(end)*u(end),rho(end)*u(end),1/(obj.gamma-1)];
                %Use change of variables to ensure derivatives w.r.t. first
                %three variables are w.r.t. PRIMITIVE
                droeF(:,4:6,end)=droeF(:,4:6,end)*dUdV;
            end
            
            %Get forcing term at each volume
            if nargout == 3
                [Q,dQ,dQdS] = obj.forceTermGNAT(u,P);
            else
                [Q,dQ] = obj.forceTermGNAT(u,P);
            end
            
            %Compute residual
            R2 = bsxfun(@rdivide,(bsxfun(@times,roeF(:,obj.iarrayFaceI+1),obj.S(obj.iarrayFace(1+obj.ind1:end-obj.indN)+1)) ...
                -  bsxfun(@times,roeF(:,obj.iarrayFaceI),obj.S(obj.iarrayFace(1+obj.ind1:end-obj.indN)))),...
                (obj.SVol(1+obj.ind1:end-obj.indN).*obj.dx(1+obj.ind1:end-obj.indN))) - Q(:,1+obj.ind1:end-obj.indN);
         
            J2=zeros(3*9*(obj.nVolMask-obj.ind1-obj.indN),1);

            for k = 1:obj.nVolMask-obj.ind1-obj.indN
                
                temp = [-(1./(obj.SVol(k+obj.ind1)*obj.dx(k+obj.ind1)))*obj.S(obj.iarrayFace(k+obj.ind1))*droeF(:,1:3,obj.iarrayFaceI(k)),...
                         (1./(obj.SVol(k+obj.ind1)*obj.dx(k+obj.ind1)))*(obj.S(obj.iarrayFace(k+obj.ind1)+1)*droeF(:,1:3,obj.iarrayFaceI(k)+1) - obj.S(obj.iarrayFace(k+obj.ind1))*droeF(:,4:6,obj.iarrayFaceI(k)))-[zeros(1,3);dQ(1,:,k+obj.ind1);zeros(1,3)],...
                         (1./(obj.SVol(k+obj.ind1)*obj.dx(k+obj.ind1)))*obj.S(obj.iarrayFace(k+obj.ind1)+1)*droeF(:,4:6,obj.iarrayFaceI(k)+1)]';
                 J2(27*(k-1)+1:27*k)=temp(:);
            end
            dR2dp=[];
            if nargout == 3
                n=length(obj.p);
                dR2dp=zeros(obj.nVolMask,n);
                
                for k = 1:obj.nVolMask-obj.ind1-obj.indN
                    
                    ind = obj.iarraySFull(k+obj.ind1):obj.iarraySFull(k+obj.ind1)+2;
                    dR2dp(3*(k-1)+1:3*k,:) = ([-(1/(obj.SVol(k+obj.ind1)*obj.dx(k+obj.ind1)))*roeF(:,obj.iarrayFaceI(k)),...
                                               -(1/(obj.SVol(k+obj.ind1)^2*obj.dx(k+obj.ind1)))*(roeF(:,obj.iarrayFaceI(k)+1)*obj.S(obj.iarrayFace(k+obj.ind1)+1)-roeF(:,obj.iarrayFaceI(k))*obj.S(obj.iarrayFace(k+obj.ind1))),...
                                               (1/(obj.SVol(k+obj.ind1)*obj.dx(k+obj.ind1)))*roeF(:,obj.iarrayFaceI(k)+1)] ...
                                               - [zeros(1,3);dQdS(1,:,k+obj.ind1);zeros(1,3)])*obj.dSdp(ind,:);
                end
            end
        end
        
        function  [cfl] = computeCFL_GNAT(obj,U)
            
            %Reshape U from a stacked vector ([rho_1;rho_1*u_1;e_1;...]) of
            %size 3nVol x 1 to a matrix of size 3 x nVol where U(1,:) =
            %rho, U(2,:) = rho.*u, and U(3,:) = e
            U = reshape(U,3,obj.nVolSamp);
            U = U(:,obj.iarrayMask);
            n=length(obj.iarrayMask);
            %Extract primitive variables.  Take into account that the first
            %three and last three variables are PRIMITIVE VARIABLES while
            %all other are conservative variables.
            startInd = 1 + (obj.ind1&&obj.iarrayMask(1)==1);
            endInd   = n - (obj.indN&&obj.iarrayMask(end)==obj.nVolSamp);
            [~,u,~,c] = obj.conservativeToPrimitive(U(:,startInd:endInd));
            if startInd > 1
                u   = [U(2,1),u]; %Velocity
                c   = [sqrt(obj.gamma*U(3,1)/U(1,1)),c]; %Speed of sound
            end
            if endInd < n
                u   = [u,U(2,end)]; %Velocity
                c   = [c,sqrt(obj.gamma*U(3,end)/U(1,end))]; %Speed of sound
            end

            cfl = min(obj.dx./(u+c));
            
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
            
            probGNAT = quasi1dEuler([],obj);
            
            temp = probGNAT.config.inFunc;
            temp2 = probGNAT.config.form;
            
            probGNAT.config = [];
            probGNAT.config.inFunc = temp;
            probGNAT.config.form=temp2;
            probGNAT.ndim = [];
            probGNAT.Jstruct = [];
            probGNAT.dx = obj.dx(gnat.sampleNodes);
            
            %Determine sample mesh and faces surrounding volumes in mask.
            sampleMesh = union(gnat.sampleNodes,...
                       setdiff(gnat.sampleNodes,1)-1);
            sampleMesh = union(sampleMesh,...
                       setdiff(gnat.sampleNodes,obj.nVol)+1);
            sampleFaces=union(gnat.sampleNodes,gnat.sampleNodes+1);

            %Compute number of volumes in mask and sample mesh
            probGNAT.nVolSamp=length(sampleMesh);
            probGNAT.nVolMask=length(gnat.sampleNodes);
            
            %Compute positions of masked volumes and faces.
            probGNAT.mesh.node=probGNAT.mesh.node(gnat.sampleNodes);
            probGNAT.mesh.gridpt=probGNAT.mesh.gridpt(sampleFaces);
            probGNAT.mesh.conn=[];
            probGNAT.mesh.dbc=[];
                        
            %Compute logical scalars indicating whether the first and last
            %nodes are in the mask.
            probGNAT.ind1=~isempty(intersect(gnat.sampleNodes,1));
            probGNAT.indN=~isempty(intersect(gnat.sampleNodes,obj.nVol));
             
            %The index iarray is meant to index into rho, u, P, c, e stored
            %on the SAMPLE MESH.  rho(iarray) will return the ith rho values
            %needed to compute the roe fluxes F_{i+1/2}, i.e.
            %\hat{rho}_{i+1/2} = sqrt(rho(iarray)).*sqrt(rho(iarray+1)).
            temp= setdiff(gnat.sampleNodes,[1,obj.nVol]);
            iarraySM=union(temp,temp-1);
            probGNAT.iarray=zeros(1,length(iarraySM));
            for i = 1:length(iarraySM)
                probGNAT.iarray(i) = find(sampleMesh == iarraySM(i));
            end
            
            %This gives an index into an array that is stored on the sample
            %mesh (with re-indexing).  Extracts mask from sample mesh.
            probGNAT.iarrayMask=zeros(1,length(gnat.sampleNodes));
            for i = 1:length(gnat.sampleNodes)
                probGNAT.iarrayMask(i) = find(sampleMesh == gnat.sampleNodes(i));
            end
            
            %This array indexes into an array stored on the FACES
            %surrounding the sample mesh volumes.  It is used to compute the
            %residual of governing equation (R2) by indexing into S,
            %i.e. S_{i+1/2} = S(:,iarrayFace)
            probGNAT.iarrayFace=zeros(1,length(gnat.sampleNodes));
            for i = 1:length(gnat.sampleNodes)
                probGNAT.iarrayFace(i) = find(sampleFaces == gnat.sampleNodes(i));
            end

            %This array indexes into an array stored on the FACES
            %surrounding the mask volumes.  It is used to compute the
            %residual of governing equation (R2) by indexing into rhoF,
            %i.e. roeF_{i+1/2} = rhoF(:,iarrayFaceI)
            n = length(gnat.sampleNodes) - probGNAT.ind1 - probGNAT.indN;
            probGNAT.iarrayFaceI=zeros(1,n);
            probGNAT.iarrayFaceI(1)=1;
            cnt=1;
            for i = 2:n
                add1 = isempty(intersect(gnat.sampleNodes,gnat.sampleNodes(i+probGNAT.ind1)-1));
                cnt=cnt+1+add1;
                probGNAT.iarrayFaceI(i) = cnt;
            end
            
            %This array indexes into an array stored on the on the full S
            %array (S evaluated at faces and cell centers).  It is used to
            %index into dSdp to compute dR2dp.  Using this array, 
            %[dS_{i-1/2}dp,dS_{i}dp,dS_{i+1/2}dp] = dSdp(iarraySFull(i):iarraySFull+2)
            n = length(gnat.sampleNodes);
            probGNAT.iarraySFull=zeros(1,n);
            probGNAT.iarraySFull=1;
            cnt=1;
            for i = 2:n
                add1 = isempty(intersect(gnat.sampleNodes,gnat.sampleNodes(i)-1));
                cnt=cnt+2+add1;
                probGNAT.iarraySFull(i) = cnt;
            end
            
            %Re-evaluate S, SVol, xVol, dSdp on mask
            [probGNAT.S,probGNAT.SVol,probGNAT.xVol,probGNAT.dSdp] = probGNAT.computeNozzleAreaGNAT();
            probGNAT.S=(1/probGNAT.nondim.L^2)*probGNAT.S;
            probGNAT.SVol=(1/probGNAT.nondim.L^2)*probGNAT.SVol;
            probGNAT.xVol=(1/probGNAT.nondim.L)*probGNAT.xVol;
            
            %Compute B, G, ic on mask and sample mesh.
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
            %probGNAT - new "reduced" OneDBurgers object
            %--------------------------------------------------------------
            
            for kk = 1:gnat.nBases
                if isempty(gnat.sampleInd{kk,1})
                    continue;
                end
                
                probGNAT(kk) = quasi1dEuler([],obj);
                
                temp = probGNAT(kk).config.inFunc;
                temp2 = probGNAT(kk).config.form;
                
                probGNAT(kk).config = [];
                probGNAT(kk).config.inFunc = temp;
                probGNAT(kk).config.form=temp2;
                probGNAT(kk).ndim = [];
                probGNAT(kk).Jstruct = [];
                probGNAT(kk).dx = obj.dx(gnat.sampleNodes{kk});
                
                %Determine sample mesh and faces surrounding volumes in mask.
                sampleMesh = union(gnat.sampleNodes{kk},...
                    setdiff(gnat.sampleNodes{kk},1)-1);
                sampleMesh = union(sampleMesh,...
                    setdiff(gnat.sampleNodes{kk},obj.nVol)+1);
                sampleFaces=union(gnat.sampleNodes{kk},gnat.sampleNodes{kk}+1);
                
                %Compute number of volumes in mask and sample mesh
                probGNAT(kk).nVolSamp=length(sampleMesh);
                probGNAT(kk).nVolMask=length(gnat.sampleNodes{kk});
                
                %Compute positions of masked volumes and faces.
                probGNAT(kk).mesh.node=probGNAT(kk).mesh.node(gnat.sampleNodes{kk});
                probGNAT(kk).mesh.gridpt=probGNAT(kk).mesh.gridpt(sampleFaces);
                probGNAT(kk).mesh.conn=[];
                probGNAT(kk).mesh.dbc=[];
                
                %Compute logical scalars indicating whether the first and last
                %nodes are in the mask.
                probGNAT(kk).ind1=~isempty(intersect(gnat.sampleNodes{kk},1));
                probGNAT(kk).indN=~isempty(intersect(gnat.sampleNodes{kk},obj.nVol));
                
                %The index iarray is meant to index into rho, u, P, c, e stored
                %on the SAMPLE MESH.  rho(iarray) will return the ith rho values
                %needed to compute the roe fluxes F_{i+1/2}, i.e.
                %\hat{rho}_{i+1/2} = sqrt(rho(iarray)).*sqrt(rho(iarray+1)).
                temp= setdiff(gnat.sampleNodes{kk},[1,obj.nVol]);
                iarraySM=union(temp,temp-1);
                probGNAT(kk).iarray=zeros(1,length(iarraySM));
                for i = 1:length(iarraySM)
                    probGNAT(kk).iarray(i) = find(sampleMesh == iarraySM(i));
                end
                
                %This gives an index into an array that is stored on the sample
                %mesh (with re-indexing).  Extracts mask from sample mesh.
                probGNAT(kk).iarrayMask=zeros(1,length(gnat.sampleNodes{kk}));
                for i = 1:length(gnat.sampleNodes{kk})
                    probGNAT(kk).iarrayMask(i) = find(sampleMesh == gnat.sampleNodes{kk}(i));
                end
                
                %This array indexes into an array stored on the FACES
                %surrounding the sample mesh volumes.  It is used to compute the
                %residual of governing equation (R2) by indexing into S,
                %i.e. S_{i+1/2} = S(:,iarrayFace)
                probGNAT(kk).iarrayFace=zeros(1,length(gnat.sampleNodes{kk}));
                for i = 1:length(gnat.sampleNodes{kk})
                    probGNAT(kk).iarrayFace(i) = find(sampleFaces == gnat.sampleNodes{kk}(i));
                end
                
                %This array indexes into an array stored on the FACES
                %surrounding the mask volumes.  It is used to compute the
                %residual of governing equation (R2) by indexing into rhoF,
                %i.e. roeF_{i+1/2} = rhoF(:,iarrayFaceI)
                n = length(gnat.sampleNodes{kk}) - probGNAT(kk).ind1 - probGNAT(kk).indN;
                probGNAT(kk).iarrayFaceI=zeros(1,n);
                probGNAT(kk).iarrayFaceI(1)=1;
                cnt=1;
                for i = 2:n
                    add1 = isempty(intersect(gnat.sampleNodes{kk},gnat.sampleNodes{kk}(i+probGNAT(kk).ind1)-1));
                    cnt=cnt+1+add1;
                    probGNAT(kk).iarrayFaceI(i) = cnt;
                end
                
                %This array indexes into an array stored on the on the full S
                %array (S evaluated at faces and cell centers).  It is used to
                %index into dSdp to compute dR2dp.  Using this array,
                %[dS_{i-1/2}dp,dS_{i}dp,dS_{i+1/2}dp] = dSdp(iarraySFull(i):iarraySFull+2)
                n = length(gnat.sampleNodes{kk});
                probGNAT(kk).iarraySFull=zeros(1,n);
                probGNAT(kk).iarraySFull=1;
                cnt=1;
                for i = 2:n
                    add1 = isempty(intersect(gnat.sampleNodes{kk},gnat.sampleNodes{kk}(i)-1));
                    cnt=cnt+2+add1;
                    probGNAT(kk).iarraySFull(i) = cnt;
                end
                
                %Re-evaluate S, SVol, xVol, dSdp on mask
                [probGNAT(kk).S,probGNAT(kk).SVol,probGNAT(kk).xVol,probGNAT(kk).dSdp] = probGNAT(kk).computeNozzleAreaGNAT();
                probGNAT(kk).S=(1/probGNAT(kk).nondim.L)*probGNAT(kk).S;
                probGNAT(kk).SVol=(1/probGNAT(kk).nondim.L)*probGNAT(kk).SVol;
                probGNAT(kk).xVol=(1/probGNAT(kk).nondim.L)*probGNAT(kk).xVol;
                
                %Compute B, G, ic on mask and sample mesh.
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
             
        %Nozzle Functions
        function [S,SVol,xVol,dSdp] = computeNozzleArea(obj)
            
            x    = obj.mesh.gridpt';
            xVol = obj.mesh.node';
            
            switch obj.Atype
                case {1,2}
                    [S,dSIdp]    = obj.NozzleAreaSpline(x,obj.p);
                    [SVol,dSVdp] = obj.NozzleAreaSpline(xVol,obj.p);
                    dSdp=zeros(size(dSIdp,1)+size(dSVdp,1),length(obj.p)-2*(obj.Atype==1));
                case 3
                    [S,dSIdp]    = obj.NozzleArea1D(x,obj.p);
                    [SVol,dSVdp] = obj.NozzleArea1D(xVol,obj.p);
                    dSdp=zeros(size(dSIdp,1)+size(dSVdp,1),1);
                case 4
                    [S,dSIdp]    = obj.NozzleAreaAIAA(x,obj.p);
                    [SVol,dSVdp] = obj.NozzleAreaAIAA(xVol,obj.p);
                    dSdp=zeros(size(dSIdp,1)+size(dSVdp,1),length(obj.p));
            end
            
            dSdp(1:2:end,:) = dSIdp;
            dSdp(2:2:end,:) = dSVdp;
        end  
        
        function [S,SVol,xVol,dSdp] = computeNozzleAreaGNAT(obj)
            
            x    = obj.mesh.gridpt';
            xVol = obj.mesh.node';
            
            switch obj.Atype
                case 1
                    [S,dSIdp]    = obj.NozzleAreaSpline(x,obj.p);
                    [SVol,dSVdp] = obj.NozzleAreaSpline(xVol,obj.p);
                    dSdp = zeros(obj.iarraySFull(end)+2,length(obj.p)-2);
                case 2
                    [S,dSIdp]    = obj.NozzleAreaSpline(x,obj.p);
                    [SVol,dSVdp] = obj.NozzleAreaSpline(xVol,obj.p);
                    dSdp = zeros(obj.iarraySFull(end)+2,length(obj.p));
                case 3
                    [S,dSIdp]    = obj.NozzleArea1D(x,obj.p);
                    [SVol,dSVdp] = obj.NozzleArea1D(xVol,obj.p);
                    dSdp = zeros(obj.iarraySFull(end)+2,length(obj.p));
                case 4
                    [S,dSIdp]    = obj.NozzleAreaAIAA(x,obj.p);
                    [SVol,dSVdp] = obj.NozzleAreaAIAA(xVol,obj.p);
                    dSdp = zeros(obj.iarraySFull(end)+2,length(obj.p));
            end
            
            dSdp(obj.iarraySFull,:) = dSIdp(obj.iarrayFace,:);
            dSdp(obj.iarraySFull+1,:) = dSVdp;
            dSdp(obj.iarraySFull+2,:) = dSIdp(obj.iarrayFace+1,:);
        end

        function  [A,dAdp,dAdx] = NozzleArea1D(obj,xi,p)
            
            A = 0.6*(xi-1).^2 + p;
            dAdp = ones(size(A'));
            dAdx = 1.2*(xi-1);
            
        end
        
        function  [A,dAdp,dAdx] = NozzleAreaSpline(obj,xi,p)
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
            A = zeros(1,n);
            dAdx = zeros(1,n);
            dAdp = zeros(n,obj.nSpPts);
            for i = 1:n
                ind = find(obj.SpPts < xi(i),1,'last');
                if isempty(ind)
                    ind=find(obj.SpPts >= xi(i),1,'first');
                end
                
                %A(i) = a(ind)*(xi(i)-obj.SpPts(ind))^3 + b(ind)*(xi(i)-obj.SpPts(ind))^2 + c(ind)*(xi(i)-obj.SpPts(ind)) + d(ind);
                %if i == 1 || i == n
                %  continue;
                %end
                %Want the sensitivity of A with respect to the inlet and
                %outlet area to be zero because I want these fixed
                A(i) = a(ind)*(xi(i)-obj.SpPts(ind))^3 + b(ind)*(xi(i)-obj.SpPts(ind))^2 + c(ind)*(xi(i)-obj.SpPts(ind)) + d(ind);
                dAdx(i) = 3*a(ind)*(xi(i)-obj.SpPts(ind))^2 + 2*b(ind)*(xi(i)-obj.SpPts(ind)) + c(ind);
                dAdp(i,:) = dadp(ind,:)*(xi(i)-obj.SpPts(ind))^3 + dbdp(ind,:)*(xi(i)-obj.SpPts(ind))^2 + dcdp(ind,:)*(xi(i)-obj.SpPts(ind)) + dddp(ind,:);
            end
            
            if obj.Atype == 1
                dAdp = dAdp(:,2:end-1);
            end
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
            
            A=zeros(1,n);
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
        
        function  [f,df] = NozzleObj_MatchArea(obj,xi,p,Ades)
            [A,dAdp] = obj.NozzleAreaSpline(xi,p);
            f = 0.5*norm(A - Ades,2)^2;
            df = dAdp'*(A-Ades)';
        end
        
        function  [dAdx,ceq,d2Adxdp,dceq] = NozzleConstr_DerivsInOut(obj,xi,p)
            ceq = []; dceq = [];
            
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

            dadp = repmat(1./(3*Dx),1,obj.nSpPts).*(dsig(2:end,:) - dsig(1:end-1,:));
            dbdp = dsig(1:end-1,:);
     
            dcdp = zeros(obj.nSpPts-1,obj.nSpPts);
            dcdp(:,1:end-1) = -diag(1./Dx);
            dcdp(:,2:end) = dcdp(:,2:end) + diag(1./Dx);
            dcdp = dcdp + repmat(-Dx/3,1,obj.nSpPts).*(dsig(2:end,:) + 2*dsig(1:end-1,:));

            n = length(xi);
            dAdx = zeros(n,1);
            d2Adxdp = zeros(n,obj.nSpPts);
            for i = [1,2,n-1,n]
                ind = find(obj.SpPts < xi(i),1,'last');
                if isempty(ind)
                    ind=find(obj.SpPts >= xi(i),1,'first');
                end
                %A(i) = a(ind)*(xi(i)-obj.SpPts(ind))^3 + b(ind)*(xi(i)-obj.SpPts(ind))^2 + c(ind)*(xi(i)-obj.SpPts(ind)) + d(ind);
                %if i == 1 || i == n
                %  continue;
                %end
                %Want the sensitivity of A with respect to the inlet and
                %outlet area to be zero because I want these fixed
                dAdx(i) = 3*a(ind)*(xi(i)-obj.SpPts(ind))^2 + 2*b(ind)*(xi(i)-obj.SpPts(ind)) + c(ind);
                d2Adxdp(i,:) = 3*dadp(ind,:)*(xi(i)-obj.SpPts(ind))^2 + 2*dbdp(ind,:)*(xi(i)-obj.SpPts(ind)) + dcdp(ind,:);
            end
            if obj.Atype == 1
                d2Adxdp = d2Adxdp(:,2:end-1);
            end
            
            dAdx = [dAdx(1:2);-dAdx(end-1:end)];
            d2Adxdp = [d2Adxdp(1:2,:)',-d2Adxdp(end-1:end,:)'];
        end
        
        function  [d2Adx2,ceq,d3Adx2dp,dceq] = NozzleConstr_ConvexSpline(obj,xi,p)
            ceq = []; dceq = [];
            
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
            
            dadp = repmat(1./(3*Dx),1,obj.nSpPts).*(dsig(2:end,:) - dsig(1:end-1,:));
            dbdp = dsig(1:end-1,:);

            n = length(xi);
            d2Adx2 = zeros(n,1);
            d3Adx2dp = zeros(n,obj.nSpPts);
            for i = 1:n
                ind = find(obj.SpPts < xi(i),1,'last');
                if isempty(ind)
                    ind=find(obj.SpPts >= xi(i),1,'first');
                end

                %Want the sensitivity of A with respect to the inlet and
                %outlet area to be zero because I want these fixed
                d2Adx2(i)    = 6*a(ind)*(xi(i)-obj.SpPts(ind)) + 2*b(ind);
                d3Adx2dp(i,:)= 6*dadp(ind,:)*(xi(i)-obj.SpPts(ind)) + 2*dbdp(ind,:);
            end
            if obj.Atype == 1
                d3Adx2dp = d3Adx2dp(:,2:end-1);
            end
            
            d2Adx2 = -d2Adx2;
            d3Adx2dp = -d3Adx2dp';
        end

        function  [c,ceq,dc,dceq] = NozzleConstr_AreaBounds(obj,p,Al,Au)
            ceq=[]; dceq=[];
            
            if isempty(Al)
                [c,dc] = obj.NozzleAreaSpline(obj.mesh.gridpt,p);
                c = Al-c;  dc = -dc'; %Lower bound
            elseif isempty(Au)
                [c,dc] = obj.NozzleAreaSpline(obj.mesh.gridpt,p);
                c = c-Au; %Upper bound
                dc = dc';
            else
                [A,dA] = obj.NozzleAreaSpline(obj.mesh.gridpt,p);
                c = [Al-A';A'-Au];  dc = [-dA;dA]'; %upper and lower bounds
            end
            
        end
        
        function  [c,ceq,dc,dceq] = NozzleConstr_Splines(obj,p,xi,Al,Au,flag,model)   
            n=length(xi);
            ceq=[]; dceq = [];
            [c1,~,dc1,~] = obj.NozzleConstr_AreaBounds(p,Al,Au);
            [c2,~,dc2,~] = obj.NozzleConstr_ConvexSpline(xi,p);
            [c3,~,dc3,~] = obj.NozzleConstr_DerivsInOut(xi,p);
            [A,dAdp] = obj.NozzleAreaSpline(xi([1,ceil(n/2),n]),p);
            c4=[A(2)-A(1);A(2)-A(3)]; dc4=[dAdp(2,:)'-dAdp(1,:)',dAdp(2,:)'-dAdp(3,:)'];
                        
            c5=[]; dc5=[];
            if nargin >= 7 && ~isempty(model)
                if norm(p - model.curr_param) > 0
                    model.prob.updateParameters(p);
                    
                    R(1)=norm(model.prob.ResJac(model.svdUpdateData.wRef));
                    R(2)=inf;%norm(model.prob.ResJac(model.sv(:,end)));
                    for i = 1:size(model.fom_samples,2)
                        R(2+i)=norm(model.prob.ResJac(model.fom_samples(:,i)));
                    end
                    [~,min_ind] = min(R);
                    
                    if (min_ind==2), model.modifyICandTranslateBasisRk1(model.sv(:,end)); end;
                    if (min_ind>2) , model.modifyICandTranslateBasisRk1(model.fom_samples(:,min_ind-2)); end;
                    
                    model.executeModel();
                end
                
                if model.killflag
                    c5=1e6; dc5=zeros(length(p),1);
                    %fprintf('last constraint = %f \n',nan);
                else
                    w=model.sv(:,end);
                    [Res,~,dRdp] = model.prob.ResSens(w,p);
                    
                    c5 = 0.5*(Res'*Res)-0.5*200*200;
                    dc5=dRdp'*Res;
                end
                
                fprintf('last constraint = %f \n',c5);
            end
            
            if nargin < 6 || isempty(flag)
                c=[c1;c2;c3;c4;c5];
                dc=[dc1,dc2,dc3,dc4,dc5];
                return;
            elseif strcmpi(flag,'sizes')
                %for ipopt, VARIABLE NAMES DONT MAKE SENSE
                c = length(c5) + length(c4) + length(c3) + length(xi) + length(Al) + length(Au); %Number of inequality constraints
                ceq = 0; %Number of equality constraints
            elseif strcmpi(flag,'constraint')
                c = [c1;c2;c3;c4;c5]; %Nonlinear constraints [inequality; equality]
            elseif strcmpi(flag,'jacobian')
                c = sparse([dc1,dc2,dc3,dc4,dc5]');
            elseif strcmpi(flag,'jacobianstructure')
                c = sparse(ones(length(c5) + length(c4) + length(c3) + length(xi) + length(Al) + length(Au),obj.nSpPts));
            end
        end
        
        function  [c,ceq,dc,dceq] = NozzleNLConstraintsAIAA(obj,p,flag)
            if nargin < 3 || isempty(flag)
                flag='both';
            end
            
            if strcmpi(flag,'jacobianstructure')
                c=sparse(ones(4,10));
                return;
            end
            
            ceq = [];
            dceq = [];
            
            Ad = p(1); Ai = p(2);    Ae = p(3);
            k  = p(5); alpha = p(6); xstar = p(7);
            xb = p(8); xf = p(9);
            
            L = obj.mesh.gridpt(end);
            
            if strcmpi(flag,'both')
                c = [Ad + k*(xb-xstar)^2 - Ai; ...
                    Ad + k*(xf-xstar)^2 - Ae; ...
                    Ai - Ad - k*xstar^2     ;  ...
                    Ae - Ad - k*(L-xstar)^2];
                dc = [1,-1,0,0,(xb-xstar)^alpha,0,-k*alpha*(xb-xstar)^(alpha-1),k*alpha*(xb-xstar)^(alpha-1),0,0;...
                    1,0,-1,0,(xf-xstar)^alpha,0,-k*alpha*(xf-xstar)^(alpha-1),0,k*alpha*(xf-xstar)^(alpha-1),0;...
                    -1,1,0,0,-xstar^2,0,-2*k*xstar,0,0,0;...
                    -1,0,1,0,-(L-xstar)^2,0,2*k*(L-xstar),0,0,0];
            elseif strcmpi(flag,'constraint')
                c = [Ad + k*(xb-xstar)^2 - Ai; ...
                    Ad + k*(xf-xstar)^2 - Ae; ...
                    Ai - Ad - k*xstar^2     ;  ...
                    Ae - Ad - k*(L-xstar)^2];
                dc=[];
            elseif strcmpi(flag,'jacobian')
                c = [1,-1,0,0,(xb-xstar)^alpha,0,-k*alpha*(xb-xstar)^(alpha-1),k*alpha*(xb-xstar)^(alpha-1),0,0;...
                     1,0,-1,0,(xf-xstar)^alpha,0,-k*alpha*(xf-xstar)^(alpha-1),0,k*alpha*(xf-xstar)^(alpha-1),0;...
                    -1,1,0,0,-xstar^2,0,-2*k*xstar,0,0,0;...
                    -1,0,1,0,-(L-xstar)^2,0,2*k*(L-xstar),0,0,0];
                dc=[];
            end
        end
        
        function  [A,b] = NozzleLinConstraintsAIAA(obj)
            A=zeros(3,10);
            b=zeros(3,1);
            
            A(1,1)=1;  A(1,2)=-1;
            A(2,1)=1;  A(2,3)=-1;
            A(3,1)=-1; A(3,4)=1;
      %      A(4,7)=-1; A(4,8)=1; A(4,10)=0.5;
      %      A(5,7)=1;  A(5,9)=-1; A(5,10)=0.5;
        end
        
        function  [p] = findSplineFromPts(obj,N,xi,Alow,Aup)
            %nSp - number of spline points
            %N - number of points to select on screen


            newObj = quasi1dEuler([],obj);
            newObj.nSpPts = N+2;
            
            [x,y] = ginput(N);
            [X,I] = sort(x);
            
            newObj.SpPts = [obj.mesh.gridpt(1);X;obj.mesh.gridpt(end)];
            newObj.p     = interp1(X,y(I),newObj.SpPts,'linear','extrap');
            
            Ades = newObj.NozzleAreaSpline(xi,newObj.p);
            
            constr = @(p) obj.NozzleConstr_ConvexSpline(xi,p);
            
            options = optimset('Algorithm','sqp',...
                    'Display','iter-detailed',...
                    'GradObj','on',...
                    'GradConstr','on',...
                    'Hessian','bfgs',...
                    'TolFun',1e-6,...
                    'TolCon',1e-12,...
                    'TolX',1e-8);
                
            FUN = @(p) obj.NozzleObj_MatchArea(xi,p,Ades);
            pGuess = interp1(X,y(I),obj.SpPts,'linear','extrap');
            p=fmincon(FUN,pGuess,[],[],[],[],Alow,Aup,constr,options);
        end
        
        %Postprocessing Functions
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
        end %Not Done
        
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
            
        end %Not Done
        
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
        end %Not Done
        
        %%%%%%%%%%%%%%%Optimization functions%%%%%%%%%%%%%%%%%
        function  [R,dRdw,dRdp,Cp,Cpp,dCp,dCpp] = ResSens(obj,U,~)
            
            %Reshape U from a stacked vector ([rho_1;rho_1*u_1;e_1;...]) of
            %size 3nVol x 1 to a matrix of size 3 x nVol where U(1,:) =
            %rho, U(2,:) = rho.*u, and U(3,:) = e.  Except first three and
            %last three values are primitive variables while others are
            %conservative.
            U = reshape(U,3,obj.nVol);
            %Extract primitive variables.  Take into account that the first
            %three and last three variables are PRIMITIVE VARIABLES while
            %all other are conservative variables.
            [rho,u,P,c] = obj.conservativeToPrimitive(U(:,2:end-1));
            rho = [U(1,1),rho,U(1,end)]; %Density
            u   = [U(2,1),u,U(2,end)]; %Velocity
            P   = [U(3,1),P,U(3,end)]; %Pressure
            c   = [sqrt(obj.gamma*P(1)/rho(1)),c,sqrt(obj.gamma*P(end)/rho(end))]; %Speed of sound
            e   = [P(1)/(obj.gamma-1)+rho(1)*u(1)^2/2,U(3,2:end-1),P(end)/(obj.gamma-1)+rho(end)*u(end)^2/2]; %Energy
            %Derivative of sound speed w.r.t. primitive variables.
            %dc_prim = [dcdrho_1, ..., dcdrho_nVol; ...
            %           dcdu_1  , ..., dcdu_nVol; ...
            %           dcdP_1  , ..., dcdP_nVol]
            dc_prim = [-0.5*obj.gamma*P./(c.*rho.*rho);zeros(1,obj.nVol);0.5*obj.gamma./(c.*rho)]';
            %Derivative of sound speed w.r.t. conservative variables.
            %dc_cons = [dcdrho_1 , ..., dcdrho_nVol; ...
            %           dcdrhou_1, ..., dcdrhou_nVol; ...
            %           dcde_1  , ..., dcde_nVol]
            dc_cons = [0.5*obj.gamma./(c.*rho).*(0.5*(obj.gamma-1)*u.*u - P./rho);...
                      -0.5*obj.gamma*(obj.gamma-1)*u./(rho.*c);...
                       0.5*obj.gamma*(obj.gamma-1)./(rho.*c)]';
            
            %Get residual and Jacobian contributions from inlet, governing
            %equations (interior), and outlet.  Also compute matrices
            %multiplying time derivatives.
            [R1,dR1dw,Cpp,dCpp,dR1dp]=obj.fullImplicitInletBC(rho(1:2),u(1:2),P(1:2),c(1:2),dc_prim(1:2,:));
            [R2,dR2dw,dR2dp]=obj.governEqn(rho,u,P,c,e,dc_cons); %R2 is a 3 x nVol-2
            [R3,dR3dw,Cp,dCp,dR3dp]=obj.fullImplicitOutletBC(rho(end-1:end),u(end-1:end),P(end-1:end),c(end-1:end),dc_prim(end-1:end,:));
            
            %Form full residual and Jacobian.
            R=[R1;R2(:);R3];
            dRdw=[dR1dw,zeros(3,3*(obj.nVol-2)); dR2dw; zeros(3,3*(obj.nVol-2)), dR3dw];
            dRdp=[dR1dp;dR2dp;dR3dp];
        end
        
        function  [R,dRdw,dRdp] = ResSensGNAT(obj,U,~)
            
            %Reshape U from a stacked vector ([rho_1;rho_1*u_1;e_1;...]) of
            %size 3nVol x 1 to a matrix of size 3 x nVol where U(1,:) =
            %rho, U(2,:) = rho.*u, and U(3,:) = e.  Except first three and
            %last three values are primitive variables while others are
            %conservative.
            U = reshape(U,3,obj.nVolSamp);
            %Extract primitive variables.  Take into account that the first
            %three and last three variables are PRIMITIVE VARIABLES while
            %all other are conservative variables.
            startInd = 1 + obj.ind1;
            endInd   = obj.nVolSamp - obj.indN;
            [rho,u,P,c] = obj.conservativeToPrimitive(U(:,startInd:endInd));
            e=U(3,startInd:endInd);
            if obj.ind1
                rho = [U(1,1),rho]; %Density
                u   = [U(2,1),u]; %Velocity
                P   = [U(3,1),P]; %Pressure
                c   = [sqrt(obj.gamma*P(1)/rho(1)),c]; %Speed of sound
                e   = [P(1)/(obj.gamma-1)+rho(1)*u(1)^2/2,e]; %Energy
            end
            if obj.indN
                rho = [rho,U(1,end)]; %Density
                u   = [u,U(2,end)]; %Velocity
                P   = [P,U(3,end)]; %Pressure
                c   = [c,sqrt(obj.gamma*P(end)/rho(end))]; %Speed of sound
                e   = [e,P(end)/(obj.gamma-1)+rho(end)*u(end)^2/2]; %Energy
            end
            %Derivative of sound speed w.r.t. primitive variables.
            %dc_prim = [dcdrho_1, ..., dcdrho_nVol; ...
            %           dcdu_1  , ..., dcdu_nVol; ...
            %           dcdP_1  , ..., dcdP_nVol]
            dc_prim = [-0.5*obj.gamma*P./(c.*rho.*rho);zeros(1,obj.nVolSamp);0.5*obj.gamma./(c.*rho)]';
            %Derivative of sound speed w.r.t. conservative variables.
            %dc_cons = [dcdrho_1 , ..., dcdrho_nVol; ...
            %           dcdrhou_1, ..., dcdrhou_nVol; ...
            %           dcde_1  , ..., dcde_nVol]
            dc_cons = [0.5*obj.gamma./(c.*rho).*(0.5*(obj.gamma-1)*u.*u - P./rho);...
                      -0.5*obj.gamma*(obj.gamma-1)*u./(rho.*c);...
                       0.5*obj.gamma*(obj.gamma-1)./(rho.*c)]';
            
            %Get residual and Jacobian contributions from inlet, governing
            %equations (interior), and outlet.
            R1=[]; dR1dw=[]; dR1dp=[]; R3=[]; dR3dw=[]; dR3dp=[];
            if obj.ind1
                [R1,dR1dw,~,~,dR1dp]=obj.fullImplicitInletBC(rho(1:2),u(1:2),P(1:2),c(1:2),dc_prim(1:2,:));
            end
            [R2,dR2dw,dR2dp]=obj.governEqnGNAT(rho,u,P,c,e,dc_cons); %R2 is a 3 x nVol-2
            if obj.indN
                [R3,dR3dw,~,~,dR3dp]=obj.fullImplicitOutletBC(rho(end-1:end),u(end-1:end),P(end-1:end),c(end-1:end),dc_prim(end-1:end,:));
            end

            %Form residual and Jacobian.
            R=[R1;R2(:);R3];
            dR1dw=dR1dw'; dR3dw=dR3dw';
            dRdw=[dR1dw(:); dR2dw(:); dR3dw(:)];
            dRdp=[dR1dp;dR2dp;dR3dp];
        end
        
        function  [f,dfdw,dfdp] = OptSampleROM(obj,w,p)
            
            %Gamma=50./(ColumnwiseNorm(obj.pbar,2).^2);
            Gamma=0;
            
            obj.updateParameters(p);
            [Res,dRdw,dRdp] = obj.ResSens(w,[]);
            
            f = -0.5*(Res'*Res);% - sum(Gamma.*ColumnwiseNorm(bsxfun(@minus,obj.pbar,p),2).^2);
            dfdw = -dRdw'*Res;
            %The + sign of the second term comes from the fact that the
            %term inside the bsxfun is p_i - p (acutal derivative needs p - p_i)
            dfdp = -dRdp'*Res;% + 2*sum(bsxfun(@times,bsxfun(@minus,obj.pbar,p),Gamma),2);
            
            f(real(f)~=f)=nan;
            dfdw(real(dfdw)~=dfdw)=nan;
            dfdp(real(dfdp)~=dfdp)=nan;
            
            if ~isreal(f) || ~isreal(dfdw) || ~isreal(dfdp) || (sum(isnan(f)) + sum(isnan(dfdw)) + sum(isnan(dfdp)) > 0)
                f = nan;
            end
            
%             f = -0.5*(Res'*Res);
%             dfdw = -dRdw'*Res;
%             dfdp = -dRdp'*Res;
            
            fprintf('Objective function = %f\n',f);
            fprintf('Norm of Residual = %f\n',norm(Res));
            respart=100*abs((0.5*(Res'*Res))/f);
            fprintf('Residual part of objective = %f%% \n',respart);
            fprintf('Other part of objective = %f%%\n',100-respart);
        end

        function  [] = setPBAR(obj,newp)
            obj.pbar = [obj.pbar,newp];
        end
        
        function  [f,dfdw,dfdp] = OptSampleGNAT(obj,w,p)
            
            [RHat,dRdwHat,dRdpHat] = obj.ResSensGNAT(w,[]);
            obj.Jhat(obj.reconstJhatInd) = dRdwHat;
            
            f = -0.5*(RHat'*RHat);
            dfdw = -obj.Jhat'*RHat;
            dfdp = -dRdpHat'*RHat;         
        end
        
        function  [f,dfdw,dfdp] = SimOptSampleROM(obj,w,p)
            [Res,dRdw,dRdp] = obj.ResSens(w,[]);
            [fo,dfodw,dfodp] = obj.ParameterEstimation(w,p);

            f    = -0.5*obj.Gamma*(Res'*Res) + fo;
            dfdw = -obj.Gamma*dRdw'*Res + dfodw;
            dfdp = -obj.Gamma*dRdp'*Res + dfodp;
            
            if ~isreal(f) || ~isreal(dfdw) || ~isreal(dfdp) || (sum(isnan(f)) + sum(isnan(dfdw)) + sum(isnan(dfdp)) > 0)
                f = nan;
            end
            
            fprintf('Original Objective function = %f\n',fo);
            fprintf('Norm of Residual = %f\n',norm(Res));
            fprintf('Gamma*0.5*norm(Residual,2)^2 = %f\n',obj.Gamma*0.5*norm(Res)^2);
            respart=100*abs((0.5*obj.Gamma*(Res'*Res))/f);
            fprintf('Residual part of objective = %f%% \n',respart);
            fprintf('Original obj. part of objective = %f%%\n',100-respart);
        end

        function  [f,dfdw,dfdp] = SimOptSampleGNAT(obj,w,p)
            
            [RHat,dRdwHat,dRdpHat] = obj.ResSensGNAT(w,[]);
            obj.Jhat(obj.reconstJhatInd) = dRdwHat;
            [fo,dfodw,dfodp] = obj.ParameterEstimationROM(w,p);
            
            f = -0.5*obj.Gamma*(RHat'*RHat) + fo;
            dfdw = -obj.Gamma*obj.Jhat'*RHat + dfodw;
            dfdp = -obj.Gamma*dRdpHat'*RHat + dfodp;         
        end
        
        function  [f,dfdw,dfdp] = ParameterEstimation(obj,w,p)
            f = 0.5*sum((w-obj.targetSoln).^2);
            dfdp = zeros(size(p,1),1);
            dfdw = w - obj.targetSoln;
            fprintf('Objective Function Value = %f\n',f);
        end %Not Done
        
        function  [f,dfdw_hat,dfdp] = ParameterEstimationGNAT(obj,w_hat,p)
            f = 0.5*sum((w_hat-obj.targetSoln).^2);
            dfdp = zeros(size(p,1),1);
            dfdw_hat = w_hat - obj.targetSoln;
        end %Not Done
        
        function  [] = setTargetSolution(obj,wstar)
            obj.targetSoln = wstar;
        end %Not Done
        
        function  [] = updateParameters(obj,p)
            
            obj.p=p;
            
            [obj.S,obj.SVol,obj.xVol,obj.dSdp] = obj.computeNozzleArea();
            obj.S=(1/obj.nondim.L^2)*obj.S;
            obj.SVol=(1/obj.nondim.L^2)*obj.SVol;
            obj.xVol=(1/obj.nondim.L)*obj.xVol;
            obj.dSdp=(1/obj.nondim.L)*obj.dSdp;

%             %Update the parameters
%             if obj.Atype == 1
%                 %Nozzle Spline, Fixed Ends
%                 [obj.A, obj.dAdp] = obj.NozzleAreaSpline(0.5*(obj.mesh.node(2:end,1) + obj.mesh.node(1:end-1,1)),[obj.A_Bnd(1);p(:);obj.A_Bnd(2)]);
%             elseif obj.Atype == 2
%                 %Nozzle Spline, Adjustable Ends
%                 [obj.A, obj.dAdp] = obj.NozzleAreaSpline(0.5*(obj.mesh.node(2:end,1) + obj.mesh.node(1:end-1,1)),p(:));
%             else
%                 [obj.A, obj.dAdp] = obj.NozzleArea1D(0.5*(obj.mesh.node(2:end,1) + obj.mesh.node(1:end-1,1)),p);
%             end
        end %Not Done
        
        function  [] = updateParametersGNAT(obj,p)
            
            obj.p=p;
            
            [obj.S,obj.SVol,obj.xVol,obj.dSdp] = obj.computeNozzleAreaGNAT();
            obj.S=(1/obj.nondim.L^2)*obj.S;
            obj.SVol=(1/obj.nondim.L^2)*obj.SVol;
            obj.xVol=(1/obj.nondim.L)*obj.xVol;
            obj.dSdp=(1/obj.nondim.L)*obj.dSdp;

        end

        %Debugging Functions
        function  [rho,u,P,c,e,dc_prim,dc_cons] = getVariables(obj,U,flag)
            
            %Reshape U from a stacked vector ([rho_1;rho_1*u_1;e_1;...]) of
            %size 3nVol x 1 to a matrix of size 3 x nVol where U(1,:) =
            %rho, U(2,:) = rho.*u, and U(3,:) = e
            U = reshape(U,3,obj.nVol);
            %Extract primitive variables
            [rho,u,P,c] = obj.conservativeToPrimitive(U(:,2:end-1));
            rho = [U(1,1),rho,U(1,end)];
            u   = [U(2,1),u,U(2,end)];
            P   = [U(3,1),P,U(3,end)];
            c   = [sqrt(obj.gamma*P(1)/rho(1)),c,sqrt(obj.gamma*P(end)/rho(end))]; 
            e   = [P(1)/(obj.gamma-1)+rho(1)*u(1)^2/2,U(3,2:end-1),P(end)/(obj.gamma-1)+rho(end)*u(end)^2/2];
            dc_prim = [-0.5*obj.gamma*P./(c.*rho.*rho);zeros(1,obj.nVol);0.5*obj.gamma./(c.*rho)]';
            dc_cons = [0.5*obj.gamma./(c.*rho).*(0.5*(obj.gamma-1)*u.*u - P./rho);...
                      -0.5*obj.gamma*(obj.gamma-1)*u./(rho.*c);...
                       0.5*obj.gamma*(obj.gamma-1)./(rho.*c)]';

            if (nargin > 2) && ~isempty(flag) && strcmpi(flag,'dimensional')
                rho=rho*obj.nondim.rho;
                u=u*obj.nondim.u;
                P=P*obj.nondim.p;
                c=c*obj.nondim.u;
                e=e*obj.nondim.p;
                dc_prim=dc_prim*obj.nondim.u;
                dc_cons=dc_cons*obj.nondim.u;
            end
            
        end
        
        function  [drhoH,drhoHfd,duH,duHfd,dcH,dcHfd] = checkRoeTermsWithFD(obj,U)
                        
            eps=1e-4;
            
            [rho,u,P,c,e,dc_prim,dc_cons] = obj.getVariables(U);
            [rhoH,uH,cH,drhoH,duH,dcH] = obj.roeTerms(rho,u,P,e);
            
            drhoHfd=zeros(obj.nVol-1,obj.nVol,3);
            duHfd=zeros(obj.nVol-1,obj.nVol,3);
            dcHfd=zeros(obj.nVol-1,obj.nVol,3);
            
            ind=repmat(1:obj.nVol,3,1); ind=ind(:);
            for i = 1:3*obj.nVol
                ei=zeros(3*obj.nVol,1);
                ei(i)=1;
                
                [rhop,up,Pp,cp,ep,dc_primp,dc_consp] = obj.getVariables(U+eps*ei);
                [rhoHp,uHp,cHp] = obj.roeTerms(rhop,up,Pp,ep);
                
                [rhom,um,Pm,cm,em,dc_primm,dc_consm] = obj.getVariables(U-eps*ei);
                [rhoHm,uHm,cHm] = obj.roeTerms(rhom,um,Pm,em);
                
                k=mod(i-1,3)+1;
                drhoHfd(:,ind(i),k) = (1./(2*eps))*(rhoHp - rhoHm);
                duHfd(:,ind(i),k) = (1./(2*eps))*(uHp - uHm);
                dcHfd(:,ind(i),k) = (1./(2*eps))*(cHp - cHm);
            end
            tmp=drhoHfd;
            drhoHfd=zeros(obj.nVol-1,2,3);
            drhoHfd(:,1,1)=diag(tmp(:,:,1));
            drhoHfd(:,1,2)=diag(tmp(:,:,2));
            drhoHfd(:,1,3)=diag(tmp(:,:,3));
            drhoHfd(:,2,1)=diag(tmp(:,:,1),1);
            drhoHfd(:,2,2)=diag(tmp(:,:,2),1);
            drhoHfd(:,2,3)=diag(tmp(:,:,3),1);
            
            tmp=duHfd;
            duHfd=zeros(obj.nVol-1,2,3);
            duHfd(:,1,1)=diag(tmp(:,:,1));
            duHfd(:,1,2)=diag(tmp(:,:,2));
            duHfd(:,1,3)=diag(tmp(:,:,3));
            duHfd(:,2,1)=diag(tmp(:,:,1),1);
            duHfd(:,2,2)=diag(tmp(:,:,2),1);
            duHfd(:,2,3)=diag(tmp(:,:,3),1);
            
            tmp=dcHfd;
            dcHfd=zeros(obj.nVol-1,2,3);
            dcHfd(:,1,1)=diag(tmp(:,:,1));
            dcHfd(:,1,2)=diag(tmp(:,:,2));
            dcHfd(:,1,3)=diag(tmp(:,:,3));
            dcHfd(:,2,1)=diag(tmp(:,:,1),1);
            dcHfd(:,2,2)=diag(tmp(:,:,2),1);
            dcHfd(:,2,3)=diag(tmp(:,:,3),1);
        end
        
        function  [droeF,droeFfd,droeFdiff] = checkRoeFluxWithFD(obj,U)
                        
            eps=1e-6;
            
            [rho,u,P,c,e,dc_prim,dc_cons] = obj.getVariables(U);
            [roeF,droeF] = obj.roeFlux(rho,u,P,c,e,dc_cons);
            dUdV=[1,0,0;u(1),rho(1),0;0.5*u(1)*u(1),rho(1)*u(1),1/(obj.gamma-1)];
            droeF(:,1:3,1)=droeF(:,1:3,1)*dUdV;
            dUdV=[1,0,0;u(end),rho(end),0;0.5*u(end)*u(end),rho(end)*u(end),1/(obj.gamma-1)];
            droeF(:,4:6,end)=droeF(:,4:6,end)*dUdV;
            
            
            droeFfd = zeros(3,obj.nVol-1,3*obj.nVol);
            
            ind=repmat(1:obj.nVol,3,1); ind=ind(:);
            for i = 1:3*obj.nVol
                ei=zeros(3*obj.nVol,1);
                ei(i)=1;
                
                [rhop,up,Pp,cp,ep,dc_primp,dc_consp] = obj.getVariables(U+eps*ei);
                [roeFp] = obj.roeFlux(rhop,up,Pp,cp,ep,dc_consp);
                
                [rhom,um,Pm,cm,em,dc_primm,dc_consm] = obj.getVariables(U-eps*ei);
                [roeFm] = obj.roeFlux(rhom,um,Pm,cm,em,dc_consm);
                
                %k=mod(i-1,3)+1;
                n=ind(i);
                
                droeFfd(:,:,i) = (1./(2*eps))*(roeFp - roeFm);
%                 for j = 1:obj.nVol
%                     for k = 1:3
%                         drhoHfd(:,j,k) = 
%                     end
%                 end
            end
            
            droeFdiff=zeros(size(droeF,3),1);
            for i = 1:size(droeF,3)
                droeFdiff(i)=max(max(squeeze(droeFfd(:,i,3*(i-1)+1:3*(i-1)+6))-droeF(:,:,i)));
            end
        end
        
        function  [dS,dSi,dCA,dCAi,dSfd,dSifd,dCAfd,dCAifd] = checkEulerJacDerWithFD(obj,U)
            eps=1e-4;
            
            [rho,u,P,c,e,dc_prim,dc_cons] = obj.getVariables(U);
            
            dS=zeros(3,3,obj.nVol,2);
            dSi=zeros(3,3,obj.nVol,2);
            dCA=zeros(3,3,obj.nVol,2);
            dCAi=zeros(3,3,obj.nVol,2);
            
            dSfd=zeros(3,3,obj.nVol,3);
            dSifd=zeros(3,3,obj.nVol,3);
            dCAfd=zeros(3,3,obj.nVol,3);
            dCAifd=zeros(3,3,obj.nVol,3);
            
            for i = 1:obj.nVol

                
                [~,~,~,~,~,dS(:,:,i,:),dSi(:,:,i,:),dCA(:,:,i,:),dCAi(:,:,i,:)] = obj.eulerJac(rho(i),u(i),c(i));
                                              
                [Sp,Sip,CAp,CAip] = obj.eulerJac(rho(i)+eps,u(i),c(i));
                [Sm,Sim,CAm,CAim] = obj.eulerJac(rho(i)-eps,u(i),c(i));
                
                dSfd(:,:,i,1) = (1./(2*eps))*(Sp - Sm);
                dSifd(:,:,i,1) = (1./(2*eps))*(Sip - Sim);
                dCAfd(:,:,i,1) = (1./(2*eps))*(CAp - CAm);
                dCAifd(:,:,i,1) = (1./(2*eps))*(CAip - CAim);
                
                
                [Sp,Sip,CAp,CAip] = obj.eulerJac(rho(i),u(i)+eps,c(i));
                [Sm,Sim,CAm,CAim] = obj.eulerJac(rho(i),u(i)-eps,c(i));
                
                dSfd(:,:,i,2) = (1./(2*eps))*(Sp - Sm);
                dSifd(:,:,i,2) = (1./(2*eps))*(Sip - Sim);
                dCAfd(:,:,i,2) = (1./(2*eps))*(CAp - CAm);
                dCAifd(:,:,i,2) = (1./(2*eps))*(CAip - CAim);
                
                
                [Sp,Sip,CAp,CAip] = obj.eulerJac(rho(i),u(i),c(i)+eps);
                [Sm,Sim,CAm,CAim] = obj.eulerJac(rho(i),u(i),c(i)-eps);
                
                dSfd(:,:,i,3) = (1./(2*eps))*(Sp - Sm);
                dSifd(:,:,i,3) = (1./(2*eps))*(Sip - Sim);
                dCAfd(:,:,i,3) = (1./(2*eps))*(CAp - CAm);
                dCAifd(:,:,i,3) = (1./(2*eps))*(CAip - CAim);
            end
        end

        function  [A,Afd] = checkEulerJacWithFD(obj,U)
            eps=1e-4;
            
            [rho,u,P,c,e,dc_prim,dc_cons] = obj.getVariables(U);
            
            A = zeros(3,3,obj.nVol);
            Afd = zeros(3,3,obj.nVol);
            
            for i = 1:obj.nVol

                
                [S,Si,CA,CAi,Lam] = obj.eulerJac(rho(i),u(i),c(i));
                A(:,:,i)=Si*CAi*Lam*CA*S;
                
                for j = 1:3
                    ei=zeros(3,obj.nVol);
                    ei(j,:)=1;
                    
                    U=reshape(U,3,obj.nVol);
                    Up=reshape(U+eps*ei,3*obj.nVol,1);
                    Um=reshape(U-eps*ei,3*obj.nVol,1);
                    
                    [rhop,up,Pp,cp,ep] = obj.getVariables(Up);
                    [rhom,um,Pm,cm,em] = obj.getVariables(Um);
                    
                    Fp = obj.eulerNonlin(rhop,up,ep,Pp);
                    Fm = obj.eulerNonlin(rhom,um,em,Pm);
                    
                    Afd(:,j,:)=(0.5/eps)*(Fp-Fm);
                end
            end
        end
    end
    end

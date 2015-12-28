classdef quasiEuler1DImplicitPseudoTime < handle
    
    properties (SetAccess = private, GetAccess = public)
        %Input Properties
        dt;
        ndof;
        
        explicit = false;
        
        uniqueJROWind;
        jdiagHat;
        jdiag;
        computeJPHIhat;
    end
    
    methods
        %Constructor
        function  [obj] = quasiEuler1DImplicitPseudoTime(MODEL)
            %This is the constructor of the TimeIntegrate class.  It reads
            %the appropriate entries from the input file and stores them in
            %the class instance.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %FOM - FOM object
            %
            %Outputs:
            %--------
            %obj - instance of the TimeIntegrate object that was
            %      constructed
            %--------------------------------------------------------------
            
            obj.dt = MODEL.time.dt;
            %obj.dt = MODEL.time.dt/MODEL.prob.nondim.t;
            if ~(strcmpi(class(MODEL),'GNAT') || strcmpi(class(MODEL),'locGNAT'))
                obj.ndof = MODEL.prob.config.ndof;
            end
            
        end
        
        function  [] = addGNAT(obj,GNAT)
            
            obj.uniqueJROWind = GNAT.uniqueJROWind;
            obj.jdiagHat = GNAT.jdiagHat;
            obj.jdiag = GNAT.jdiag;
            obj.computeJPHIhat = @(J) GNAT.computeJPHIhat(J);
            
        end
        
        function  [] = setDT(obj,dti)
            obj.dt=dti;
        end
        
        %Entire functions
        function  [R,J,prevsv] = TimeIntNLFunc(obj,prob,U,Uprev,t)
            %This function sets up parameters to be used in Backward Euler calculation
            %dX/dt = R(X) where R is the residual of the nonlinear function
            %--------------------------------------------------------------------------
            %Inputs:
            %-------
            %modelobj     - Model object
            %t            - 1 x (nstep+1) vector containing the times at which the
            %               solution is to be computed (first entry is the time of the
            %               initial condition)
            %dt           - scalar indicating the time step size
            %Outputs:
            %--------
            %R            - the residual of the nonlinear function, adjusted for the
            %               Backward Euler scheme (wrapped in the nonlinear solver)
            %J            - the jacobian of the nonlinear function, adjusted for the
            %               Backward Euler scheme (wrapped in the nonlinear solver)
            %--------------------------------------------------------------------------
            prevsv=[];
%             %Reshape U
%             N = length(U);
%             U = reshape(U,3,N/3)';
%             % Boundary Conditions
%             [Ainlet,Binlet,dUinlet,AdjustInlet] = prob.fullyImplicitInletBC(U);
%             % at the outlet
%             [Aoutlet,Coutlet,dUoutlet,AdjustOutlet] = prob.fullyImplicitOutletBC(U);
%             % Build the Fluxes
%             F = prob.fluxes(U);
%             % Build Roe matrices
%             Jac = prob.roeMat(U);
%             % Build Roe fluxes
%             roeF = prob.roeFlux(U,Jac,F);
%             %Extract the pressure from U
%             P = prob.extractP(U);
%             % Build the RHS
%             R = obj.dt*prob.buildRHS(dUinlet,dUoutlet,roeF,P)';
%             R = R(:);
%             % Build the LHS
%             J = prob.buildLHS(U,Jac,Ainlet,Binlet,Aoutlet,Coutlet);
% %             J = prob.buildLHS(U,Jac,obj.dt*Ainlet - AdjustInlet,Binlet,obj.dt*Aoutlet-AdjustOutlet,Coutlet);
%             J(3:end,:) = obj.dt*J(3:end,:);
%             J = J - blkdiag(AdjustInlet,eye(N-6),AdjustOutlet);
% %             J([1:N-3,N],:) = obj.dt*J([1:N-3,N],:);
% %             J = J - blkdiag(AdjustOutlet,eye(N-6),AdjustInlet);
            
            N=length(U);
            %Determine the Residual and Jacobian (function handles) of the continuous
            %system of equations
            [R,J,Cp,Cpp,dCp,dCpp] = prob.ResJac(U,t);
            
            dU=U-Uprev;
            %Add contribution to the residual for the Backward Euler scheme)
            R(1:3)     = R(1:3) + Cpp*dU(1:3)/obj.dt;
            R(4:end-3) = R(4:end-3) + dU(4:end-3)/obj.dt;
            R(end-2:end)=R(end-2:end) + Cp*dU(end-2:end)/obj.dt;
            
            %Set up the appropriate jacobian function for the
            %backward euler integration scheme used (see notes for
            %details)
%             J=obj.dt*J;

            J(1:3,1:3) = J(1:3,1:3) + (Cpp + [dCpp(:,:,1)*dU(1:3),dCpp(:,:,2)*dU(1:3),dCpp(:,:,3)*dU(1:3)])/obj.dt;
            J(4:end-3,4:end-3)=J(4:end-3,4:end-3) + speye(N-6)/obj.dt;
            J(end-2:end,end-2:end) = J(end-2:end,end-2:end) + (Cp + [dCp(:,:,1)*dU(end-2:end),dCp(:,:,2)*dU(end-2:end),dCp(:,:,3)*dU(end-2:end)])/obj.dt;
        end
        
        function  [R,J,dRdp,prevsv] = TimeIntNLFuncSens(obj,prob,U,Uprev,t)
            %This function sets up parameters to be used in Backward Euler calculation
            %dX/dt = R(X) where R is the residual of the nonlinear function
            %--------------------------------------------------------------------------
            %Inputs:
            %-------
            %modelobj     - Model object
            %t            - 1 x (nstep+1) vector containing the times at which the
            %               solution is to be computed (first entry is the time of the
            %               initial condition)
            %dt           - scalar indicating the time step size
            %Outputs:
            %--------
            %R            - the residual of the nonlinear function, adjusted for the
            %               Backward Euler scheme (wrapped in the nonlinear solver)
            %J            - the jacobian of the nonlinear function, adjusted for the
            %               Backward Euler scheme (wrapped in the nonlinear solver)
            %--------------------------------------------------------------------------
            prevsv=[];

            N=length(U);
            %Determine the Residual and Jacobian (function handles) of the continuous
            %system of equations
            [R,J,dRdp,Cp,Cpp,dCp,dCpp] = prob.ResSens(U,t);
            
            dU=(1/obj.dt)*(U-Uprev);
            %Add contribution to the residual for the Backward Euler scheme)
            R(1:3)     = R(1:3) + Cpp*dU(1:3);
            R(4:end-3) = R(4:end-3) + dU(4:end-3);
            R(end-2:end)=R(end-2:end) + Cp*dU(end-2:end);
            
            %Set up the appropriate jacobian function for the
            %backward euler integration scheme used (see notes for
            %details)
%             J=obj.dt*J;
            J(1:3,1:3) = J(1:3,1:3) + (1/obj.dt)*Cpp + [dCpp(:,:,1)*dU(1:3),dCpp(:,:,2)*dU(1:3),dCpp(:,:,3)*dU(1:3)];
            J(4:end-3,4:end-3)=J(4:end-3,4:end-3) + (1/obj.dt)*speye(N-6);
            J(end-2:end,end-2:end) = J(end-2:end,end-2:end) + (1/obj.dt)*Cp + [dCp(:,:,1)*dU(end-2:end),dCp(:,:,2)*dU(end-2:end),dCp(:,:,3)*dU(end-2:end)];
      
%             %Sensitivity
%             dRdp = dRdp;
        end
        
        function  [Rhat,JVhat] = TimeIntNLFuncGNAT(obj,probGNAT,partialU,partialUprev,t)
            %This function sets up parameters to be used in Backward Euler calculation
            %dX/dt = R(X) where R is the residual of the nonlinear function
            %--------------------------------------------------------------------------
            %Inputs:
            %-------
            %modelobj     - Model object
            %t            - 1 x (nstep+1) vector containing the times at which the
            %               solution is to be computed (first entry is the time of the
            %               initial condition)
            %dt           - scalar indicating the time step size
            %Outputs:
            %--------
            %R            - the residual of the nonlinear function, adjusted for the
            %               Backward Euler scheme (wrapped in the nonlinear solver)
            %J            - the jacobian of the nonlinear function, adjusted for the
            %               Backward Euler scheme (wrapped in the nonlinear solver)
            %--------------------------------------------------------------------------
          
            [Rhat,Jhat,Cp,Cpp,dCp,dCpp] = probGNAT.ResJacGNAT(partialU,t);
            dU=partialU-partialUprev;
             
%             Rhat=obj.dt*Rhat;
%             Jhat=obj.dt*Jhat;
            
            if probGNAT.ind1
                Rhat(1:3) = Rhat(1:3) + Cpp*dU(1:3)/obj.dt;
                J1time = Cpp/obj.dt + [dCpp(:,:,1)*dU(1:3),dCpp(:,:,2)*dU(1:3),dCpp(:,:,3)*dU(1:3)]/obj.dt;
                J1time=J1time';
                Jhat([1:3,7:9,13:15]) = Jhat([1:3,7:9,13:15]) + J1time(:);
            end
            
            dUmask=dU(obj.jdiagHat);
            Rhat(1+3*probGNAT.ind1:end-3*probGNAT.indN) = Rhat(1+3*probGNAT.ind1:end-3*probGNAT.indN) + dUmask(1+3*probGNAT.ind1:end-3*probGNAT.indN)/obj.dt;
            Jhat(1+18*probGNAT.ind1:end-18*probGNAT.indN) = Jhat(1+18*probGNAT.ind1:end-18*probGNAT.indN) + obj.jdiag(1+18*probGNAT.ind1:end-18*probGNAT.indN)/obj.dt;
            if probGNAT.indN
                Rhat(end-2:end) = Rhat(end-2:end) + Cp*dU(end-2:end)/obj.dt;
                J3time = Cp/obj.dt + [dCp(:,:,1)*dU(end-2:end),dCp(:,:,2)*dU(end-2:end),dCp(:,:,3)*dU(end-2:end)]/obj.dt;
                J3time=J3time';
                Jhat([end-14:end-12,end-8:end-6,end-2:end]) = Jhat([end-14:end-12,end-8:end-6,end-2:end]) + J3time(:);
                %Jhat([end-17:end-15,end-11:end-9,end-5:end-3]) = Jhat([end-17:end-15,end-11:end-9,end-5:end-3]) + J3time(:);
            end
            JVhat = obj.computeJPHIhat(Jhat);
            %keyboard
	end
        
        function  [Rhat,JVhat] = TimeIntNLFuncGNATloc(obj,probGNAT,partialU,partialUprev,t,cLoc)
            %This function sets up parameters to be used in Backward Euler calculation
            %dX/dt = R(X) where R is the residual of the nonlinear function
            %--------------------------------------------------------------------------
            %Inputs:
            %-------
            %modelobj     - Model object
            %t            - 1 x (nstep+1) vector containing the times at which the
            %               solution is to be computed (first entry is the time of the
            %               initial condition)
            %dt           - scalar indicating the time step size
            %Outputs:
            %--------
            %R            - the residual of the nonlinear function, adjusted for the
            %               Backward Euler scheme (wrapped in the nonlinear solver)
            %J            - the jacobian of the nonlinear function, adjusted for the
            %               Backward Euler scheme (wrapped in the nonlinear solver)
            %--------------------------------------------------------------------------
            
            %This function sets up parameters to be used in Backward Euler calculation
            %dX/dt = R(X) where R is the residual of the nonlinear function
            %--------------------------------------------------------------------------
            %Inputs:
            %-------
            %modelobj     - Model object
            %t            - 1 x (nstep+1) vector containing the times at which the
            %               solution is to be computed (first entry is the time of the
            %               initial condition)
            %dt           - scalar indicating the time step size
            %Outputs:
            %--------
            %R            - the residual of the nonlinear function, adjusted for the
            %               Backward Euler scheme (wrapped in the nonlinear solver)
            %J            - the jacobian of the nonlinear function, adjusted for the
            %               Backward Euler scheme (wrapped in the nonlinear solver)
            %--------------------------------------------------------------------------
          
            [Rhat,Jhat,Cp,Cpp,dCp,dCpp] = probGNAT(cLoc).ResJacGNAT(partialU,t);
            dU=partialU-partialUprev;
            
%             Rhat=obj.dt*Rhat;
%             Jhat=obj.dt*Jhat;
            
            if probGNAT(cLoc).ind1
                Rhat(1:3) = Rhat(1:3) + Cpp*dU(1:3)/obj.dt;
                J1time = Cpp/obj.dt + [dCpp(:,:,1)*dU(1:3),dCpp(:,:,2)*dU(1:3),dCpp(:,:,3)*dU(1:3)]/obj.dt;
                J1time=J1time';
                Jhat([1:3,7:9,13:15]) = Jhat([1:3,7:9,13:15]) + J1time(:);
            end
            
            dUmask=dU(obj.jdiagHat{cLoc});
            Rhat(1+3*probGNAT(cLoc).ind1:end-3*probGNAT(cLoc).indN) = Rhat(1+3*probGNAT(cLoc).ind1:end-3*probGNAT(cLoc).indN) + dUmask(1+3*probGNAT(cLoc).ind1:end-3*probGNAT(cLoc).indN)/obj.dt;
            Jhat(1+18*probGNAT(cLoc).ind1:end-18*probGNAT(cLoc).indN) = Jhat(1+18*probGNAT(cLoc).ind1:end-18*probGNAT(cLoc).indN) + obj.jdiag{cLoc}(1+18*probGNAT(cLoc).ind1:end-18*probGNAT(cLoc).indN)/obj.dt;
            if probGNAT(cLoc).indN
                Rhat(end-2:end) = Rhat(end-2:end) + Cp*dU(end-2:end)/obj.dt;
                J3time = Cp/obj.dt + [dCp(:,:,1)*dU(end-2:end),dCp(:,:,2)*dU(end-2:end),dCp(:,:,3)*dU(end-2:end)]/obj.dt;
                J3time=J3time';
                Jhat([end-14:end-12,end-8:end-6,end-2:end]) = Jhat([end-14:end-12,end-8:end-6,end-2:end]) + J3time(:);
                %Jhat([end-17:end-15,end-11:end-9,end-5:end-3]) = Jhat([end-17:end-15,end-11:end-9,end-5:end-3]) + J3time(:);
            end
            JVhat = obj.computeJPHIhat(Jhat);
            
%             %Determine the Residual and Jacobian (function handles) of the continuous
%             %system of equations
%             [Rhat,Jhat] = probGNAT(cLoc).ResJacGNAT(partialU(obj.uniqueJROWind{cLoc}),t);
%             %Add contribution to the residual for the Backward Euler scheme)
%             Rhat = obj.dt*Rhat - (partialU(obj.jdiagHat{cLoc},1)-partialUprev(obj.jdiagHat{cLoc} ,1));
%             %Set up the appropriate jacobian function for the
%             %backward euler integration scheme used (see notes for
%             %details)
%             JVhat = obj.computeJPHIhat(obj.dt*Jhat - obj.jdiag{cLoc});
        end
        
        %Directional functions
        function  [phi,dphi] = TimeIntNLFuncDir(obj,prob,U,Uprev,t,alpha,v)
            
            [R,J] = TimeIntNLFunc(obj,prob,U + alpha*v,Uprev,t);
            
            phi = 0.5*norm(R,2)^2;
            dphi = v'*(J'*R);
        end
        
        function  [phi,dphi] = TimeIntNLFuncGNATDir(obj,probGNAT,U,Uprev,t,alpha,v)
            
            [R,J] = TimeIntNLFuncGNAT(obj,probGNAT,U + alpha*v,Uprev,t);
            
            phi = 0.5*norm(R,2)^2;
            dphi = v'*(J'*R);
        end
        
        function  [phi,dphi] = TimeIntNLFuncGNATlocDir(obj,probGNAT,U,Uprev,t,cLoc,alpha,v)
            
            [R,J] = TimeIntNLFuncGNATloc(obj,probGNAT,U + alpha*v,Uprev,t,cLoc);
            
            phi = 0.5*norm(R,2)^2;
            dphi = v'*(J'*R);
        end
    end
end

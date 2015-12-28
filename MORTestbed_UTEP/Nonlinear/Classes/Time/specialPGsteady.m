classdef specialPGsteady < handle
    
    properties (SetAccess = private, GetAccess = public)
        %Input Properties
        ndof;
        modelType;
        
        rob;
        
        uniqueJROWind;
        computeJPHIhat;
        
        dt;
    end
    
    methods
        %Constructor
        function  [obj] = specialPGsteady(MODEL)
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
            
            obj.ndof = MODEL.prob.config.ndof;
            switch class(MODEL)
                case 'FOM'
                    obj.modelType = 'FOM';
                case 'ROM'
                    obj.modelType = MODEL.Rtype;
            end
        end
        
        function  [] = addGNAT(obj,GNAT)
            
            obj.uniqueJROWind = GNAT.uniqueJROWind;
            obj.computeJPHIhat = @(J) GNAT.computeJPHIhat(J);
            
        end
        
        function  [] = setDT(obj,dti)
            obj.dt=dti;
        end
        
        %Entire functions
        function  [R,J,iprev] = TimeIntNLFunc(obj,prob,U,Up,t)
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
            
            %Determine the Residual and Jacobian (function handles) of the continuous
            %system of equations
            [r,j] = prob.ResJac(U,t);
            
            R = (1/obj.dt)*(U-Up) + j'*r;
            J = (1/obj.dt)*speye(obj.ndof) + j'*j;
            
            iprev=[];
        end
        
        function  [Rhat,JVhat] = TimeIntNLFuncGNAT(obj,probGNAT,partialU,~,t)
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
            
            %Determine the Residual and Jacobian (function handles) of the continuous
            %system of equations
%             [Rhat,Jhat] = probGNAT.ResJacGNAT(partialU(obj.uniqueJROWind),t);
%             JVhat = obj.computeJPHIhat(Jhat);
        end
        
        function  [Rhat,JVhat] = TimeIntNLFuncGNATloc(obj,probGNAT,partialU,~,t,cLoc)
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
            
            %Determine the Residual and Jacobian (function handles) of the continuous
            %system of equations
            [Rhat,Jhat] = probGNAT(cLoc).ResJacGNAT(partialU(obj.uniqueJROWind{cLoc}),t);
            JVhat = obj.computeJPHIhat(Jhat);
        end
        
        %Directional functions
        function  [phi,dphi] = TimeIntNLFuncDir(obj,prob,U,~,t,alpha,v,rom)
                        
            switch lower(obj.modelType(1))
                case 'f'
                    [R,J] = TimeIntNLFunc(obj,prob,U + alpha*v,[],t);

                    phi = 0.5*norm(R,2)^2;
                    dphi = v'*(J'*R);
                case 'g'
                    
                case 'p'
                    w = U + alpha*(rom.phi*v);
                    [R,J] = TimeIntNLFunc(obj,prob,w,[],t);
                    
                    phi = 0.5*norm(R)^2;
                    if sum(isnan(R))>0 || ~isreal(R)
                        phi=inf;
                    end
                    dphi = v'*(rom.phi'*(J'*R));
%                     
%                     dphi= v'*(Jr'*Jr + d2Rdu2_r)*Rr;
%                     w = U + alpha*(rom.phi*v);
%                     [R,J] = TimeIntNLFunc(obj,prob,w,[],t);
%                     
%                     Rr = rom.phi'*(J'*R);
%                     phi = 0.5*norm(Rr)^2;
%                     if nargout == 1, dphi=[]; return; end;
%                     
%                     Jr = J*rom.phi;
%                     %d2Rdu2_r = SecondDerivsPG(rom,R,J,[],rom.phi,1e-5,t,w,2);
%                     d2Rdu2_r=zeros(rom.nY);
%                     
%                     dphi= v'*(Jr'*Jr + d2Rdu2_r)*Rr;
            end
        end
        
        function  [phi,dphi] = TimeIntNLFuncGNATDir(obj,probGNAT,U,~,t,alpha,v)
            
            [R,J] = TimeIntNLFuncGNAT(obj,probGNAT,U + alpha*v,[],t);
            
            phi = 0.5*norm(R,2)^2;
            dphi = v'*(J'*R);
        end
        
        function  [phi,dphi] = TimeIntNLFuncGNATlocDir(obj,probGNAT,U,~,t,alpha,v)
            
            [R,J] = TimeIntNLFuncGNATloc(obj,probGNAT,U + alpha*v,[],t);
            
            phi = 0.5*norm(R,2)^2;
            dphi = v'*(J'*R);
        end
    end
end
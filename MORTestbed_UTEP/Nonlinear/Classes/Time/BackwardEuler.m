classdef BackwardEuler < handle
    
    properties (SetAccess = private, GetAccess = public)
        %Input Properties
        dt;
        ndof;
        
        explicit = false;
        
        uniqueJROWind;
        jdiagHat;
        jdiag;
        jrow;
        computeJPHIhat;
    end
    
    methods
        %Constructor
        function  [obj] = BackwardEuler(MODEL)
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
            if ~(strcmpi(class(MODEL),'GNAT') || strcmpi(class(MODEL),'locGNAT'))
                obj.ndof = MODEL.prob.config.ndof;
            end
            
        end
        
        function  [] = addGNAT(obj,GNAT)
            
            obj.uniqueJROWind = GNAT.uniqueJROWind;
            obj.jdiagHat = GNAT.jdiagHat;
            obj.jdiag = GNAT.jdiag;
            obj.jrow = GNAT.jrow;
            obj.computeJPHIhat = @(J) GNAT.computeJPHIhat(J);
            
        end
        
        function  [] = setDT(obj,dti)
            obj.dt=dti;
        end

        %Entire functions
        function  [R,J,iprevsv] = TimeIntNLFunc(obj,prob,U,Uprev,t)
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
            iprevsv=[];
            
            %Determine the Residual and Jacobian (function handles) of the continuous
            %system of equations
            [R,J] = prob.ResJac(U,t);
            %Add contribution to the residual for the Backward Euler scheme) 
            %Set up the appropriate jacobian function for the
            %backward euler integration scheme used (see notes for
            %details)
            switch prob.config.form
                case 'non_descriptor'
                    R = obj.dt*R - (U - Uprev);
                    J = obj.dt*J - speye(obj.ndof);
                otherwise
                    [RT,JT] = prob.ResJacTimeDer(U);
                    [RTprev,~] = prob.ResJacTimeDer(Uprev);
                    R = obj.dt*R - (RT - RTprev);
                    J = obj.dt*J - JT;
            end
            
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
            
            [Rhat,Jhat] = probGNAT.ResJacGNAT(partialU(obj.uniqueJROWind),t);
            switch probGNAT.config.form
                case 'non_descriptor'
                    %Add contribution to the residual for the Backward Euler scheme)
                    Rhat = obj.dt*Rhat - (partialU(obj.jdiagHat ,1)-partialUprev(obj.jdiagHat ,1));
                    %Set up the appropriate jacobian function for the
                    %backward euler integration scheme used (see notes for
                    %details)
                    JVhat = obj.computeJPHIhat(obj.dt*Jhat - obj.jdiag);
                otherwise
                    [RThat,JThat] = probGNAT.ResJacGNATTimeDer(partialU(obj.jdiagHat ,1),[],[],[],[],obj.jdiag);
                    [RThatprev,~] = probGNAT.ResJacGNATTimeDer(partialUprev(obj.jdiagHat ,1),[],[],[],[],obj.jdiag);
                    
                    
                    Rhat = obj.dt*Rhat - (RThat - RThatprev);
                    JVhat = obj.computeJPHIhat(obj.dt*Jhat - JThat); 
            end
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
            
            %Determine the Residual and Jacobian (function handles) of the continuous
            %system of equations
            [Rhat,Jhat] = probGNAT(cLoc).ResJacGNAT(partialU(obj.uniqueJROWind{cLoc}),t);
            switch probGNAT(cLoc).config.form
                case 'non_descriptor'
                    %Add contribution to the residual for the Backward Euler scheme)
                    Rhat = obj.dt*Rhat - (partialU(obj.jdiagHat{cLoc},1)-partialUprev(obj.jdiagHat{cLoc} ,1));
                    %Set up the appropriate jacobian function for the
                    %backward euler integration scheme used (see notes for
                    %details)
                    JVhat = obj.computeJPHIhat(obj.dt*Jhat - obj.jdiag{cLoc});
                otherwise
                    [RThat,JThat] = probGNAT(cLoc).ResJacGNATTimeDer(partialU(obj.jdiagHat{cLoc} ,1),[],[],[],[],obj.jdiag{cLoc});
                    [RThatprev,~] = probGNAT(cLoc).ResJacGNATTimeDer(partialUprev(obj.jdiagHat{cLoc} ,1),[],[],[],[],obj.jdiag{cLoc});
                    
                    Rhat = obj.dt*Rhat - (RThat - RThatprev);
                    JVhat = obj.computeJPHIhat(obj.dt*Jhat - JThat); 
            end
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

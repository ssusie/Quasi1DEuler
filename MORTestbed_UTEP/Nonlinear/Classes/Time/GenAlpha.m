classdef GenAlpha < handle
    
    properties (SetAccess = private, GetAccess = public)
        %Input Properties
        dt;
        ndof;
        
        beta;
        gamma;
        alpha_m;
        alpha_f;
        
        explicit=false;
        
        uniqueJROWind;
        jdiagHat;
        jdiag;
        jrow;
        computeJPHIhat;
    end
    
    methods
        %Constructor
        function  [obj] = GenAlpha(MODEL)
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
            
            if length(MODEL.time.genalpha)==1
                rho = MODEL.time.genalpha;
                
                obj.alpha_f = rho/(rho+1);
                obj.alpha_m = (2*rho-1)/(rho+1);
                obj.beta    = 0.25*(1-obj.alpha_m+obj.alpha_f)^2;
                obj.gamma   = 0.5 - obj.alpha_f + obj.alpha_m;
            elseif length(MODEL.time.genalpha)==4
                obj.beta = MODEL.time.genalpha(1);
                obj.gamma = MODEL.time.genalpha(2);
                obj.alpha_m = MODEL.time.genalpha(3);
                obj.alpha_f = MODEL.time.genalpha(4);
            end
            
            if obj.beta == 0, obj.explicit = true; end
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
        
        function  [] = setGenAlpha(obj,genalpha)
            if length(genalpha)==1
                rho = genalpha;
                
                obj.alpha_f = rho/(rho+1);
                obj.alpha_m = (2*rho-1)/(rho+1);
                obj.beta    = 0.25*(1-obj.alpha_m+obj.alpha_f)^2;
                obj.gamma   = 0.5 + obj.alpha_f - obj.alpha_m;
            elseif length(genalpha)==4
                obj.beta = genalpha(1);
                obj.gamma = genalpha(2);
                obj.alpha_m = genalpha(3);
                obj.alpha_f = genalpha(4);
            end
            if obj.beta == 0, obj.explicit = true; end
        end
        
        %Entire functions
        function  [u_np1] = ExplicitStep(obj,prob,u_n,t_np1,precomp)
            
            if nargin < 5, precomp=false; end;
            
            u_np1 = u_n + obj.dt*prob.prevsv.v + 0.5*obj.dt*obj.dt*prob.prevsv.a;
            
            u_int = (1-obj.alpha_f)*u_np1 + obj.alpha_f*u_n;
            t_int = (1-obj.alpha_f)*t_np1 + obj.alpha_f*(t_np1-obj.dt);
            if precomp
                f=-prob.ResJacPrecomp(u_int,t_int,1);
            else
                f = -prob.ResJac(u_int,t_int);
            end
            
            a_np1 = (1/(1-obj.alpha_m))*(prob.M\f - obj.alpha_m*prob.prevsv.a);
            v_np1 = prob.prevsv.v + obj.dt*(1-obj.gamma)*prob.prevsv.a + obj.dt*obj.gamma*a_np1;
            
            prevstates.v = v_np1;
            prevstates.a = a_np1;
            prob.setPrevSV(prevstates);
        end
        
        function  [R,J,prevstates] = TimeIntNLFunc(obj,prob,U,Uprev,t,precomp)
            %This function sets up parameters to be used in implicit Newmark calculation
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
            %               Implicit Newmark scheme (wrapped in the nonlinear solver)
            %J            - the jacobian of the nonlinear function, adjusted for the
            %               Implicit Newmark scheme (wrapped in the nonlinear solver)
            %--------------------------------------------------------------------------
            
            if nargin < 6, precomp=false; end;
                        
            %Determine the Residual and Jacobian (function handles) of the continuous
            %system of equations.  f = fint - fext.
            u_tilde_np1 = Uprev + obj.dt*prob.prevsv.v + 0.5*obj.dt*obj.dt*(1-2*obj.beta)*prob.prevsv.a;
            v_tilde_np1 = prob.prevsv.v + (1-obj.gamma)*obj.dt*prob.prevsv.a;
            
            a_np1 = (1/(obj.beta*obj.dt*obj.dt))*(U - u_tilde_np1);
            v_np1 = v_tilde_np1 + obj.dt*obj.gamma*a_np1;
            
            t_int = (1-obj.alpha_f)*t + obj.alpha_f*(t-obj.dt);
            a_int = (1-obj.alpha_m)*a_np1 + obj.alpha_m*prob.prevsv.a;
            %v_int = (1-obj.alpha_f)*v_np1 + obj.alpha_f*prob.prevsv.v;
            u_int = (1-obj.alpha_f)*U + obj.alpha_f*Uprev;
            
            if precomp
                [f,df_internal] = prob.ResJacPrecomp(u_int,t_int,1);
            else
                [f,df_internal] = prob.ResJac(u_int,t_int);
            end

            prevstates.a = a_np1;
            prevstates.v = v_np1;
            
            R = prob.M*a_int + f;
            J = ((1-obj.alpha_m)/(obj.beta*obj.dt*obj.dt))*prob.M + (1-obj.alpha_f)*df_internal;
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
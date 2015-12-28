classdef SQPold < handle
    properties
        id = [];
        type = 'nonlinearly constrained';
        method = 'SQP';
        QNupdate = 'BFGS';
        QPsubprob = 'ActiveSet';
        MeritFunc = 'L1penalty';
        optFrame = 'NAND';
        sensMethod = 'adjoint';
        
        optStrategy = 'LS';
        LStype = 'backtrack'; %zoom or newton or backtrack
        LSprop = 0.75; %0<tau<1 for backtrack,
        
        gradTol = 1e-6;
        dpTol = 1e-6;
        objective;
        constraints;
        
        E;
        I;
        mu;
        
        PDEstartLS  = 'reset';
        PDEstartNew = 'warm';
        
        paramBounds = [];
        paramHist = [];
        solnHist = [];
        nSolverCalls = 0;
        PDECallsParam = 0;
        
        model;
        time = 0;
        printConv = true;
        dp;
    end
      
    methods
        function [obj] = SQPold(OPTfile,optid,modelobj)
            %This is the constructor of the QuasiNewton class.  It reads the
            %appropriate entries from the input file and stores them in the
            %class instance.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %OPTfile  - string indicating the filename of the opt file to
            %           use
            %optid    - integer indicating the OPT object in OPTfile to
            %           use
            %modelobj - instance of a model object 
            %
            %Outputs:
            %--------
            %obj     - instance of the OPT object that was constructed
            %--------------------------------------------------------------
            
            %obj.I = 1;
            obj.I = [1:modelobj.ndof-1]';
            obj.E = [];
            
            %Set model object as a property
            obj.model = modelobj;
            %Extract the OPT text from the cfg file
            OPTtxt = readInFile('OPT',OPTfile,optid);
            %Determine the OPT properties based on the text in the OPT file
            obj.determineOPT(OPTtxt,[modelobj.GVARtxt; readInFile('VAR',OPTfile,1)]);
        end
        
        function  [] = determineOPT(obj,OPTtext,VARtxt)
            %This function computes interprets the text in OPTtext
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - instance of OPT class
            %OPTtxt  - a string containing the text from the OPT input file
            %VARtxt  - a string containing the text from the VAR block in
            %          the OPT file
            %
            %Outputs:
            %--------
            %No outputs.  Properties are stored in OPT handle class.
            %--------------------------------------------------------------
            
            %Extract variables from VARtxt and evaluate them in the current
            %workspace.
            if ~isempty(VARtxt)
                for i = 1:length(VARtxt)
                    eval([VARtxt{i},';']);
                end
            end
            
            %Determine ID number from OPT file
            obj.id = extractInputRobust(OPTtext,'id',[]);
            if isempty(obj.id)
                error('Need to specify id in OPT file');
            end
            %Determine the type of optimization problem
            obj.type = extractInputRobust(OPTtext,'type',obj.type);
            %Determine the details of the optimization algorithm
            details = extractInputRobust(OPTtext,'details',[]);
            details(isspace(details)) = [];
            ind = regexp(details,'-');
            obj.optStrategy = details(1:ind(1)-1);
            obj.MeritFunc = details(ind(1)+1:ind(2));
            obj.QNupdate = details(ind(2)+1:end);
            %Determine the details of the optimization algorithm
            obj.optFrame = extractInputRobust(OPTtext,'optFramework',obj.optFrame);
            %Determine the optimization strategy (linesearch or trust region)
            obj.optStrategy = extractInputRobust(OPTtext,'optStrategy',obj.optStrategy);
            %Determine linesearch or trust region properties
            switch upper(obj.optStrategy)
                case 'LS'
                    obj.LStype = extractInputRobust(OPTtext,'LStype',obj.LStype);
                    obj.LSprop = extractInputRobust(OPTtext,'LSprop',obj.LSprop);
                case 'TR'
                    obj.LStype = [];
                    obj.LSprop = [];
            end
            %Determine starting point for PDE solver in different scenarios
            obj.PDEstartLS = extractInputRobust(OPTtext,'PDEstartLS',obj.PDEstartLS);
            obj.PDEstartNew = extractInputRobust(OPTtext,'PDEstartNew',obj.PDEstartNew);
            %Determine whether to use adjoint or direct method when
            %computing state sensitivities
            obj.sensMethod = extractInputRobust(OPTtext,'sensMethod',obj.sensMethod);
            %Determine gradient and search direction tolerances
            obj.gradTol = extractInputRobust(OPTtext,'gradTol',obj.gradTol);
            obj.dpTol = extractInputRobust(OPTtext,'dpTol',obj.dpTol);
            %Determine parameter bounds
            obj.paramBounds = extractInputRobust(OPTtext,'paramBounds',[]);
            %Set the objective and constraint function handles
            objectiveTxt = extractInputRobust(OPTtext,'objective','ParameterEstimation');
            obj.model.setOptObjective(objectiveTxt);
            %obj.objective = eval(['@(u,z) obj.model.prob.',objectiveTxt,'(u,z)']);
            constraintTxt = extractInputRobust(OPTtext,'constraints','none');
            obj.model.setOptConstraints(constraintTxt);
            %obj.constraints = eval(['@(u,z,E,I)',constraintTxt,'(u,z,E,I)']);
            %Determine whether or not to print convergence indicators
            obj.printConv = extractInputRobust(OPTtext,'printConv',obj.printConv);
        end
        
        function [LSfun,LSder] = LineSearchNAND(obj,alpha)
            %This function is the 1d function (of the step size: alpha)
            %used in the linesearch for NAND method
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj   - QuasiNewton object
            %alpha - step size
            %
            %Outputs:
            %--------
            %LSfun - value of this linesearch function for the given alpha
            %LSder - value of the derivative of this linesearch function
            %        for the given alpha
            %--------------------------------------------------------------
            
            p = obj.paramHist(:,end)+alpha*obj.dp;
            
            if alpha ~= 0
                %We only call the PDEsolver again if we have not yet solved
                %the PDE at this value of the parameter.  For alpha = 0, we
                %have already solved the PDE at p + alpha*dp = p on the
                %previous iteration, so we save ourselves a PDE call by
                %only calling the PDE solver if alpha ~= 0
                obj.model.prob.updateParameters(p);
                %Determine which initial guess to use in the PDE
                %solver
                switch obj.PDEstartLS
                    case 'warm'
                        obj.model.resetInitialCondition(obj.solnHist(:,end));
                    case 'reset'
                        obj.model.resetInitialCondition(obj.model.prob.ic);
                end
                obj.model.executeModel;
            end
            %Compute reduced gradient
            %[dcdw,dcdp] = obj.model.prob.ResSens(obj.model.sv(:,2));
            [DfDp,f] = obj.model.reducedGrad(p,obj.sensMethod);
            %[DfDp,f] = obj.reducedGrad(obj.model.sv(:,2),dcdw,dcdp);
            
            %Return relevant outputs
            LSfun = f;
            LSder = DfDp'*obj.dp;
        end
        
        function  [phi,Dphi] = MeritFunction(obj,p,s)
           
            %Set new parameters
            obj.model.prob.updateParameters(p);
            %Run the model to solve the constraints and store the solution
            obj.model.executeModel;
            
            [DfDp,f] = obj.model.reducedGrad(p,obj.sensMethod);
            [~,c]    = obj.model.reducedConstJac(p,obj.sensMethod);
            cbar = c;
            cbar(obj.I,1) = cbar(obj.I,1) + s;
            
            phi = f + obj.mu*norm(cbar,1);
            Dphi = DfDp'*obj.dp - obj.mu*norm(cbar,1);
        end
        
        function  [f] = callObjective(obj,p)
            
            obj.model.prob.updateParameters(p);
            obj.model.executeModel;
            [~,f] = obj.model.reducedGrad(p,obj.sensMethod);
            
        end
        
        function  [c,ceq] = callConstraints(obj,p)
            
            c = -1;
            ceq = [];
            
%             obj.model.prob.updateParameters(p);
%             obj.model.executeModel;
%             [~,cAll] = obj.model.reducedConstJac(p,obj.sensMethod);
%             c = cAll(obj.I);
%             ceq = cAll(~obj.I);
            
        end
        
        function [p] = lineSQP_NAND(obj,p0,lam0)
            %This function runs the SQP linesearch NAND optimization loop.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj - QuasiNewton object
            %p   - Initial point for optimization algorithm 
            %
            %Outputs:
            %--------
            %No outputs.  Data stored in SQP object.
            %--------------------------------------------------------------
            
            func = @(p) obj.callObjective(p);
            cons = @(p) obj.callConstraints(p);
            
            options=optimset('Algorithm','interior-point');
            [p,fval,exitflag,output,lamTemp,grad] = fmincon(func,p0,[],[],[],[],obj.paramBounds(:,1),obj.paramBounds(:,1),cons,options);
            obj.paramHist = p;
            obj.solnHist = obj.model.sv(:,end);
            return;

            tStartPDEopt = tic; %Timer
            eta = 0.05; tau = 0.5;
            
            %Set up the problem residual and jacobian with the appropriate
            %parameter
            p = p0(:);
            m = size(p,1); %Number of parameters
            nIneqCon = length(obj.I); %Number of inequality constraints
            
            obj.model.prob.updateParameters(p);
            %Run the model to solve the constraints and store the solution
            obj.model.executeModel;
            obj.solnHist = obj.model.sv(:,end);
            
            %Compute the gradient of the objective
            [DfDp,f] = obj.model.reducedGrad(p,obj.sensMethod);
            DfDpbar = [DfDp;zeros(nIneqCon,1)];
            [A,c]    = obj.model.reducedConstJac(p,obj.sensMethod);
            cbar = c;
            s = max([zeros(nIneqCon,1),-cbar],[],2);
            cbar(obj.I,1) = cbar(obj.I,1) + s;
            tmp = zeros(size(A,1),nIneqCon);
            tmp(obj.I + nIneqCon*([1:nIneqCon]'-1)) = 1; %Contributions of slack variables to jacobian
            Abar = [A,tmp];
            
            %Record the number of Solver calls and the parameter evolution
            obj.nSolverCalls = 1;
            obj.PDECallsParam = p; %Includes linesearch parameters
            obj.paramHist = p; %Does not include linesearch parameters
            
            %Initial Hessian approximation
            %Hbar = speye(size(DfDpbar,1));
            %H = Hbar(1:m,1:m);
            H = speye(m);
            Hbar = sparse([H, zeros(m,nIneqCon); zeros(nIneqCon,m), zeros(nIneqCon,nIneqCon)]);
            cnt = 0;
            
            %Make initial guess for LM
            lam = (Abar')\(Hbar*[p;s] + DfDpbar);
            
            QPsubprog.lb = [obj.paramBounds(:,1);-ones(nIneqCon,1)];
            QPsubprog.ub = [obj.paramBounds(:,2);inf(nIneqCon,1)];
            QPsubprog.solver = 'quadprog';
            QPsubprog.options = optimset('quadprog');
            while (norm(DfDp,2)>obj.gradTol)
                cnt = cnt+1;
                %Solve QP
                QPsubprog.H = Hbar;
                QPsubprog.f = DfDpbar;
                %QPsubprog.H = [H, zeros(m,nIneqCon); zeros(nIneqCon,m), zeros(nIneqCon,nIneqCon)];
                %QPsubprog.f = [DfDp;zeros(nIneqCon,1)];
                
                %QPsubprog.Aineq = A(obj.I,:);
                %QPsubprog.bineq = -c(obj.I,:);
                
                QPsubprog.Aeq   = Abar;
                QPsubprog.beq = -cbar; 
                
                %QPsubprog.Aeq   = A(obj.E,:);
                %QPsubprog.beq = -c(obj.E,:);
                %[Dx,~,~,~,lamTemp] = quadprog(QPsubprog);
                options=optimset('Algorithm','active-set');
                %[Dx,fval,exitflag,output,lamTemp,grad] = fmincon(@(z) DfDp'*z + 0.5*z'*H*z,p,A,-c,[],[],obj.paramBounds(:,1),obj.paramBounds(:,2),[],options);
%                 [pSub,fval,exitflag,output,lamTemp,grad] = fmincon(@(z) DfDpbar'*z + 0.5*z'*Hbar*z,[p;s],[],[],[],[],QPsubprog.lb,QPsubprog.ub,[],options);
                [pSub,fval,exitflag,output,lamTemp,grad] = fmincon(@(z) DfDpbar'*z + 0.5*z'*Hbar*z,[p;s],[],[],Abar,-cbar,QPsubprog.lb,QPsubprog.ub,[],options);
                Dx = pSub - p;
                
%                 [Dx,fval,exitflag,output,lamTemp,grad] = fmincon(@(z) (Hbar*[p;s] + DfDpbar)'*z + 0.5*z'*Hbar*z,QPsubprog.lb-[p;s],[],[],Abar,zeros(size(cbar)),QPsubprog.lb-[p;s],QPsubprog.ub-[p;s],[],options);
                
                Dx(abs([p;s] - [obj.paramBounds(:,1);zeros(nIneqCon,1)]) < 1e-6 | abs([p;s] - [obj.paramBounds(:,2);inf(nIneqCon,1)])<1e-6) = 0;
%                 Dx([p;s] == [obj.paramBounds(:,1);zeros(nIneqCon,1)] | [p;s] == [obj.paramBounds(:,2);inf(nIneqCon,1)]) = 0;
                
                obj.dp = Dx(1:m);
                ds     = Dx(m+1:end);
                
                lamQPsub = (Abar')\(Hbar*Dx + DfDpbar);
                %lamQPsub = lamTemp.eqlin;
                %lamQPsub = zeros(size(lam));
                %lamQPsub(obj.I) = lamTemp.ineqlin;
                %lamQPsub(obj.E) = lamTemp.eqlin;
                dLam = lamQPsub - lam;
                
                rho = 0.5;
                %obj.mu = 0;
                if ~(obj.mu >= (DfDpbar'*[obj.dp;ds] + 0.5*[obj.dp;ds]'*Hbar*[obj.dp;ds])/((1-rho)*norm(cbar,1)))
                    obj.mu = abs((DfDpbar'*[obj.dp;ds] + 0.5*[obj.dp;ds]'*Hbar*[obj.dp;ds])/((1-rho)*norm(cbar,1)));
                end
                %obj.mu = abs((DfDp'*obj.dp + 0.5*obj.dp'*H*obj.dp)/((1-rho)*norm(cbar,1)));
                
                %This selection of alpha_max ensures that
                %p_low <= p + alpha_max*dp <= p_up.
                %Rearranging this equation, we get that alpha_max must
                %satsify the following two conditions:
                %alpha_max <= ((p_up - p)./dp)_j   if dp_j > 0
                %alpha_max <= ((p - p_low)./dp)_j  if dp_j < 0
                if isempty(obj.paramBounds)
                    alpha_max = min(obj.LSprop(3),1);
                else
                    temp = zeros(size(Dx));
                    pxind = (Dx > 0);
                    ppind = (obj.dp > 0);
                    psind = (ds > 0);
                    if sum(pxind) > 0
                        temp(pxind) = 0.5*([obj.paramBounds(ppind,2);inf(sum(psind),1)] - [p(ppind);s(psind)])./Dx(pxind);
                    end
                    
                    mxind = (Dx < 0);
                    mpind = (obj.dp < 0);
                    msind = (ds < 0);
                    if sum(mxind) > 0
                        temp(mxind) = 0.5*([obj.paramBounds(mpind,1);zeros(sum(msind),1)] - [p(mpind);s(msind)])./Dx(mxind);
                    end
                    
                    zind = (abs(Dx) < 1e-6);
                    temp(zind) = 0.5*1;
                    
                    alpha_max = min(temp);
%                     alpha_max = min(...
%                        (Dx > 0).*([obj.paramBounds(:,2);inf(nIneqCon,1)] - [p;s])./Dx + ...
%                        (Dx < 0).*([obj.paramBounds(:,1);zeros(nIneqCon,1)] - [p;s])./Dx + ...
%                        (Dx == 0));
                    %alpha_max = min(...
                    %    (obj.dp > 0).*(obj.paramBounds(:,2) - p)./obj.dp + ...
                    %    (obj.dp < 0).*(obj.paramBounds(:,1) - p)./obj.dp);
                    alpha_max = min([alpha_max,obj.LSprop(3),1]);
                end
                
                maxit = 5;
                while 1
                    alpha = alpha_max;
                    [phiL,DphiL] = obj.MeritFunction(p,s);
                    if DphiL > 0
                        disp('switching');
                        obj.dp = -obj.dp;
                    end
                    phi  = obj.MeritFunction(p + alpha*obj.dp,s + alpha*ds);
                    
                    it = 0;
                    conv = true;
                    while phi > phiL + eta*alpha*DphiL
                        if it == maxit
                            obj.mu = 2*obj.mu;
                            conv = false;
                            %conv = true;
                            break;
                        end
                        it = it + 1;
                        
                        alpha = tau*alpha
                        phi  = obj.MeritFunction(p + alpha*obj.dp,s+alpha*ds);
                        [phi, phiL + eta*alpha*DphiL]
                    end
                    
                    if conv 
                        break;
                    end
                end
                
                dpConv = norm(obj.dp)/norm(p);
                if obj.printConv
                    fprintf('\n\n Iteration %d:\n',cnt);
                    disp(['Gradient Norm = ', num2str(norm(DfDp,2))]);
                    disp(['Step Length = ', num2str(alpha)]);
                    disp(['Max. Step Length = ', num2str(alpha_max)]);
                    disp(['||dp||/||p|| = ', num2str(dpConv)]);
                end
                
                Dx = alpha*Dx;
                obj.dp = alpha*obj.dp;
                ds = alpha*ds;
                dLam = alpha*dLam;
                
                p = p + obj.dp;
                s = s + ds;
                lam = lam + dLam;
                
                obj.paramHist = [obj.paramHist,p];
                
                [~,dLdpOld] = obj.Lagrangian(lam,f,DfDpbar,cbar,Abar);
                obj.model.prob.updateParameters(p);
                %Run the model to solve the constraints and store the solution
                obj.model.executeModel;
                obj.solnHist = obj.model.sv(:,end);
                %Compute the gradient of the objective
                [DfDp,f] = obj.model.reducedGrad(p,obj.sensMethod);
                DfDpbar  = [DfDp;zeros(nIneqCon,1)];
                [A,c]    = obj.model.reducedConstJac(p,obj.sensMethod);
                cbar = c;
                cbar(obj.I,1) = cbar(obj.I,1) + s;
                tmp = zeros(size(A,1),nIneqCon);
                tmp(obj.I + nIneqCon*([1:nIneqCon]'-1)) = 1; %Contributions of slack variables to jacobian
                Abar = [A,tmp];
                
                [~,dLdp] = obj.Lagrangian(lam,f,DfDpbar,cbar,Abar);
                
                ss = Dx; y = dLdp - dLdpOld;
                Hbar = updateHessian(Hbar,y,ss,obj.QNupdate);
                H = Hbar(1:m,1:m);
            end
            obj.time = toc(tStartPDEopt);
            
            
%             %Set up the problem residual and jacobian with the appropriate
%             %parameter
%             p = p(:);
%             obj.model.prob.updateParameters(p);
%             %Run the model to solve the constraints and store the solution
%             obj.model.executeModel;
%             obj.solnHist = obj.model.sv(:,end);
%             %Compute the residual sensitivities
%             %[dcdw,dcdp] = obj.model.prob.ResSens(obj.solnHist);
%             
%             %Record the number of Solver calls and the parameter evolution
%             obj.nSolverCalls = 1;
%             obj.PDECallsParam = p; %Includes linesearch parameters
%             obj.paramHist = p; %Does not include linesearch parameters
%             
%             %Compute the gradient of the objective
%             DfDp = obj.model.reducedGrad(p,obj.sensMethod);
%             %DfDp = obj.reducedGrad(obj.solnHist(:,end),dcdw,dcdp);
%             %Initial Hessian approximation
%             H = speye(size(DfDp,1));
%             cnt = 0;
%             while (norm(DfDp,2)>obj.gradTol)
%                 cnt = cnt+1;
%                 obj.dp = -H*DfDp; % QN updates on inverse Hessian
%                 
%                 %This selection of alpha_max ensures that
%                 %p_low <= p + alpha_max*dp <= p_up.
%                 %Rearranging this equation, we get that alpha_max must
%                 %satsify the following two conditions:
%                 %alpha_max <= ((p_up - p)./dp)_j   if dp_j > 0
%                 %alpha_max <= ((p - p_low)./dp)_j  if dp_j < 0
%                 if isempty(obj.paramBounds)
%                     alpha_max = min(obj.LSprop(3),1);
%                 else
%                     alpha_max = min(...
%                         (obj.dp > 0).*(obj.paramBounds(:,2) - p)./obj.dp + ...
%                         (obj.dp < 0).*(obj.paramBounds(:,1) - p)./obj.dp);
%                     alpha_max = min([alpha_max,obj.LSprop(3),1]);
%                 end
%                 
%                 switch lower(obj.LStype)
%                     case 'newton'
%                         alpha = alpha_max;
%                         %Set up the problem residual and jacobian with the appropriate
%                         %parameter
%                         obj.model.prob.updateParameters(p);
%                         %Determine which initial guess to use in the PDE
%                         %solver
%                         switch obj.PDEstartNew
%                             case 'warm'
%                                 obj.model.resetInitialCondition(obj.solnHist(:,end));
%                             case 'reset'
%                                 obj.model.resetInitialCondition(obj.model.prob.ic);
%                         end
%                         %Run the model to solve the constraints
%                         obj.model.executeModel;
%                         
%                         obj.PDECallsParam = [obj.PDeCallsParam, p+alpha*obj.dp];
%                         obj.nSolverCalls = obj.nSolverCalls + 1;
%                     case 'zoom'
%                         [alpha,~,nFunCalls,alphaCalls] = linesearch(@(alphaloc) obj.LineSearchNAND(alphaloc),[obj.LSprop(1:2);alpha_max;obj.LSprop(4:end)]);
%                         obj.PDECallsParam = [obj.PDECallsParam,repmat(p,1,size(alphaCalls,2))+repmat(obj.dp,1,size(alphaCalls,2)).*repmat(alphaCalls,size(obj.dp,1),1)];
%                         obj.nSolverCalls = obj.nSolverCalls + nFunCalls;
%                 end
%                 
%                 %Compute step, update iterate, and update parameter history
%                 obj.dp = alpha*obj.dp;
%                 p = p + obj.dp;
%                 obj.paramHist = [obj.paramHist,p];
%                 dpConv = abs(norm(obj.dp)/norm(p));
%                 
%                 if obj.printConv
%                     fprintf('\n\n Iteration %d:\n',cnt);
%                     disp(['Gradient Norm = ', num2str(norm(DfDp,2))]);
%                     disp(['Step Length = ', num2str(alpha)]);
%                     disp(['Max. Step Length = ', num2str(alpha_max)]);
%                     disp(['||dp||/||p|| = ', num2str(dpConv)]);
%                 end
%                 
%                 %Break when search direction is small compared to current
%                 %iterate (i.e. we are taking a very small step compared to
%                 %our size)
%                 if (dpConv<obj.dpTol)
%                     break;
%                 end
%                 
%                 %We don't run the model at the new value of p because this
%                 %has ALREADY BEEN DONE IN THE LINESEARCH.  It would be a
%                 %waste to do it again.  Therefore, we simply update the
%                 %values of paramHist and solnHist according to this LAST
%                 %VALUE OF alpha from the linesearch.
%                 obj.solnHist = [obj.solnHist,obj.model.sv(:,2)];
%                 %[dcdw,dcdp] = obj.model.prob.ResSens(obj.solnHist(:,end));
%                 
%                 %Record previous gradient and get the updated gradient
%                 DfDp_old = DfDp;
%                 DfDp = obj.model.reducedGrad(p,obj.sensMethod);
%                 %DfDp = obj.reducedGrad(obj.solnHist(:,end),dcdw,dcdp);
%                 %Use gradients and search directions to update the inverse
%                 %Hessian
%                 s = obj.dp;
%                 y = DfDp - DfDp_old;
%                 H = updateInverseHessian(H,y,s,obj.QNupdate);
%             end
%             obj.time = toc(tStartPDEopt);
        end
        
        function  [L,dLdp] = Lagrangian(obj,lam,f,DfDp,c,A)
            %This function returns the Lagrangian of the current problem.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj - SQP object
            %p   - parameter vector
            %lam - vector of Lagrange multipliers
            %
            %Outputs:
            %--------
            %No outputs.  Data stored in SQP object.
            %--------------------------------------------------------------
            
            %[DfDp,f] = obj.model.reducedGrad(p,obj.sensMethod);
            %[A,c] = reducedConstJac(obj,p,obj.sensMethod);
            
            L = f - lam'*c;
            dLdp = DfDp - A'*lam;
        end
                        
        function  [] = executeOptimization(obj,p)
            %This function calls the appropriate QuasiNewton
            %PDE-constrained optimization algorithm
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj - QuasiNewton object
            %
            %Outputs:
            %--------
            %No outputs.  Data stored in QuasiNewton object.
            %--------------------------------------------------------------
            
            switch obj.optFrame
                case 'NAND'
                    %NAND optimization framework
                    if strcmpi(obj.optStrategy,'LS')
                        %Line search strategy
                        lineSQP_NAND(obj,p);
                    elseif strcmpi(obj.optStrategy,'TR')
                        %Trust region strategy
                        trustSQP_NAND(obj,p);
                    end
                case 'SAND'
                    %SAND optimization framework
                    if strcmpi(obj.optStrategy,'LS')
                        %Line search strategy
                        lineSQP_SAND(obj,p);
                    elseif strcmpi(obj.optStrategy,'TR')
                        %Trust region strategy
                        trustSQP_SAND(obj,p);
                    end
            end
        end
    end
end
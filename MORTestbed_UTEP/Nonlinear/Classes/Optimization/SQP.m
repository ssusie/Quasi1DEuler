classdef SQP < handle
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
        LSprop = [10^(-4);0.9;1;1000;1000;0.05;0.1;1.2]; %default is for zoom, not needed for newton, 0<tau<1 for backtrack,
        %LSprop = 0.75; %0<tau<1 for backtrack,
        
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
        function [obj] = SQP(OPTfile,optid,modelobj)
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
%             obj.I = [1:31]';
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
            cbar(obj.I,1) = cbar(obj.I,1) - s;
            
            phi = f + obj.mu*norm(cbar,1);
            Dphi = DfDp'*obj.dp - obj.mu*norm(cbar,1);
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

            tStartPDEopt = tic; %Timer
            %eta = 0.01; tau = 0.5; %Linesearch parameters
            %eta = 0.5; tau = 0.75; %Linesearch parameters
            
            %Set up the problem residual and jacobian with the appropriate
            %parameter
            p = p0(:);
            m = size(p,1); %Number of parameters
            nIneqCon = length(obj.I); %Number of inequality constraints
            
            %Set the parameter inside the problem object
            obj.model.prob.updateParameters(p);
            %Run the model to solve the constraints and store the solution
            obj.model.executeModel;
            obj.solnHist = obj.model.sv(:,end);
            
            %Compute the gradient of the objective and constraints
            [DfDp,f] = obj.model.reducedGrad(p,obj.sensMethod);
            [A,c]    = obj.model.reducedConstJac(p,obj.sensMethod);
            %Set slack variables and include them in the objective and
            %jacobian
            cbar = c;
            s = max([zeros(nIneqCon,1),c],[],2);
            cbar(obj.I,1) = cbar(obj.I,1) - s;
            tmp = zeros(size(A,1),nIneqCon);
            tmp(obj.I + nIneqCon*([1:nIneqCon]'-1)) = -1; %Contributions of slack variables to jacobian
            
            Abar = [A,tmp];
            DfDpbar = [DfDp;zeros(nIneqCon,1)];
            
           
            
%             n = m + nIneqCon;
%             lb = [obj.paramBounds(:,1);zeros(nIneqCon,1)];
%             ub = [obj.paramBounds(:,2);inf(nIneqCon,1)];
%             Asbnd = sparse([-eye(n);eye(n)]);
%             bsbnd = [-ub;lb];
%             act = determineActiveSet([p;s],Asbnd,bsbnd);
%             %act = determineActiveSet([p;s],lb,ub);
%             nact = sum(act);
            
            %Record the number of Solver calls and the parameter evolution
            obj.nSolverCalls = 1;
            obj.PDECallsParam = p; %Includes linesearch parameters
            obj.paramHist = p; %Does not include linesearch parameters
            
            %Initial Hessian approximation
            H = eye(m);
            Hbar = [H, zeros(m,nIneqCon); zeros(nIneqCon,m), zeros(nIneqCon,nIneqCon)];
            cnt = 0;
            
            %Make initial guess for LM
            %lam = (Abar')\(Hbar*[p;s] + DfDpbar);
            %lam = ([A,zeros(size(A,1),1)]')\([H*p;0] + [DfDp;0]);
            
            
           
            
%             QPsubprog.lb = [obj.paramBounds(:,1);zeros(nIneqCon,1)];
%             QPsubprog.ub = [obj.paramBounds(:,2);inf(nIneqCon,1)];
%             QPsubprog.solver = 'quadprog';
%             QPsubprog.options = optimset('quadprog');
            
           
            Np = size(p,1);
            Neq = length(obj.E);
            Nineq = length(obj.I);
            
            lb = obj.paramBounds(:,1);
            ub = obj.paramBounds(:,2);
            Asbnd = [-eye(Np),zeros(Np,Nineq);eye(Np),zeros(Np,Nineq);zeros(Nineq,Np),eye(Nineq)];
            bsbnd = [-ub;lb;zeros(Nineq,1)];
%             Asbnd = [-eye(Np);eye(Np)];
%             bsbnd = [-ub;lb];
            lam = zeros(Neq+Nineq,1);
            [~,dLdp] = Lagrangian(obj,lam,f,DfDp,c,A);
            %Nc = length(obj.I) + length(obj.E);
            
%             lb = [obj.paramBounds(:,1);zeros(Ns,1)];
%             ub = [obj.paramBounds(:,2);inf(Ns,1)];
            
            while (norm([dLdp;c],2) > obj.gradTol)
                cnt = cnt+1;
                
                %Aeq   = A(obj.E,:);
                %beq   = -c(obj.E,1);
                %Aineq = [A(obj.I,:);Asbnd];
                %bineq = [-c(obj.I,1);bsbnd+[p;-p]];
                
                Aeq = [A(obj.E,:),zeros(Neq,Nineq);A(obj.I,:),-eye(Nineq)];
                beq = -[c(obj.E);c(obj.I)-s];
                Aineq = Asbnd;
                bineq = bsbnd-[-p;p;s];
                
                [Dx,lamEq,lamIneq,act] = obj.solveQPsub(Hbar,DfDpbar,Aeq,beq,Aineq,bineq,Np+Nineq,Neq+Nineq,2*Np+Nineq);
                %[obj.dp,lamTemp,act] = obj.solveQPsub(H,DfDp,Aeq,beq,Aineq,bineq,Np,Neq,Nineq);
                %[Dx,lamTemp,act] = solveQPsub(obj,Hbar,DfDpbar,Abar,Asbnd,-cbar,bsbnd,lb - p,ub - p,Np,Nc);
                
                
                actPrimal = logical(sum(reshape(act(1:2*Np,1),Np,2),2));
                actSlack  = act(2*Np+1:end,1);
                
                lamNL = lamEq;
                
                
%                 nact = sum(act);
%                 nactSB = sum(act(end-(2*Np-1):end) > 0);
%                 
%                 actNL = false(size(act));
%                 actNL(1:Neq+Nineq) = act(1:Neq+Nineq);
%                 
%                 lamSB = lamTemp(end-(nactSB-1):end);
%                 lamNL = lamTemp(1:end-nactSB);
                %lamSB = lamTemp(Neq+nact+1:end);
                
                
%                 alpha_max = min(...
%                       (obj.dp > 0).*(ub - p)./obj.dp + ...
%                       (obj.dp < 0).*(lb - p)./obj.dp);
%                 alpha_max = min([alpha_max,obj.LSprop(3),1]);
                
                alpha_max = min(...
                   (Dx > 0).*([ub;inf(Nineq,1)] - [p;s])./Dx + ...
                   (Dx < 0).*([lb;zeros(Nineq,1)] - [p;s])./Dx);
                alpha_max = min([alpha_max,obj.LSprop(3),1]);
                
                
%                 %Solve QP
%                 QPsubprog.H = Hbar;
%                 QPsubprog.f = DfDpbar;
%                 
%                 %QPsubprog.Aineq = A(obj.I,:);
%                 %QPsubprog.bineq = -c(obj.I,:);
%                 
%                 QPsubprog.Aeq   = Abar;
%                 QPsubprog.beq = -cbar; 
%                 
%                 %QPsubprog.Aeq   = A(obj.E,:);
%                 %QPsubprog.beq = -c(obj.E,:);
%                 [Dx,~,~,~,lamTemp] = quadprog(QPsubprog);
                
                
                
%                 %Solve QP
%                 Asbnd_hat = [Abar;Asbnd(act,:)];
%                 Asbnd_bar = Asbnd(~act,:);
%                 bsbnd_bar = bsbnd(~act,:);
%                 tmp = -[Hbar,Asbnd_hat';Asbnd_hat,zeros(nact+nIneqCon,nact+nIneqCon)]\[DfDpbar;cbar;zeros(nact,1)];
%                 Dx = tmp(1:n,1);
%                 lamNL = -tmp(n+1:n+nIneqCon);
%                 lamSB = -tmp(n+nIneqCon+1:end);
%                 %Step direction
%                 %dx = -H*df;
                
%                 gamma = zeros(n - nact,1);
%                 ind = Asbnd_bar*Dx <= 0;
%                 gamma(ind) = inf;
%                 gamma(~ind) = (bsbnd_bar(~ind,1) - Asbnd_bar(~ind,:)*[p;s])./(Asbnd_bar(~ind,:)*Dx);
%                 
%                 if sum(~ind) == 0
%                     alpha_max = inf;
%                 else
%                     alpha_max = min(abs(gamma));
%                 end
%                 alpha_max = min(1,alpha_max);
                
                %obj.dp = Dx(1:m);
                %ds     = Dx(m+1:end);
                
                
                obj.dp = Dx(1:Np);
                ds = Dx(Np+1:Np+Nineq);
                dLam = lamNL - lam;
                %dLam = zeros(Neq+Nineq,1);
                %dLam(actNL,1) = lamNL - lam(actNL,1);
                
                rho = 0.5;
                if norm(cbar,1) == 0
                    obj.mu = 1;
                else
                    obj.mu = abs((DfDpbar'*[obj.dp;ds] + 0.5*[obj.dp;ds]'*Hbar*[obj.dp;ds])/((1-rho)*norm(cbar,1)));
                end
%                 rho = 0.5;
%                 if norm(c(actNL),1) == 0
%                     obj.mu = 1;
%                 else
%                     obj.mu = abs((DfDp'*obj.dp + 0.5*obj.dp'*H*obj.dp)/((1-rho)*norm(c(actNL),1)));
%                 end
                
       
                tau = 0.5;
                eta = 1e-3;
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
%                         if it == maxit
%                             obj.mu = obj.mu/2;
%                             conv = false;
%                             %conv = true;
%                             break;
%                         end
                        it = it + 1;
                        
                        alpha = tau*alpha
                        phi  = obj.MeritFunction(p + alpha*obj.dp,s+alpha*ds);
                        [phi, phiL + eta*alpha*DphiL]
                    end
                    
                    if conv
                        break;
                    end
                end

%                 [alpha,success] = linesearch(@(alphaloc) obj.MeritFunction(p+alphaloc*obj.dp,s+alphaloc*ds),[obj.LSprop(1:2);alpha_max;obj.LSprop(4:end)]);
                
                %if ~success
                %    alpha = alpha_max/20;
                %end
                
                dpConv = norm(obj.dp)/norm(p);
                if obj.printConv
                    fprintf('\n Iteration %d:\n',cnt);
                    disp(['Gradient Norm = ', num2str(norm([dLdp;c],2))]);
                    disp(['Step Length = ', num2str(alpha)]);
                    disp(['Max. Step Length = ', num2str(alpha_max)]);
                    disp(['||dp||/||p|| = ', num2str(dpConv)]);
                    disp(['p = ', num2str(p(:)')]);
                end
                
                %Dx = alpha*Dx;
                obj.dp = alpha*obj.dp;
                %ds = alpha*ds;
                dLam = alpha*dLam;
                
                %pOld = p;
                p = p + obj.dp;
                s = s + ds;
                lam = lam + dLam;
                obj.paramHist = [obj.paramHist,p];
                
%                 %Remove constraints from active set if they have negative
%                 %lagrange multiplier (simple bounds only)
%                 if  sum(lamSB < 0) > 0
%                     act = removeConstraint(act,lamSB);
%                     nact = sum(act);
%                 end
%                 
%                 %Add constraints to the active set when we step to them
%                 %(simple bounds only)
%                 if alpha_max < 1
%                     act = addConstraint(act,gamma,alpha);
%                     nact = sum(act);
%                 end

                %Break when search direction is small compared to current
                %iterate (i.e. we are taking a very small step compared to
                %our size)
                if (dpConv<obj.dpTol)
                    break;
                end
                
                %[~,dLdpOld] = obj.Lagrangian(lam,f,DfDp,c,A);
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
                cbar(obj.I,1) = cbar(obj.I,1) - s;
                tmp = zeros(size(A,1),nIneqCon);
                tmp(obj.I + nIneqCon*([1:nIneqCon]'-1)) = -1; %Contributions of slack variables to jacobian
                Abar = [A,tmp];
                %[~,dLdp] = obj.Lagrangian(lam,f,DfDp,c,A);
                [~,dLdp] = obj.Lagrangian(lam,f,DfDpbar,cbar,Abar);
                
                ss = Dx; y = dLdp - dLdpOld;
                Hbar = updateHessian(Hbar,y,ss,obj.QNupdate);
                H = Hbar(1:m,1:m);
            end
            obj.time = toc(tStartPDEopt);
        end
        
        function  [z,lamEq,lamIneq,ActCon] = solveQPsub(obj,G,c,Aeq,beq,Aineq,bineq,Nx,Neq,Nineq)
            %Nx = number of primal (non-slack) variable
            %Nc  = number of constraints that are not simple bound
            
            
            %Determine initial point
            z = linprog(zeros(Nx,1), -Aineq, -bineq, Aeq, beq);
            %x = A\b;
            
            %Determine active set at initial point
            act = determineActiveSet(z,Aineq,bineq);
            nact = sum(act);
            
            while 1
                %Set up matrices defining the constraints
                Abar = [Aeq;Aineq(act,:)];
                %bbar = [b;bbnd(act,:)];
                
                %Solve KKT system
                tmp = [G,Abar';Abar,zeros(size(Abar,1))]\[-G*z - c;zeros(Neq+nact,1)];
                
                %Extract components of KKT vector
                q = tmp(1:Nx,1);
                lamEq = -tmp(Nx+1:Nx+Neq);
                lamIneq = -tmp(Nx+Neq+1:end);
                
                if norm(q) < 1e-6
                    %Remove constraints from active set if they have negative
                    %lagrange multiplier (simple bounds only)
                    if  sum(lamIneq < 0) > 0
                        act = removeConstraint(act,lamIneq);
                        nact = sum(act);
                    else
                        %lam = [lamEq;lamIneq];
                        ActCon = [true(Neq,1);act];
                        return;
                    end
                else
                    %Determine maximum step length
                    gamma = zeros(Nineq,1);
                    ind = Aineq*q >= 0 | act;
                    gamma(ind) = inf;
                    gamma(~ind) = (bineq(~ind,1) - Aineq(~ind,:)*z)./(Aineq(~ind,:)*q);
                    %gamma = zeros(Nineq,1);
                    %ind = Aineq(~act,:)*q >= 0;
                    %gamma(ind) = inf;
                    %gamma(~ind) = (bineq(~ind,1) - Aineq(~ind,:)*z)./(Aineq(~ind,:)*q);
                    
                    alpha = min(1,min(gamma));
                    
                    z = z + alpha*q;
                    
                    if alpha <= 1
                        act = addConstraint(act,gamma,alpha);
                        nact = sum(act);
                    end
                end
            end
        end

        
%         function  [x,lam,act] = solveQPsub(obj,H,g,A,Abnd,b,bbnd,lb,ub,x,Nx,Ns,Nc)
%             %Nx = number of primal (non-slack) variable
%             %Ns = number of slack variables (=nIneqCon)
%             %Nc  = number of constraints that are not simple bound
%             
%             
%             %Determine initial point
%             p = linprog(zeros(Nx+Ns,1), [], [], A, b, lb - x, ub - x);
%             %x = A\b;
%             
%             %Determine active set at initial point
%             act = determineActiveSet(p,Abnd,bbnd);
%             nact = sum(act);
%             
%             while 1
%                 %Set up matrices defining the constraints
%                 Abar = [A;Abnd(act,:)];
%                 %bbar = [b;bbnd(act,:)];
%                 
%                 %Solve KKT system
%                 tmp = [H,Abar';Abar,zeros(size(Abar,1))]\[-g - H*p;b - A*p;bbnd(act,:)];
%                 
%                 %Extract components of KKT vector
%                 v = tmp(1:Nx+Ns,1);
%                 lam = -tmp(Nx+Ns+1:Nx + Ns + Nc);
%                 lamSB = -tmp(Nx + Ns + Nc + 1:end);
%                 
%                 if norm(v) == 0
%                     %Remove constraints from active set if they have negative
%                     %lagrange multiplier (simple bounds only)
%                     if  sum(lamSB < 0) > 0
%                         act = removeConstraint(act,lamSB);
%                         nact = sum(act);
%                     else
%                         return;
%                     end
%                 else
%                     %Determine maximum step length
%                     gamma = zeros(2*(Nx + Ns),1);
%                     ind = Abnd(~act,:)*v >= 0;
%                     gamma(ind) = inf;
%                     gamma(~ind) = (bbnd(~ind,1) - Abnd(~ind,:)*p)./(Abnd(~ind,:)*v);
%                     
%                     alpha = min(1,min(gamma));
%                     
%                     p = p + alpha*v;
%                     
%                     if alpha <= 1
%                         act = addConstraint(act,gamma,alpha);
%                         nact = sum(act);
%                     end
%                 end
%             end
%         end
        
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
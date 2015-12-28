classdef BoxConstMin < handle
    properties
        id = [];
        type = 'bound constrained';
        optStrategy = 'LS';
        method = 'QN';
        details = 'wolfe';
        secondOrder = 'BFGS';
        optFrame = 'NAND';
        sensMethod = 'adjoint';
        meritFunc = 'L1';
        
        LSprop = [10^(-4);0.9;1;1000;1000;0.05;0.1;1.2]; %default is for zoom, not needed for newton, 0<tau<1 for backtrack,
        
        gradTol = 1e-6;
        dpTol = 1e-6;
        %objective;
        %constraints;
        
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
        dps;
        
        numModelEval=0;
    end
    
    methods
        function [obj] = BoxConstMin(OPTfile,optid,modelobj)
            %This is the constructor of the BoxConstMin class.  It reads the
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
            %Determine the optimization strategy (linesearch or trust region)
            obj.optStrategy = extractInputRobust(OPTtext,'optStrategy',obj.optStrategy);
            %Determine the method that will be used to get a search
            %direction
            obj.method = extractInputRobust(OPTtext,'method',obj.method);
            %Determine the details of the optimization algorithm
            obj.details = extractInputRobust(OPTtext,'details',obj.details);
            %Determine what second order information to use
            obj.secondOrder = extractInputRobust(OPTtext,'secondOrder',obj.secondOrder);
            %Determine the merit function to use (if applicable)
            obj.meritFunc = extractInputRobust(OPTtext,'meritFunc',obj.meritFunc);
            %             details(isspace(details)) = [];
            %             ind = regexp(details,'-');
            %             obj.optStrategy = details(1:ind-1);
            %             obj.QNupdate = details(ind+1:end);
            %Determine the details of the optimization algorithm
            obj.optFrame = extractInputRobust(OPTtext,'optFramework',obj.optFrame);
            %Determine linesearch or trust region properties
            switch upper(obj.optStrategy)
                case 'LS'
                    %obj.LStype = extractInputRobust(OPTtext,'LStype',obj.LStype);
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
%             constraintTxt = extractInputRobust(OPTtext,'constraints','none');
%             obj.model.setOptConstraints(constraintTxt);
            %obj.constraints = eval(['@(u,z,E,I)',constraintTxt,'(u,z,E,I)']);
            %Determine whether or not to print convergence indicators
            obj.printConv = extractInputRobust(OPTtext,'printConv',obj.printConv);
        end
        
        function [LSfun,LSder] = LineSearchReducedSpace(obj,w,alpha)
            %This function is the 1d function (of the step size: alpha)
            %used in the linesearch for NAND method
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj   - BoxConstMin object
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
                        %obj.model.resetInitialCondition(obj.solnHist(:,end));
                        obj.model.resetInitialCondition(w);
                    case 'reset'
                        obj.model.resetInitialCondition(obj.model.prob.ic);
                end
                obj.model.executeModel;
                w = obj.model.sv(:,end);
            end
            %Compute reduced gradient
            %[~,dcdw,dcdp] = obj.model.prob.ResSens(obj.model.sv(:,2));
            [DfDp,f] = obj.model.reducedGrad(w,p,obj.sensMethod);
            %[DfDp,f] = obj.reducedGrad(obj.model.sv(:,2),dcdw,dcdp);
            
            %Return relevant outputs
            LSfun = f;
            LSder = DfDp'*obj.dp;
        end
        
        function [LSfun,LSder] = LineSearchFullSpace(obj,p,w,lambda,dw,dlambda,alpha,phi)
            %This function is the 1d function (of the step size: alpha)
            %used in the linesearch for NAND method
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj   - BoxConstMin object
            %alpha - step size
            %
            %Outputs:
            %--------
            %LSfun - value of this linesearch function for the given alpha
            %LSder - value of the derivative of this linesearch function
            %        for the given alpha
            %--------------------------------------------------------------
            
            lambda_new =  lambda+alpha*dlambda;
            p_new = p+alpha*obj.dp;
            w_new = w + alpha*dw;
            
            obj.model.prob.updateParameters(p_new);
            [f,dfdw,dfdp] = obj.model.objective(w_new, p_new);
            [c,dcdw,dcdp] = obj.model.PDEconstraint(w_new, p_new);
            
            switch obj.meritFunc
                case 'QuadPen'
                    LSfun =  f +  phi*norm(c,2)^2;
                    LSder = dfdw'*dw + dfdp'*obj.dp + 2*phi*c'*(dcdw*dw + dcdp*obj.dp);
                case 'AugLag'
                    LSfun =  f + lambda_new'*c + phi*norm(c,2)^2;
                    LSder = (dfdw'*dw + lambda_new'*dcdw*dw)  +...
                        (dfdp' + lambda_new'*dcdp)*obj.dp +...
                        c'*dlambda + 2*phi*c'*(dcdw*dw  +dcdp*obj.dp);
                case 'L1'
                    LSfun =  f +  phi*norm(c,1);
                    LSder = dfdw'*dw + dfdp'*obj.dp - phi*norm(c,1);
                case 'L2'
                    LSfun = f + phi*norm(c,2);
                    LSder = dfdw'*dw + dfdp'*obj.dp - phi*norm(c,2);
            end
        end
        
        function [] = lineBoxConstMinNAND(obj,p)
            %This function runs the BoxConstMin NAND optimization loop.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj - BoxConstMin object
            %p   - Initial point for optimization algorithm
            %
            %Outputs:
            %--------
            %No outputs.  Data stored in BoxConstMin object.
            %--------------------------------------------------------------
            
            %Set up the problem residual and jacobian with the appropriate
            %parameter
            p = p(:);
            obj.model.prob.updateParameters(p);
            
            n = length(p);
            if size(obj.paramBounds,1) == 1
                obj.paramBounds = repmat(obj.paramBounds,n,1);
            elseif size(obj.paramBounds,1) ~= n
                error('paramBounds does not have the appropriate number of rows');
            end
            
            tStartPDEopt = tic; %Timer
            %Run the model to solve the constraints and store the solution
            obj.model.executeModel;
            w = obj.model.sv(:,end);
            obj.solnHist = w;
            
            %Record the number of Solver calls and the parameter evolution
            obj.nSolverCalls = 1;
            obj.PDECallsParam = p; %Includes linesearch parameters
            obj.paramHist = p; %Does not include linesearch parameters
                        
            [Il,Iu,J,nact]=obj.determineActSet(p);
            
            %Compute the gradient of the objective
            DfDp = obj.model.reducedGrad(w,p,obj.sensMethod);
            %Initial Hessian approximation
            B = obj.reducedHessFD(p,DfDp,J,nact,n);

            cnt = 0;
            while 1
                while 1
                    if norm(DfDp(J),2) < obj.gradTol
                        break;
                    end
                    
                    cnt = cnt+1;
                    obj.dp = zeros(n,1);
                    obj.dp(J) = -B\DfDp(J); % QN updates on inverse Hessian
                    %obj.dp = -H*DfDp; % QN updates on inverse Hessian
                    
                    %Determine maximum step
                    alpha_max = obj.determineMaxStep(p,Il,Iu);
                    
                    switch lower(obj.details)
                        case 'newton'
                            alpha = alpha_max;
                            %Set up the problem residual and jacobian with the appropriate
                            %parameter
                            obj.model.prob.updateParameters(p);
                            %Determine which initial guess to use in the PDE
                            %solver
                            switch obj.PDEstartNew
                                case 'warm'
                                    obj.model.resetInitialCondition(obj.solnHist(:,end));
                                case 'reset'
                                    obj.model.resetInitialCondition(obj.model.prob.ic);
                            end
                            %Run the model to solve the constraints
                            obj.model.executeModel;
                            %Update parameters and increment solver
                            obj.PDECallsParam = [obj.PDECallsParam, p+alpha*obj.dp];
                            obj.nSolverCalls = obj.nSolverCalls + 1;
                        case 'wolfe'
                            %Linesearch
                            %[alpha,~,nFunCalls,alphaCalls] = linesearchBacktracking(@(alphaloc) obj.LineSearchReducedSpace(w,alphaloc),[obj.LSprop(1:2);alpha_max;obj.LSprop(4:end)]);
                            [alpha,~,nFunCalls,alphaCalls] = linesearch(@(alphaloc) obj.LineSearchReducedSpace(w,alphaloc),[obj.LSprop(1:2);alpha_max;obj.LSprop(4:end)]);
                            %[alpha,~,nFunCalls,alphaCalls] = linesearchOLD(@(alphaloc) obj.LineSearchReducedSpace(w,alphaloc),[obj.LSprop(1:2);alpha_max;obj.LSprop(4:end)]);
                            %Update parameters and increment solver
                            obj.PDECallsParam = [obj.PDECallsParam,repmat(p,1,size(alphaCalls,2))+repmat(obj.dp,1,size(alphaCalls,2)).*repmat(alphaCalls,size(obj.dp,1),1)];
                            obj.nSolverCalls = obj.nSolverCalls + nFunCalls;
                    end
                    
                    %Compute step, update iterate, and update parameter history
                    obj.dp = alpha*obj.dp;
                    p = p + obj.dp;
                    obj.model.prob.updateParameters(p);
                    obj.paramHist = [obj.paramHist,p];
                    
                    dpConv = abs(norm(obj.dp)/norm(p));
                    
                    if obj.printConv
                        fprintf('\n\n Iteration %d:\n',cnt);
                        disp(['Reduced Gradient Norm = ', num2str(norm(DfDp(J),2))]);
                        disp(['Step Length = ', num2str(alpha)]);
                        disp(['Max. Step Length = ', num2str(alpha_max)]);
                        disp(['||dp||/||p|| = ', num2str(dpConv)]);
                        f = obj.model.objective(w,p);
                        disp(['Objective Value = ',num2str(f)]);
                        disp(['Iterate : p = ',num2str(p')]);
                    end
                    
                    %If alpha_max is restricted by a constraint and we
                    %acheive a step of alpha_max, we must be on a new
                    %constraint
                    if alpha_max < 1 && (abs(alpha - alpha_max) < 1e-6)
                        [B,Il,Iu,J,nact] = obj.add2ActSet(B,p,J);
                    end
                    
                    %Determine whether to remove constraint
                    if sum(Il)+sum(Iu) > 0
                        [B,Il,Iu,J,nact] = obj.removeFromActSet(DfDp,Il,Iu,J,n,nact,B);
                    end                        
                    %If our step is very small, but we have not converged,
                    %re-compute the Hessian (FD)
                    
                    if (dpConv<obj.dpTol)
                        B = obj.reducedHessFD(p,DfDp,J,nact,n);
                    end

                    %We don't run the model at the new value of p because this
                    %has ALREADY BEEN DONE IN THE LINESEARCH.  It would be a
                    %waste to do it again.  Therefore, we simply update the
                    %values of paramHist and solnHist according to this LAST
                    %VALUE OF alpha from the linesearch.
                    w = obj.model.sv(:,end);
                    obj.solnHist = [obj.solnHist,w];
                    
                    %Record previous gradient and get the updated gradient
                    DfDp_old = DfDp;
                    DfDp = obj.model.reducedGrad(w,p,obj.sensMethod);
                    %DfDp = obj.reducedGrad(obj.solnHist(:,end),dcdw,dcdp);
                    %Use gradients and search directions to update the inverse
                    %Hessian
                    s = obj.dp;
                    y = DfDp - DfDp_old;
                    B = updateHessian(B,y(J),s(J),obj.secondOrder,1e-12);            
                end
                
                [conv,B,Il,Iu,J,nact] = checkConvOrRemoveConst(obj,DfDp,Il,Iu,J,n,nact,B);
                if conv
                    break;
                end
            end

            obj.time = toc(tStartPDEopt);
        end
        
        function  [] = lineSQPBoxConstMinFullSAND(obj,p)
            %This function runs the BoxConstMin full SAND optimization
            %loop.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj - BoxConstMin object
            %p   - Initial point for optimization algorithm
            %
            %Outputs:
            %--------
            %No outputs.  Data stored in BoxConstMin object.
            %--------------------------------------------------------------
            
            tStartPDEopt = tic; %Timer
            
            %Set up the problem residual and jacobian with the appropriate
            %parameter
            p = p(:);
            obj.model.prob.updateParameters(p);
            %Run the model to solve the constraints and store the solution
            obj.model.executeModel;
            obj.solnHist = obj.model.sv(:,end);
            w = obj.solnHist;
            
            %Extract model information
            [~,dfdw,dfdp] = obj.model.objective(w,p);
            [R,dRdw,dRdp] = obj.model.PDEconstraint(w,p);
            
            %Form the full gradient of f
            DLDx = [dfdw; dfdp];
            DfDx = DLDx;
            DfDx0 = DfDx;
            
            %Tolerances
            eps = 1e-4;
            eps_DfDp = 0.15;
            
            %Initial Hessian guess (SPD approx)
            H = speye(size(dRdw,1)+size(p,1));
            H(end - (size(p,1)-1),end - (size(p,1)-1)) = 0;
            
            %Store iterate and initial Solver counter and PDEsolver counter
            obj.paramHist = p;
            obj.nSolverCalls = 1;
            obj.PDECallsParam = p;
            cnt = 0;
            alpha = 0;
            alpha_max = 0;
            phi0 = 0;
            obj.dp = 0;
            lambda = 0;
            
            while 1
                cnt = cnt + 1;
                
                normR = norm(R);
                normDfDx = norm(DfDx);
                dpConv = abs(norm(obj.dp)/norm(p));
                if obj.printConv
                    fprintf('\n\n Iteration %d:\n',cnt);
                    disp(['Gradient Norm = ', num2str(normDfDx)]);
                    disp(['Step Length = ', num2str(alpha)]);
                    disp(['Max. Step Length = ', num2str(alpha_max)]);
                    disp(['||dp||/||p|| = ', num2str(dpConv)]);
                    disp(['Residual 2-norm = ',num2str(normR)]);
                    f = obj.model.objective(w,p);
                    disp(['Objective Value = ',num2str(f)]);
                    disp(['norm(lambda) = ',num2str(norm(lambda))]);
                end
                
                %Check for convergence
                if (normDfDx<eps*norm(DfDx0) && normR<eps)
                    disp('Converged!');
                    break;
                end
                
                %Setup and solve KKT system
                dR = [dRdw, dRdp];
                K   = [H , dR';
                    dR, zeros(size(dRdw,1))];
                rhs = [-dfdw; -dfdp; -R];
                sol = K\rhs;
                dw = sol(1:size(dRdw,1));
                obj.dp = sol(1+size(dRdw,1):size(dRdw,1)+size(dRdp,2));
                lambda = -sol(1+size(dRdw,1)+size(dRdp,2):end);
                obj.dps = [obj.dps, obj.dp];
                
                %Determine alpha_max0 based on size of gradient
                if (normDfDx<eps_DfDp)
                    alpha_max0 = 1;
                else
                    alpha_max0  = 0.5;
                end
                
                %This selection of alpha_max ensures that
                %p_low <= p + alpha_max*dp <= p_up.
                %Rearranging this equation, we get that alpha_max must
                %satsify the following two conditions:
                %alpha_max <= ((p_up - p)./dp)_j   if dp_j > 0
                %alpha_max <= ((p - p_low)./dp)_j  if dp_j < 0
                if isempty(obj.paramBounds)
                    alpha_max = min(obj.LSprop(3),1);
                else
                    alpha_max = min(...
                        (obj.dp > 0).*(obj.paramBounds(:,2) - p)./obj.dp + ...
                        (obj.dp < 0).*(obj.paramBounds(:,1) - p)./obj.dp);
%                     alpha_max = min([alpha_max,obj.LSprop(3),1]);
                end
                
                alpha_max = min(alpha_max,alpha_max0);
                dlambda = 0;
                
                %Linesearch
                LSiter = 0;
                LSiterMax = 5;
                LSsuccess = 0;
                while (~LSsuccess)
                    if(LSiter == 0)
                        if(phi0 < 1)
                            phi = 0;
                        else
                            phi = phi0;
                        end
                    else
                        if(phi0 == 0)
                            phi = 1*10^(LSiter-1);
                        else
                            phi = phi0*10^(LSiter);
                        end
                    end
                    disp(['current phi is ', num2str(phi)]);
                    LineSearchFunQN = @(alpha_loc) obj.LineSearchFullSpace(p,w,lambda,dw,dlambda,alpha_loc,phi);
                    
                    %[alpha,LSsuccess] = linesearch(LineSearchFunQN,[obj.LSprop(1:2);alpha_max;obj.LSprop(4:end)]);
                    %[alpha,LSsuccess] = linesearch(LineSearchFunQN,[1e-4;0.95;alpha_max;5;10;0.05;0.1;1.2]);
                    %[alpha,LSsuccess] = linesearch(LineSearchFunQN,[1e-7;0.99;alpha_max;5;10;0.05;0.1;1.2]);
                    [alpha,LSsuccess] = linesearchBacktracking(LineSearchFunQN,[1e-7;0.99;alpha_max;50;10;0.05;0.1;1.2]);
                    
                    if(LSsuccess)
                        if(LSiter == 0)
                            phi0 = phi/10;
                        else
                            phi0 = phi;
                        end
                        break;
                    end
                    
                    LSiter = LSiter + 1;
                    if (LSiter > LSiterMax && ~LSsuccess)
                        phi0 = phi;
                        break;
                    end
                end
                
                %Update variables
                obj.dp = alpha*obj.dp;
                dw = alpha*dw;
                
                w = w + dw;
                p = p + obj.dp;
                obj.solnHist = [obj.solnHist, w];
                obj.paramHist = [obj.paramHist, p];
                
                %Compute terms for next iteration
                obj.model.prob.updateParameters(p);
                [R,dRdw,dRdp] = obj.model.prob.ResSens(w,p);
                [~,dfdw,dfdp] = obj.model.objective(w,p);
                DLDx_old = DLDx;
                DfDx = [dfdw; dfdp];
                DLDx = [dRdw'; dRdp']*lambda;

%                 if mod(cnt,20) == 0 || condest(H) > 1e8
%                     H = speye(size(dRdw,1)+size(p,1));
%                 else
%                     s = [dw; obj.dp];
%                     y = DLDx - DLDx_old;
%                     H = dampedBFGS(H,y,s,obj.secondOrder);
%                 end
                
%                 s = [dw; obj.dp];
%                 y = DLDx - DLDx_old;
%                 H = updateHessian(H,y,s,obj.secondOrder,1e-5);
                
%                 s = dw; y = DLDx(1:length(R),1) - DLDx_old(1:length(R),1);
%                 H(1:length(R),1:length(R)) = updateHessian(H(1:length(R),1:length(R)), y, s, obj.secondOrder,1e-6);
                
                
%                 s = [dw; zeros(size(obj.dp))]; 
% %                 s = [dw; obj.dp];
%                 y = DLDx - DLDx_old; y(end-length(obj.dp)+1:end,1) = 0;
%                 H = updateInverseHessian(H,y,s,obj.secondOrder,1e-3);
                
                disp(['Iterate : p = ',num2str(p')]);
            end
            obj.time = toc(tStartPDEopt);
        end
        
        function  [] = lineSQPBoxConstMinRedSAND(obj,p)
            %This function runs the BoxConstMin reduced SAND optimization
            %loop.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj - BoxConstMin object
            %p   - Initial point for optimization algorithm
            %
            %Outputs:
            %--------
            %No outputs.  Data stored in BoxConstMin object.
            %--------------------------------------------------------------
            
            tStartPDEopt = tic; %Timer
            
            %Set up the problem residual and jacobian with the appropriate
            %parameter
            p = p(:);
            obj.model.prob.updateParameters(p);
            %Run the model to solve the constraints and store the solution
            obj.model.executeModel;
            obj.solnHist = obj.model.sv(:,end);
            w = obj.solnHist;
            
            %Evalutate the objective and its gradient
            [DfDp,f,lambda,R,dRdw,dRdp] = obj.model.reducedGrad(w,p,obj.sensMethod);
            %Lagrange multipliers
            dlambda = lambda;
            lambda = zeros(size(lambda));
            
            %Store iterate and initial Solver counter and PDEsolver counter
            obj.paramHist = p;
            obj.nSolverCalls = 1;
            obj.PDECallsParam = p;
            cnt = 0;
            alpha = 0;
            alpha_max = 0;
            phi0 = 0;
            obj.dp = 0;
            normR = norm(R);
            
            %Tolerances
            eps = 1e-5;
            eps_dp = 1e-6;
            eps_const = 5;
%             eps_const = 100;
            eps_DfDp = 0.15;
            eps_dlambda = 1000000000.01;
            dlambda_max = 1;
            
            %Initial Hessian guess (SPD approx)
            H = speye(size(DfDp,1));
            
            while 1
                cnt = cnt + 1;
                
                %Print information to user
                normDfDp = norm(DfDp);
                dpConv = norm(obj.dp)/norm(p);
                if obj.printConv
                    fprintf('\n\n Iteration %d:\n',cnt);
                    disp(['Gradient Norm = ', num2str(normDfDp)]);
                    disp(['Step Length = ', num2str(alpha)]);
                    disp(['Max. Step Length = ', num2str(alpha_max)]);
                    disp(['||dp||/||p|| = ', num2str(dpConv)]);
                    disp(['DfDp = ', num2str(normDfDp)]);
                    disp(['Residual 2-norm = ',num2str(normR)]);
                    disp(['Iterate  = ',num2str(p')]);
                end
                
                %Check for convergence
                if (normDfDp<eps && normR <eps)
                    disp('Converged!');
                    break;
                end

                % solve for decision variable
                obj.dp = -H*DfDp;
                
                % solve for state variable
                switch class(obj.model)
                    case 'FOM'%{'HFM','GALERKIN_ROM','MINRES_ROM_FD'}
                        dw = -dRdw\(dRdp*obj.dp+R);
                    case 'ROM'%'MINRES_ROM'
                        dw = -dRdw'\(dRdp*obj.dp+R);
                        dw = dRdw\dw;
                end
                
                %Determine alpha_max0 based on size of gradient
                if (normDfDp<eps_DfDp)
                    alpha_max0 = 1;
                else
                    alpha_max0  = 0.5;
                end
                
%                 %Make sure you don't overstep the boundaries
%                 if (obj.dp>0)
%                     alpha_max = min(min((obj.paramBounds(2)-p)/obj.dp),alpha_max0);
%                 elseif (obj.dp==0)
%                     break;
%                 else
%                     alpha_max = min(min((obj.paramBounds(1)-p)/obj.dp),alpha_max0);
%                 end
%                 
                %This selection of alpha_max ensures that
                %p_low <= p + alpha_max*dp <= p_up.
                %Rearranging this equation, we get that alpha_max must
                %satsify the following two conditions:
                %alpha_max <= ((p_up - p)./dp)_j   if dp_j > 0
                %alpha_max <= ((p - p_low)./dp)_j  if dp_j < 0
                if isempty(obj.paramBounds)
                    alpha_max = min(obj.LSprop(3),1);
                else
                    alpha_max = min(...
                        (obj.dp > 0).*(obj.paramBounds(:,2) - p)./obj.dp + ...
                        (obj.dp < 0).*(obj.paramBounds(:,1) - p)./obj.dp);
                    alpha_max = min([alpha_max,obj.LSprop(3),1]);
                end
                
                
                alpha_max = min(alpha_max,alpha_max0);
                
                %Linesearch
                LSiter = 0;
                LSiterMax = 5;
                LSsuccess = 0;
                while (~LSsuccess)
                    if(LSiter == 0)
                        if(phi0 < 1)
                            phi = 0;
                        else
                            phi = phi0;
                        end
                    else
                        if(phi0 == 0)
                            phi = 1*10^(LSiter-1);
                        else
                            phi = phi0*10^(LSiter);
                        end
                    end
                    disp(['current phi is ', num2str(phi)]);
                    LineSearchFunQN = @(alpha_loc) obj.LineSearchFullSpace(p,w,lambda,dw,dlambda,alpha_loc,phi);
                    %[alpha,LSsuccess] = linesearch(LineSearchFunQN,[1e-4;0.9;alpha_max;5;10;0.05;0.1;1.2]);
                    %[alpha,LSsuccess] = linesearch2(LineSearchFunQN,[1e-4;0.9;alpha_max;5;10;0.05;0.1;1.2]);
                    [alpha,LSsuccess] = linesearchBacktracking(LineSearchFunQN,[10^(-7);0.99;alpha_max;50;10;0.05;0.1;1.2]);
                    
                    if(LSsuccess)
                        if(LSiter == 0)
                            phi0 = phi/10;
                        else
                            phi0 = phi;
                        end
                        break;
                    end
                    
                    LSiter = LSiter + 1;
                    if (LSiter > LSiterMax && ~LSsuccess)
                        phi0 = phi;
                        break;
                    end
                end
                
                %Update variables
                obj.dp = alpha*obj.dp;
                dw = alpha*dw;
                dlambda = alpha*dlambda;
                
                w = w + dw;
                p = p + obj.dp;
                lambda = lambda + dlambda;
                obj.solnHist = [obj.solnHist, w];
                obj.paramHist = [obj.paramHist, p];

                %Store previous gradient
                DfDp_old = DfDp;
  
                %Evalutate the objective and its gradient
                obj.model.prob.updateParameters(p);
                [DfDp,f,lambda_new,R,dRdw,dRdp] = obj.model.reducedGrad(w,p,obj.sensMethod);
                normR = norm(R,2);
                %Determine Lagrange multiplier with largest magnitude
                max_new_lambda = max(abs(lambda_new));
                disp(['maximum absolute value for new lambda is ', num2str(max_new_lambda)] );
                
                if (normR > eps_const)
                    obj.nSolverCalls = obj.nSolverCalls + 1;
                    obj.PDECallsParam = [obj.PDECallsParam,p];
                    
                    %Run the model to solve the constraints and store the solution
                    obj.model.prob.updateParameters(p);
                    obj.model.executeModel;
                    w = obj.model.sv(:,end);
                    obj.solnHist = [obj.solnHist,w];
                    
                    %Evalutate the objective and its gradient
                    [DfDp,f,lambda_new,R,dRdw,dRdp] = obj.model.reducedGrad(w,p,obj.sensMethod);
                    normR = norm(R,2);
                    %Determine Lagrange multiplier with largest magnitude
                    max_new_lambda = max(abs(lambda_new));
                    disp(['constraint violation is ', num2str(normR)] );
                    disp(['maximum absolute value for new lambda after PDE solver call is ', num2str(max_new_lambda)] );
                end
                switch class(obj.model)
                    case 'FOM'%{'HFM','GALERKIN_ROM'}
                        dlambda = lambda_new-lambda;
                    case 'ROM'%{'MINRES_ROM','MINRES_ROM_FD'}
                        dlambda = lambda_new-JV\(JVold*lambda);
                        JVold = JV;
                end
                s = obj.dp;
                y = DfDp - DfDp_old;
                H = updateInverseHessian(H,y,s,obj.secondOrder,1e-3);
            end
            obj.time = toc(tStartPDEopt);
        end
        
        function  [x,fval,exitflag,output,lambda] = builtinFMINCON(obj,p,w0,lowlim,uplim,options)
            
            if nargin < 6 || isempty(options)
                options = optimset('Algorithm','interior-point',...
                    'Display','iter-detailed',...
                    'GradConstr','on',...
                    'GradObj','on',...
                    'TolFun',1e-6,...
                    'TolCon',1e-6,...
                    'TolX',1e-8);
            end
                           
            func = @(u) obj.model.ObjFMINCON_SAND(u);
            nlconst = @(u) obj.model.NLConstFMINCON_SAND(u);
            
            if strcmpi(class(obj.model),'GNAT')
                Aineq=[]; bineq=[];
%                 [Aineq,bineq] = obj.model.probGNAT.linIneqGNAT();
%                 bineq = bineq - A*obj.ic;
%                 Aineq = Aineq(unique(obj.model.jrow),:);
%                 Aineq = Aineq*obj.model.phi;
%                 Aineq = [Aineq,zeros(obj.model.ndof-1,length(p))];
            elseif strcmpi(class(obj.model),'ROM')
                [Aineq,bineq] = obj.model.prob.linIneqROM();
                Aineq = Aineq*obj.model.phi;
                Aineq = [Aineq,zeros(obj.model.ndof-1,length(p))];
            elseif strcmpi(class(obj.model),'FOM')
                [Aineq,bineq] = obj.model.prob.linIneqROM();
                Aineq = [Aineq,zeros(obj.model.ndof-1,length(p))];
            end
                            
            switch class(obj.model)
                case 'FOM'
                    if nargin == 3
                        ic = w0;
                    else
                        ic = obj.model.prob.ic;
                    end
                    lb =  -inf(obj.model.ndof,1);
                    ub = inf(obj.model.ndof,1);
                case 'ROM'
                    if nargin == 3
                        ic = obj.model.phi'*(w0 - obj.model.prob.ic);
                    else
                        ic = zeros(obj.model.nY,1);
                    end
                    lb =  -inf(obj.model.nY,1);
                    ub = inf(obj.model.nY,1);
                case 'GNAT'
                    if nargin == 3
                        ic = obj.model.phiYhat'*(w0(unique(obj.model.jrow)) - obj.model.probGNAT.ic);
                    else
                        ic = zeros(obj.model.nY,1);
                    end
                    lb =  -inf(obj.model.nY,1);
                    ub = inf(obj.model.nY,1);
            end
            
            tStart = tic;
            [x,fval,exitflag,output,lambda]=fmincon(func,[ic;p],Aineq,bineq,[],[],[lb;lowlim],[ub;uplim],nlconst,options);
            obj.time = toc(tStart);
        end
        
        function  [x,fval,exitflag,output,lambda] = fmincon_sand(obj,p0,w0,options,A,b,nlconstraints)  
            
            %Get initial number of model evaluations
            InitNumModelEval=obj.model.numExecute;
            
            if nargin < 4 || isempty(options)
                options = optimset('Algorithm','sqp',...
                    'Display','iter-detailed',...
                    'GradConstr','on',...
                    'GradObj','on',...
                    'Hessian','bfgs',...
                    'TolFun',1e-6,...
                    'TolCon',1e-6,...
                    'TolX',1e-8);
            end
            if nargin < 5, A=[]; end;
            if nargin < 6, b=[]; end;
            if nargin < 7, nlconstraints=[]; end;
            
            func = @(u) obj.model.ObjFMINCON_SAND(u);
%             if iscell(nlconstraints)
%                 nlconst = @(u) obj.model.NLConstFMINCON_SAND(u,nlconstraints{:});
%             else
%                 nlconst = @(u) obj.model.NLConstFMINCON_SAND(u,nlconstraints);
%             end
 nlconst = @(u) obj.model.NLConstFMINCON_SAND(u);
%             if iscell(nlconstraints)
%                 nlconst = @(u) obj.concatNLConstraints(u,@(x)obj.model.NLConstFMINCON_SAND(x),nlconstraints{:});
%             else
%                 nlconst = @(u) obj.concatNLConstraints(u,@(x)obj.model.NLConstFMINCON_SAND(x),nlconstraints);
%             end
            
            switch class(obj.model)
                case 'FOM'
                    ic = [w0;p0];
                    lowlim = [-inf(size(w0));obj.paramBounds(:,1)];
                    uplim  = [inf(size(w0));obj.paramBounds(:,2)];
                case 'ROM'
                    ic = [w0;p0];
                    lowlim = [-inf(size(w0));obj.paramBounds(:,1)];
                    uplim  = [inf(size(w0));obj.paramBounds(:,2)];
                case 'GNAT'
                    ic = [w0;p0];
                    lowlim = [-inf(size(w0));obj.paramBounds(:,1)];
                    uplim  = [inf(size(w0));obj.paramBounds(:,2)];
            end
            
            tStart = tic;
            [x,fval,exitflag,output,lambda]=fmincon(func,ic,A,b,[],[],lowlim,uplim,nlconst,options);
            obj.time = toc(tStart);
            
            %Determine total number of model evaluations
            obj.numModelEval=obj.model.numExecute - InitNumModelEval;
        end

        function  [x,fval,exitflag,output,lambda] = fmincon_sand_romvar(obj,p0,w0,options,A,b,nlconstraints)  
            
            if nargin < 4 || isempty(options)
                options = optimset('Algorithm','sqp',...
                    'Display','iter-detailed',...
                    'GradConstr','on',...
                    'GradObj','on',...
                    'Hessian','bfgs',...
                    'TolFun',1e-6,...
                    'TolCon',1e-6,...
                    'TolX',1e-8);
            end
            if nargin < 5, A=[]; end;
            if nargin < 6, b=[]; end;
            if nargin < 7, nlconstraints=[]; end;
            
            func = @(u) obj.model.ObjFMINCON_SAND_ROMVAR(u);
%             if iscell(nlconstraints)
%                 nlconst = @(u) obj.model.NLConstFMINCON_SAND(u,nlconstraints{:});
%             else
%                 nlconst = @(u) obj.model.NLConstFMINCON_SAND(u,nlconstraints);
%             end
 nlconst = @(u) obj.model.NLConstFMINCON_SAND(u);
%             if iscell(nlconstraints)
%                 nlconst = @(u) obj.concatNLConstraints(u,@(x)obj.model.NLConstFMINCON_SAND(x),nlconstraints{:});
%             else
%                 nlconst = @(u) obj.concatNLConstraints(u,@(x)obj.model.NLConstFMINCON_SAND(x),nlconstraints);
%             end
            
            switch class(obj.model)
                case 'FOM'
                    ic = [w0;p0];
                    lowlim = [-inf(size(w0));obj.paramBounds(:,1)];
                    uplim  = [inf(size(w0));obj.paramBounds(:,2)];
                case 'ROM'
                    ic = [w0;p0];
                    lowlim = [-inf(size(w0));obj.paramBounds(:,1)];
                    uplim  = [inf(size(w0));obj.paramBounds(:,2)];
                case 'GNAT'
                    ic = [w0;p0];
                    lowlim = [-inf(size(w0));obj.paramBounds(:,1)];
                    uplim  = [inf(size(w0));obj.paramBounds(:,2)];
            end
            
            tStart = tic;
            [x,fval,exitflag,output,lambda]=fmincon(func,ic,A,b,[],[],lowlim,uplim,nlconst,options);
            obj.time = toc(tStart);
        end

        function  [x,fval,exitflag,output,lambda] = fmincon_nand(obj,p0,options,A,b,Aeq,beq,nlconstraints,objOverwrite)
                     
            %Get initial number of model evaluations
            InitNumModelEval=obj.model.numExecute;
            
            if  isempty(options)
                options = optimset('Algorithm','sqp',...
                    'Display','iter-detailed',...
                    'GradObj','on',...
                    'GradConstr','on',...
                    'Hessian','bfgs',...
                    'TolFun',1e-6,...
                    'TolCon',1e-12,...
                    'TolX',1e-6);
            end
            
            if strcmpi(obj.sensMethod,'findiff') && (isempty(options.GradObj) || strcmpi(options.GradObj,'on'))
                warning('Sensitivity Method is finite differences, so changing fmincon option ''GradObj'' to ''off''');
                options.GradObj = 'off';
            end
            
            if nargin < 4, A=[]; end;
            if nargin < 5, b=[]; end;
            if nargin < 6, Aeq=[]; end;
            if nargin < 7, beq=[]; end;
            if nargin < 8, nlconstraints=[]; end;
            if nargin < 9 || isempty(objOverwrite)
                func = @(u) obj.model.ObjFMINCON_NAND(u,obj.sensMethod);
            
                %Run model at initial condition
                switch class(obj.model)
                    case 'GNAT'
                        obj.model.probGNAT.updateParametersGNAT(p0);
                    otherwise
                        if ~strcmpi(class(obj.model.prob),'structuralFEM')
                            obj.model.prob.updateParameters(p0);
                        else
                            n1 = obj.model.prob.PrecompROM.n_lam; n2=obj.model.prob.PrecompROM.n_mu; n3=obj.model.prob.PrecompROM.n_rho;

                            obj.model.prob.PrecompFromRedMatPrecomp(p0(1:n1),p0(n1+1:n1+n2),p0(n1+n2+1:end));
                        end
                end
                obj.model.executeModel;
            else
                func = @(u) obj.predefObj(u,objOverwrite);
            end
            
            %Constraints
            if ~iscell(nlconstraints),  nlconstraints={nlconstraints}; end;
            nlconst = @(u) obj.concatNLConstraints(u,nlconstraints{:});
                
            tStart = tic;
            [x,fval,exitflag,output,lambda]=fmincon(func,p0,A,b,Aeq,beq,obj.paramBounds(:,1),obj.paramBounds(:,2),nlconst,options);
            obj.time = toc(tStart);
            
            %Determine total number of model evaluations
            obj.numModelEval=obj.model.numExecute - InitNumModelEval;
        end
        
        function  [x,info] = ipopt_nand(obj,p0,tol,nlconst,outputfunc,ipoptOptions)
               
            %Get initial number of model evaluations
            InitNumModelEval=obj.model.numExecute;
            
            % Set the IPOPT options.
%             options.ipopt.accept_after_max_steps = 10;
%             options.ipopt.tiny_step_tol    = 1e-2;
            if nargin < 6 || isempty(ipoptOptions)
                options.ipopt.mu_strategy      = 'adaptive';
                options.ipopt.max_iter         = 5;
                options.ipopt.tol              = 1e-6;
                options.ipopt.hessian_approximation = 'limited-memory';
            else
                options.ipopt=ipoptOptions;
            end
            
            % Set bounds on parameters and constraints
            options.lb = obj.paramBounds(:,1)';
            options.ub = obj.paramBounds(:,2)';

            [nc,nceq] = nlconst(p0,'sizes');
            options.cl = [-inf(nc,1);zeros(nceq,1)]';
            options.cu = zeros(nc+nceq,1)';
            
            % The callback functions.
            funcs.objective         = @(p) obj.model.ObjFMINCON_NAND(p,obj.sensMethod,'objective');
            funcs.gradient          = @(p) obj.model.ObjFMINCON_NAND(p,obj.sensMethod,'gradient');
            funcs.constraints       = @(p) nlconst(p,'constraint');
            funcs.jacobian          = @(p) nlconst(p,'jacobian');
            funcs.jacobianstructure = @() nlconst(p0,'jacobianstructure');
            
            %Run model at initial condition
            switch class(obj.model)
                case 'GNAT'
                    obj.model.probGNAT.updateParametersGNAT(p0);
                otherwise
                    obj.model.prob.updateParameters(p0);
            end
            obj.model.executeModel;
            
            % Run IPOPT.
            [x,info] = ipopt(p0,funcs,options);
            
            %Determine total number of model evaluations
            obj.numModelEval=obj.model.numExecute - InitNumModelEval;
        end

        function  [x_k, xstate, xmul, F, Fmul, Fstate, Inform, nS, nInf, sInf] = snopt_nand(obj,x_0,tol,nlconst,outputfunc,snoptOptions)

            %Get initial number of model evaluations
            InitNumModelEval=obj.model.numExecute;
            
            % Set bounds on parameters and constraints
            x_L = obj.paramBounds(:,1);
            x_U = obj.paramBounds(:,2);

            % The callback functions.
            funcs.objective         = @(p) obj.model.ObjFMINCON_NAND(p,obj.sensMethod,'objective');
            funcs.gradient          = @(p) obj.model.ObjFMINCON_NAND(p,obj.sensMethod,'gradient');
            funcs.constraints       = @(p) nlconst(p,'constraint');
            funcs.jacobian          = @(p) nlconst(p,'jacobian');

            Func = @(x,status,Prob) obj.combineObjConstr(funcs,x,status,Prob);           
            objRow=1;

            % optPar similar to snopt except: optPar(14:17) are not available in
            % snoptA. Refer to snoptTL.m 
            optPar = -999*ones(72,1);
            optPar(1) = 1;
            optPar(2) = 1;
            optPar(9) = 1e-4; %Major iteration tolerance
            optPar(13) = -1; %Derivative checking
%             optPar(22) = 0.5;
            optPar(35) = 100; %Maximum number of major iterations
            optPar(40) = 0; %Linsearch (0=nonderivative linesearch)
            %optPar(41) = 1e-3;
            optPar(63)=1;
%             optPar(41) = 1;

            [nc,nceq] = nlconst(x_0,'sizes');
            b_L = [-inf;-inf(nc,1);zeros(nceq,1)];
            b_U = [inf;zeros(nc+nceq,1)];

            A = sparse(zeros(nc+nceq+1,length(x_0)));
            G = sparse([ones(1,length(x_0),1); nlconst(x_0,'jacobianstructure')]);
            
            Warmstart=0;
            xState=[];
            fState=[];
            nS=[];
 
            Prob = struct('P',1); %Not used by snoptA mex, can be anything
            PrintFile='snprint.txt';
            SummFile='snsum.txt';
            PriLev=1;

            obj.model.prob.updateParameters(x_0);
            obj.model.executeModel();
            [x_k, xstate, xmul, F, Fmul, Fstate, Inform, nS, nInf, sInf] = ...
                    snoptA( x_L, x_U, x_0, b_L, b_U, A, G, Func, objRow, optPar, ...
                          Warmstart, xState, fState, nS, Prob, ...
                               PrintFile, SummFile, PriLev );
                           
            %Determine total number of model evaluations
            obj.numModelEval=obj.model.numExecute - InitNumModelEval;
        end

        function [f,G,mode] = combineObjConstr(obj,funcs,x,status,Prob)

           mode  = status(1); %1 = first call
           needf = status(2);
           needG = status(3);
             
           f=[]; G=[];
           if needf>0
              f1=funcs.objective(x);
              f2=funcs.constraints(x);
              f=[f1;f2];
              if sum(isnan(f))>0 || ~isreal(f)
                 mode=-1;
              end
           end

           if needG>0
              G1=funcs.gradient(x);
              G2=funcs.jacobian(x);
              G=[G1';G2];
              if sum(sum(isnan(G)))>0 || ~isreal(G)
                 mode=-1;
              end
           end
        end

        function  [c,dc,ceq,dceq] = concatNLConstraints(obj,u,varargin)
            
            c=[]; dc=[]; ceq=[]; dceq=[];
            for i = 1:length(varargin)
                if isempty(varargin{i}), continue; end;
                
                [tmp,dtmp,tmpeq,dtmpeq] = varargin{i}(u);
                
                c=[c;tmp];
                ceq=[ceq;tmpeq];
                
                dc=[dc,dtmp];
                dceq=[dceq,dtmpeq];
            end
        end
        
        function  [f,df] = predefObj(obj,p,objOverwrite)
            
            if strcmpi(objOverwrite.objective,'feasibility')
                f  = 0;
                df = zeros(length(p),1);
            elseif strcmpi(objOverwrite.objective,'mindist')
                f  = 0.5*norm(p-objOverwrite.point)^2;
                df = p-objOverwrite.point;
            end
        end
        
        function  [Il,Iu,J,nact] = determineActSet(obj,p)
            Il   = (abs(p - obj.paramBounds(:,1)) < 1e-8);
            Iu   = (abs(p - obj.paramBounds(:,2)) < 1e-8);
            J    = ~(Il|Iu);
            nact = sum(Il)+sum(Iu);
        end
        
        function  [alpha_max] = determineMaxStep(obj,p,Il,Iu)
            %This selection of alpha_max ensures that
            %p_low <= p + alpha_max*dp <= p_up.
            %Rearranging this equation, we get that alpha_max must
            %satsify the following two conditions:
            %alpha_max <= ((p_up - p)./dp)_j   if dp_j > 0
            %alpha_max <= ((p - p_low)./dp)_j  if dp_j < 0
            
            if isempty(obj.paramBounds)
                alpha_max = min(obj.LSprop(3),1);
            else
                
                gammaL = -(p(~Il) - obj.paramBounds(~Il,1))./obj.dp(~Il);
                gammaL(obj.dp(~Il) >= 0) = inf;
                
                gammaU = (obj.paramBounds(~Iu,2)-p(~Iu))./obj.dp(~Iu);
                gammaU(obj.dp(~Iu) <= 0) = inf;
                
                alpha_max = min([gammaU;gammaL;1]);
                %alpha_max = min(...
                %    (obj.dp > 0).*(obj.paramBounds(:,2) - p)./obj.dp + ...
                %    (obj.dp < 0).*(obj.paramBounds(:,1) - p)./obj.dp);
                %alpha_max = min([alpha_max,obj.LSprop(3),1]);
            end
        end
        
        function  [B,Il,Iu,J,nact] = removeFromActSet(obj,DfDp,Il,Iu,J,n,nact,B)
            
            lamL = inf(n,1);
            lamU = inf(n,1);
            
            lamL(Il) = DfDp(Il);
            lamU(Iu) = -DfDp(Iu);
            
            [sigmaL,indL] = min(lamL);
            [sigmaU,indU] = min(lamU);
            [sig,tmp] = min([sigmaL(1);sigmaU(1)]);
            ind = zeros(n,1);
            if tmp==1, ind(indL(1))=true; else ind(indU(1))=true; end;
            
            beta=0.5;
            if sig < 0 && norm(DfDp(J)) < beta*norm(DfDp(J|ind))
                if sigmaL(1) < sigmaU(1)
                    Il(indL(1)) = false;
                else
                    Iu(indU(1)) = false;
                end
                newJ = ~(Il|Iu);
                ind = (J(newJ) ~= (~(Il(newJ) | Iu(newJ))));
                
                J = newJ;
                nact = nact - 1;
                tmp=B;
                B = zeros(sum(J));
                B(~ind,~ind) = tmp;
                B(ind,ind) = 1;
            end
        end
        
        function  [conv,B,Il,Iu,J,nact] = checkConvOrRemoveConst(obj,DfDp,Il,Iu,J,n,nact,B)
            lamL = inf(n,1);
            lamU = inf(n,1);
            
            lamL(Il) = DfDp(Il);
            lamU(Iu) = -DfDp(Iu);
            
            [sigmaL,indL] = min(lamL);
            [sigmaU,indU] = min(lamU);
            sig = min([sigmaL(1);sigmaU(1)]);
            
            if sig > 0
                conv=true;
            else
                conv=false;
                if sigmaL(1) < sigmaU(1)
                    Il(indL(1)) = false;
                else
                    Iu(indU(1)) = false;
                end
            end
            newJ = ~(Il|Iu);
            ind = (J(newJ) ~= (~(Il(newJ) | Iu(newJ))));
            
            J = newJ;
            nact = nact - 1;
            tmp=B;
            B = zeros(sum(J));
            B(~ind,~ind) = tmp;
            B(ind,ind) = 1;
        end
        
        function  [B,Il,Iu,J,nact] = add2ActSet(obj,B,p,J)
            [Il,Iu]=obj.determineActSet(p);
            ind = (Il(J) == J(J)) | (Iu(J) == J(J));
            [~,~,J,nact]=obj.determineActSet(p);
            B = B(~ind,~ind);
        end
        
        function  [B] = reducedHessFD(obj,p,DfDp,J,nact,n)
            h=1e-6;
            B = zeros(n-nact);
            ind = find(J);
            for i = 1:sum(J)
                padj = p; padj(ind(i)) = padj(ind(i))+h;
                obj.model.prob.updateParameters(padj);
                obj.model.executeModel;
                wadj = obj.model.sv(:,end);
                
                tmp = (1/h)*(obj.model.reducedGrad(wadj,padj,obj.sensMethod) - DfDp);
                B(:,i) = tmp(J);
            end
            B = 0.5*(B+B');
            R = modifiedCholesky(B);
            B = R'*R;
            obj.model.prob.updateParameters(p);
        end
        
        function  [] = executeOptimization(obj,p)
            %This function calls the appropriate BoxConstMin
            %PDE-constrained optimization algorithm
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj - BoxConstMin object
            %
            %Outputs:
            %--------
            %No outputs.  Data stored in BoxConstMin object.
            %--------------------------------------------------------------
            
            switch obj.optFrame
                case 'NAND'
                    %NAND optimization framework
                    if strcmpi(obj.optStrategy,'LS')
                        %Line search strategy
                        lineBoxConstMinNAND(obj,p);
                    elseif strcmpi(obj.optStrategy,'TR')
                        %Trust region strategy
                        trustBoxConstMinNAND(obj,p);
                    end
                case 'fullSAND'
                    %Full SAND optimization framework
                    if strcmpi(obj.optStrategy,'LS')
                        %Line search strategy
                        lineSQPBoxConstMinFullSAND(obj,p);
                    elseif strcmpi(obj.optStrategy,'TR')
                        %Trust region strategy
                        trustBoxConstMinSAND(obj,p);
                    end
                case 'reducedSAND'
                    %Reduced SAND optimization framework
                    if strcmpi(obj.optStrategy,'LS')
                        %Line search strategy
                        lineSQPBoxConstMinRedSAND(obj,p);
                    elseif strcmpi(obj.optStrategy,'TR')
                        %Trust region strategy
                        trustBoxConstRedSAND(obj,p);
                    end
            end
        end
    end
end

classdef QP < handle
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
        
        PDEstartLS  = 'reset';
        PDEstartNew = 'warm';
        
        G;
        c;
        A;
        I; %indicies of inequality constrains
        m; %total number of constraints
        b;
        
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
        function [obj] = QP(G,c,A,b,I)
            %This is the constructor of the QuasiNewton class.  It reads the
            %appropriate entries from the input file and stores them in the
            %class instance.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %G - matrix defining quadratic term in QP subproblem
            %c - vector defining linear term in QP subproblem
            %A - matrix defining all constraints (linear and nonlinear)
            %
            %Outputs:
            %--------
            %obj     - instance of the OPT object that was constructed
            %--------------------------------------------------------------
            
            obj.G = G;
            obj.c = c;
            obj.A = A;
            obj.b = b;
            obj.I = I;
            
            LHS = [G, -A';A zeros(size(A,1))];
            RHS = [-c;zeros(size(A,1),1)];
            %Set model object as a property
            %obj.model = modelobj;
            %Extract the OPT text from the cfg file
            %OPTtxt = readInFile('OPT',OPTfile,optid);
            %Determine the OPT properties based on the text in the OPT file
            %obj.determineOPT(OPTtxt,[modelobj.GVARtxt; readInFile('VAR',OPTfile,1)]);
            
            
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
            obj.optStrategy = details(1:ind-1);
            obj.QNupdate = details(ind+1:end);
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
        
        function [] = activeSetMethod(obj,x0)
            %This function runs the QuasiNewton NAND optimization loop.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj - QuasiNewton object
            %p   - Initial point for optimization algorithm 
            %
            %Outputs:
            %--------
            %No outputs.  Data stored in QuasiNewton object.
            %--------------------------------------------------------------
            
            tStartPDEopt = tic; %Timer
            
            n = length(x0);
            x = computeFeasible(obj,x0);
            actInd = determineActive(obj,x);
            
            
            
            
            while 1
                
                [~,ActIntWork,~] = intersect(actInd,obj.I);
                kkt = [obj.G -obj.A(actInd,:)'; obj.A(actInd,:), zeros(length(actInd))]\[-obj.G*x+obj.c;zeros(length(actInd),1)];
                p = kkt(1:n);
                lam = kkt(n+1:end);
                
                if norm(p) < 1e-14
                    
                    posLM = lam >= 0;
                    
                    if sum(posLM(ActIntWork,1)) == length(ActIntWork)
                        break;
                    else
                        [~,minInd] = min(lam(ActIntWork));
                        actInd = unique([actInd;minInd(1)]);
                    end
                else
                    
                    val = inf;
                    blkconst = [];
                    for i = setdiff(1:obj.m,actInd(:)')
                        aTp = obj.A(i,:)*p;
                        if aTp < 0
                            val = min(val,(b(i) - obj.A(i,:)*x)/aTp);
                        end
                        
                        if val < 1
                            blkconst = [blkconst;i];
                        end
                    end
                    alpha = min(val,1);
                    
                    x = x + alpha*p;
                    
                    actInd = unique([actInd;blkconst]);
                end                
            end
            
            obj.time = toc(tStartPDEopt);

%             %Set up the problem residual and jacobian with the appropriate
%             %parameter
%             p = p0(:);
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
        end
                
        function  [p0] = computeFeasible(obj,pEst)
            %This function computes a feasible point of the QP subproblem.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj - QP object
            %
            %Outputs:
            %--------
            %p0  - feasible starting point
            %--------------------------------------------------------------
            
            p0 = 0;
        end
        
        function  [aI] = determineActive(obj,p)
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
                        lineQuasiNewtonNAND(obj,p);
                    elseif strcmpi(obj.optStrategy,'TR')
                        %Trust region strategy
                        trustQuasiNewtonNAND(obj,p);
                    end
                case 'SAND'
                    %NAND optimization framework
                    if strcmpi(obj.optStrategy,'LS')
                        %Line search strategy
                        lineQuasiNewtonSAND(obj,p);
                    elseif strcmpi(obj.optStrategy,'TR')
                        %Trust region strategy
                        trustQuasiNewtonSAND(obj,p);
                    end
            end
        end
    end
end
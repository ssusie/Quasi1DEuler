classdef QuasiNewtonOLD < handle
    properties
        method = 'QN - BFGS';
        LSflag = 'linesearch';
        model = [];
        paramBounds = [];
        paramHist = [];
        solnHist = [];
        nSolverCalls = 0;
        PDECallsParam = 0;
        
        dp = []; %Store the current search direction of the parameter (for linesearch purposes)
        
        sensSolnMethod = 'adjoint'; %or 'direct'
        time = 0;
    end
    
    
    methods
   
        function [obj] = QuasiNewtonOLD(method,model,paramBounds,pstar)
            obj.model = model;
            obj.method = method;
            obj.paramBounds = paramBounds;
            
            obj.model.prob.updateParameters(pstar(:));
            obj.model.executeModel;
            obj.model.prob.setTargetSolution(obj.model.sv(:,2));
        end
        
        function  [] = determineOPT(obj,OPTtext,VARtxt)
            %This function computes interprets the text in OPTtext
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - instance of FOM class
            %OPTtxt  - a string containing the text from the OPT input file
            %VARtxt  - a string containing the text from the VAR block in
            %          the OPT file
            %
            %Outputs:
            %--------
            %No outputs.  Properties are stored in FOM handle class.
            %--------------------------------------------------------------
            
            %Extract variables from VARtxt and evaluate them in the current
            %workspace.
            if ~isempty(VARtxt)
                for i = 1:length(VARtxt)
                    eval([VARtxt{i},';']);
                end
            end
            
            
            obj.id = extractModelPropMultInput(OPTtext,1,1,'id',[]);
            obj.type = extractModelPropMultInput(OPTtext,1,1,'type',[]);
            
            
            %Determine folder
            obj.folder = extractModelPropMultInput(FOMtext,1,1,'folder',obj.folder);
            %Determine whether to save the state vector at all nonlinear
            %iterations (useful for collecting steady snapshots)
            obj.saveAllIt = extractModelPropMultInput(FOMtext,1,1,'saveAllIt',false); %false is default
            %Determine the upper and lower bounds on the time interval
            obj.time.T = extractModelPropMultInput(FOMtext,1,1,'T',obj.time.T);
            %Determine the time step size
            obj.time.dt = extractModelPropMultInput(FOMtext,1,1,'dt',obj.time.dt);
            %Determine the number of time steps
            obj.time.nstep = extractModelPropMultInput(FOMtext,1,1,'nstep',obj.time.nstep);
            %Determine whether or not to time progress
            obj.time.quiet = extractModelPropMultInput(FOMtext,1,1,'timeQuiet',obj.time.quiet);
            
            %Compute missing time quantity
            if isempty(obj.time.dt) && isempty(obj.time.nstep)
                %User needs to specify either nstep or dt
                error('Must specify either time increment or number of time steps');
            elseif isempty(obj.time.dt) && ~isempty(obj.time.nstep) %Calculate dt from nstep
                obj.time.dt = (obj.time.T(2) - obj.time.T(1))/obj.time.nstep;
            elseif isempty(obj.time.nstep) && ~isempty(obj.time.dt) %Calculate nstep from dt
                obj.time.nstep = (obj.time.T(2) - obj.time.T(1))/obj.time.dt;
            else %if both are specified, return an error
                if obj.time.dt ~= (obj.time.T(2) - obj.time.T(1))/obj.time.nstep
                    error('time increment (dt) and number of time steps (nstep) both specified and are not consistent (dt must equal (T(2) - T(1))/nstep)');
                end
            end
            
            %Determine the time stepping scheme
            TimeSchemeString = extractModelPropMultInput(FOMtext,1,1,'TimeScheme',[]);
            obj.TimeScheme = determineTimeScheme(TimeSchemeString,obj);
            
            %Determine the maximum number(s) of newton iterations
            obj.newt.maxIter = extractModelPropMultInput(FOMtext,1,1,'maxIter',obj.newt.maxIter);
            %Determine the absolute tolerance(s) for newton convergence
            obj.newt.eps = extractModelPropMultInput(FOMtext,1,1,'eps',obj.newt.eps);
            %Determine whether or not to print Newton warnings
            obj.newt.quiet = extractModelPropMultInput(FOMtext,1,1,'newtQuiet',obj.newt.quiet);
            %Initialize vector that will hold the number of newton
            %iterations required at each time step
            obj.newt.iter = zeros(1,obj.time.nstep+1);
            
            %Determine the linesearch properties
            obj.newt.linesrch.status = extractModelPropMultInput(FOMtext,1,1,'linesrch',false);
            obj.newt.linesrch.prop   = extractModelPropMultInput(FOMtext,1,1,'linesrchProp',[]);
            if isempty(obj.newt.linesrch.prop)
                obj.newt.linesrch.prop = [1e-4,0.9,50,5,10,0.05,0.1,1.2];
            end
            
            %Determine the convergence criterion to use
            if length(obj.newt.eps) == 1
                obj.newt.converge = 1;
            elseif length(obj.newt.eps) == 2
                obj.newt.converge = 2;
            elseif length(obj.newt.eps) == 3
                if obj.newt.eps(1) == 0
                    obj.newt.converge = 2;
                elseif obj.newt.eps(2) == 0 && obj.newt.eps(3) == 0
                    obj.newt.converge = 1;
                else
                    obj.newt.converge = 3;
                end
            else
                error(['eps can be either 1 x 1 (use relative residual as newton convergence criteria) ',...
                    '1 x 2 (use absolute residual and absolute iterate distance, respectively, as newton convergence criteria ',...
                    '1 x 3 (use which ever of these two convergence criteria is encountered first (format in the case = [eps_res_rel, eps_res_abs, eps_it_abs']);
            end
            
            %Determine which residual and jacobian snapshots to save (if
            %any)
            obj.saveNL = extractModelPropMultInput(FOMtext,1,1,'saveNL',obj.saveNL);
            
            
            %Determine and store the number of nodes in the FOM
            obj.ndof = extractModelPropMultInput(FOMtext,1,1,'ndof',obj.ndof);
            if length(obj.ndof) > 1
                error('ndof must be a scalar indicating the number of unknowns in the problem');
            end
            %Set up state vector matrix and store initial condition
            obj.sv = zeros(obj.ndof,obj.time.nstep+1);
        end
        
        
        function  [DfDp,f,lambda] = reducedGrad(obj,w,dcdw,dcdp)
            
            [f,dfdw,dfdp] = obj.model.prob.objective(w);
            switch obj.sensSolnMethod
                case 'direct'
                    DwDp = -dcdw\dcdp;
                    DfDp = dfdp + DwDp'*dfdw;
                    lambda = [];
                case 'adjoint'
                    lambda = -dcdw'\dfdw;
                    DfDp = dcdp'*lambda;
            end
        end
        
        function [LSfun,LSder] = LineSearchNANDQN(obj,alpha)
            
            if alpha ~= 0
                obj.model.prob.updateParameters(obj.paramHist(:,end)+alpha*obj.dp);
                obj.model.resetInitialCondition(obj.model.prob.ic);
                %obj.model.resetInitialCondition(obj.model.sv(:,2));
                obj.model.executeModel;
            end
            [dcdw,dcdp] = obj.model.prob.ResSens(obj.model.sv(:,2));
            [DfDp,f] = obj.reducedGrad(obj.model.sv(:,2),dcdw,dcdp);
            
            LSfun = f;
            LSder = DfDp'*obj.dp;
        end
        
        function [] = QuasiNewtonNAND(obj,p)
            tStartPDEopt = tic;
            %Set up the problem residual and jacobian with the appropriate
            %parameter
            p = p(:);
            obj.model.prob.updateParameters(p);
            %Run the model to solve the constraints
            obj.model.executeModel;
            obj.solnHist = obj.model.sv(:,2);
            [dcdw,dcdp] = obj.model.prob.ResSens(obj.solnHist);
            
            %Record the number of Solver calls and the parameter evolution
            obj.nSolverCalls = 1;
            obj.PDECallsParam = p; %Includes linesearch parameters
            
            obj.paramHist = p; %Does not include linesearch parameters
            
            % quasi-Newton
            eps = 1e-6;
            eps_dp = 1e-6;
            %alpha_max0 = 1;
            DfDp = obj.reducedGrad(obj.solnHist(:,end),dcdw,dcdp);
            H = speye(size(DfDp,1));
            %p_iterates = p;
            while (norm(DfDp,2)>eps )
                obj.dp = -H*DfDp; % BFGS on inverse Hessian
                disp(['gradient norm = ', num2str(norm(DfDp,2))]);
                % do not allow iterates larger than pmax or smaller than pmin
                %                 if (length(p)>1)
                %                     disp('*** Error: step limiter in NAND QN is only valid for 1-d parameter space: extend to box bounds ***');
                %                     return;
                %                 end
                
                %This selection of alpha_max ensures that
                %p_low <= p + alpha_max*dp <= p_up.
                %Rearranging this equation, we get that alpha_max must
                %satsify the following two conditions:
                %alpha_max <= ((p_up - p)./dp)_j   if dp_j > 0
                %alpha_max <= ((p - p_low)./dp)_j  if dp_j < 0
                alpha_max = min(...
                    (obj.dp > 0).*(obj.paramBounds(:,2) - p)./obj.dp + ...
                    (obj.dp < 0).*(obj.paramBounds(:,1) - p)./obj.dp);
                alpha_max = min(alpha_max,1);
                %                 if (obj.dp>0)
                %                     alpha_max = min((obj.paramBounds(1,2)-p)/obj.dp,alpha_max0);
                %                 elseif (obj.dp==0)
                %                     break;
                %                 else
                %                     alpha_max = min((obj.paramBounds(1,1)-p)/obj.dp,alpha_max0);
                %                 end
                
                switch obj.LSflag
                    case 'newton'
                        alpha = alpha_max;
                    case 'linesearch'
                        [alpha,~,nFunCalls,alphaCalls] = linesearch(@(alphaloc) obj.LineSearchNANDQN(alphaloc),[10^(-4);0.9;alpha_max;1000;1000;0.05;0.1;1.2]);
                        obj.PDECallsParam = [obj.PDECallsParam,repmat(p,1,size(alphaCalls,2))+repmat(obj.dp,1,size(alphaCalls,2)).*repmat(alphaCalls,size(obj.dp,1),1)];
                end
                obj.nSolverCalls = obj.nSolverCalls + nFunCalls;
                obj.dp = alpha*obj.dp;
                p = p + obj.dp;
                disp('Current iteration: p = ');
                disp(p);
                obj.paramHist = [obj.paramHist,p];
                if (abs(obj.dp/p)<eps_dp)
                    break;
                end
                
                %                 %Set up the problem residual and jacobian with the appropriate
                %                 %parameter
                %                 obj.model.prob.updateParameters(p);
                %                 %Run the model to solve the constraints
                %                 obj.model.resetInitialCondition(obj.solnHist(:,end));
                %                 %obj.model.resetInitialCondition(obj.model.prob.ic);
                %                 obj.model.executeModel;
                
                %We don't run the model at the new value of p because this
                %has ALREADY BEEN DONE IN THE LINESEARCH.  It would be a
                %waste to do it again.  Therefore, we simply update the
                %values of paramHist and solnHist according to this LAST
                %VALUE OF alpha from the linesearch.
                obj.solnHist = [obj.solnHist,obj.model.sv(:,2)];
                [dcdw,dcdp] = obj.model.prob.ResSens(obj.solnHist(:,end));
                
                %                 obj.PDECallsParam = [obj.PDECallsParam,p];
                %                 obj.nSolverCalls = obj.nSolverCalls + 1;
                
                DfDp_old = DfDp;
                DfDp = obj.reducedGrad(obj.solnHist(:,end),dcdw,dcdp);
                s = obj.dp;
                y = DfDp - DfDp_old;
                H = obj.updateInverseHessian(H,y,s);
            end
            obj.time = toc(tStartPDEopt);
        end
        
        function [H] = updateInverseHessian(obj,Hp,y,s)
            switch obj.method
                case 'QN - BFGS'
                    rho = 1/(y'*s);
                    if (abs(rho)>1e-16)
                        M = speye(size(Hp)) - rho*s*y';
                        H = M*Hp*M' + rho*(s*s');
                    end
                case 'QN - DFP'
                case 'QN - SR1'
            end
        end
    end
end
classdef FOM < handle
    
    properties (SetAccess = private, GetAccess = public)
        %Input properties
        id;
        fileRoot = [pwd,SlashBS];
        prob;
        ndof;
        time;
        newt;
        saveNL = 0;
        TimeScheme;
        saveAllIt = false;

        curr_param;
        objective;
        constraints;
        
        %Computed Properties
        sv = [];
        svIter = [];
        res = [];
        ontime = 0;
        cTimeIter = 0; %current time iteration number
        
        numExecute=0;
        Cnorm;
    end
    
    properties (Hidden = true, SetAccess = private, GetAccess = public)
        %Global variable text
        GVARtxt = [];
        
        printLevel=3;
        killflag=false;
        numres = 0;
    end
    
    methods
        function  [obj] = FOM(CFGfile,cfgid,probobj)
            %This is the constructor of the FOM class.  It reads the
            %appropriate entries from the input file and stores them in the
            %class instance.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %CFGfile - string indicating the filename of the cfg file to
            %          use
            %cfgid   - integer indicating the CONFIG object in CFGfile to
            %          use
            %probobj - instance of a problem object containing
            %          problem-specific information for the analysis
            %
            %Outputs:
            %--------
            %obj     - instance of the FOM object that was constructed
            %--------------------------------------------------------------
            
            if nargin == 0
                return;
            end
            
            %Store handle to problem in instance
            obj.prob = probobj;
            
            %Initialize time structure
            obj.time = probobj.config.time;
            obj.newt = struct('maxIter',[],'eps',[],'iter',[],'quiet',[]);
            
            % %%%%%%%%%%%%%%%%%%%%Set default values%%%%%%%%%%%%%%%%%%%%
            % dfname = [probobj.config.defaultFile,'.cfg'];
            % dVARtext = readInFile('VAR',dfname,1);
            % %Extract the FOM text from the cfg file
            % dFOMtext    = readInFile('FOM',dfname,1);
            % %Extract the CONFIG text with id number = cfgid from cfg file
            % dCONFIGtext = readInFile('CONFIG',dfname,1);
            % %Determine the FOM properties based on the text in the FOM and
            % %CONFIG sections of the input file
            % determineFOM(obj,[dFOMtext;dCONFIGtext],dVARtext);
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %Extract parameters from CFG file
            obj.GVARtxt = probobj.config.GVARtxt;
            VARtext = readInFile('VAR',CFGfile,1);
            %Extract the FOM text from the cfg file
            FOMtext    = readInFile('FOM',CFGfile,cfgid);
            %Extract the CONFIG text with id number = cfgid from cfg file
            CONFIGtext = readInFile('CONFIG',CFGfile,cfgid);
            
            %Determine the FOM properties based on the text in the FOM and
            %CONFIG sections of the input file
            determineFOM(obj,[FOMtext;CONFIGtext],[obj.GVARtxt;VARtext]);
            %%%%%%%%% Set up state vector matrix and store initial condition %%%%%%%%%
            if isempty(obj.ndof), obj.ndof = length(probobj.ic); end;
%             obj.sv = zeros(obj.ndof,obj.time.nstep+1);
%             %Set the first state vector to the initial condition
%             obj.sv(:,1) = probobj.ic;
        end
        
        function  [] = determineFOM(obj,FOMtext,VARtxt)
            %This function computes and stores properties of FOM class from
            %the char array FOMtext
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - instance of FOM class
            %FOMtxt  - a string containing the text from the FOM and CONFIG
            %          objects in the input file
            %VARtxt  - a string containing the text from the VAR block in
            %          the CFG file
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
            
            %%%%%%%%% Determine id %%%%%%%%%
            obj.id = extractModelPropMultInput(FOMtext,1,1,'id',[]);
            if isempty(obj.id)
                error('In CFG file, the ID must be specified.  This controls the connection between different aspects of the program.');
            end
            
            %%%%%%%%% Determine number of dofs in FOM %%%%%%%%%
            obj.ndof = extractInputRobust(FOMtext,'ndof',obj.ndof);
            if length(obj.ndof) > 1
                error('ndof must be a scalar indicating the number of unknowns in the problem');
            end
            
            %%%%%%%%% Determine fileRoot %%%%%%%%%
            obj.fileRoot = [extractInputRobust(FOMtext,'fileRoot',obj.fileRoot),'_ID',num2str(obj.id),'_'];
            
            %%%%%%%%% Determine whether to save residual snapshots %%%%%%%%%
            obj.saveNL = extractInputRobust(FOMtext,'saveNL',obj.saveNL);
            
            %%%%%%%%% Determine whether to save state at all newton iterates %%%%%%%%%
            obj.saveAllIt = extractInputRobust(FOMtext,'saveAllIt',obj.saveAllIt);
            
            %%%%%%%%% Determine time structure %%%%%%%%%
            obj.time = struct('T',[],'dt',[],'nstep',[],'quiet',[],'newmark',[]);
            obj.time.T     = extractInputRobust(FOMtext,'T',[]);
            obj.time.dt    = extractInputRobust(FOMtext,'dt',[]);
            obj.time.nstep = extractInputRobust(FOMtext,'nstep',1000);
            obj.time.quiet = extractInputRobust(FOMtext,'timeQuiet',false);
            obj.time.steadyconverge = extractInputRobust(FOMtext,'steadyconverge',false);
            obj.time.cfl   = extractInputRobust(FOMtext,'cfl',[]);
            obj.time.genalpha = extractInputRobust(FOMtext,'genalpha',[]);
            
            if isempty(obj.time.dt) && isempty(obj.time.cfl)
                obj.time.dt = (obj.time.T(2) - obj.time.T(1))/obj.time.nstep;
            end
%             %Compute missing time quantity
%             if isempty(obj.time.T) && ~isempty(obj.time.dt) && ~isempty(obj.time.nstep) %Calculate T from dt, nstep
%                 obj.time.T = [0,obj.time.nstep*obj.time.dt]; %Assume start time is zero
%             elseif ~isempty(obj.time.T) && isempty(obj.time.dt) && ~isempty(obj.time.nstep) %Calculate dt from T, nstep
%                 obj.time.dt = (obj.time.T(2) - obj.time.T(1))/obj.time.nstep;
%             elseif ~isempty(obj.time.T) && ~isempty(obj.time.dt) && isempty(obj.time.nstep) %Calculate nstep from T, dt
%                 obj.time.nstep = ceil((obj.time.T(2) - obj.time.T(1))/obj.time.dt);
%             else
%                 %User needs to exactly 2 of the 3 time quantities
%                 error('Must specify exactly 2 of the 3 following fields: T, dt, nstep');
%             end
            
            %%%%%%%%% Determine time integration %%%%%%%%%
            TimeSchemeString = extractInputRobust(FOMtext,'TimeScheme',[]);
            obj.TimeScheme   = determineTimeScheme(TimeSchemeString,obj);
            
            obj.printLevel = extractInputRobust(FOMtext,'printLevel',obj.printLevel);

            %%%%%%%%% Determine newton structure %%%%%%%%%
            obj.newt = struct('maxIter',[],'eps',[],'quiet',[],'iter',[],'linesrch',[],'converge',[],'avgIter',[]);
            obj.newt.maxIter = extractInputRobust(FOMtext,'maxIter',10);
            obj.newt.eps     = extractInputRobust(FOMtext,'eps',1e-5);
            obj.newt.quiet   = extractInputRobust(FOMtext,'newtQuiet',false);
            obj.newt.iter    = zeros(1,obj.time.nstep);
            %Determine linesearch status and properties
            obj.newt.linesrch.status = extractInputRobust(FOMtext,'linesrch',false);
            obj.newt.linesrch.prop   = extractInputRobust(FOMtext,'linesrchProp',[1e-4,0.9,50,5,10,0.05,0.1,1.2]);
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
        end
        
        function  [] = NewtonRaphson(obj)
            %This function solves systems of the form:  f(x) = 0 using
            %Newton-Raphson's method
            %--------------------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - FOM object
            %
            %Outputs:
            %--------
            %There are no outputs.  The state vector at the current time
            %(indicated by the cTimeIter property) is stored in the FOM
            %handle class.
            %--------------------------------------------------------------------------
            
            %Determine time iteration number plus 1
            itnump1 = obj.cTimeIter + 1;
            
            %Determine the residual and jacobian based on initial guess
            t = obj.time.T(1) + obj.time.dt*obj.cTimeIter;
            [R,J] = obj.TimeScheme.TimeIntNLFunc(obj.prob,obj.sv(:,itnump1-1),obj.sv(:,itnump1-1),t);
            %Write residual to file (if requested)
            obj.writeResidual(R);
                
            %Determine convergence criteria to use
            [tol,tolIt] = determineConvergeCriterion(obj,norm(R,2));
            
            indexAdj=1;
            %if obj.printLevel > 1
            %    fprintf('---- Newton Step    # --------- ||R||  ------------ ||du|| ------------- tol -------------- tol_du ---------\n');
            %end 
            for i_N = 1:obj.newt.maxIter
                %Determine step
                p = -J\R;
                %setup.type='nofill'; [L,U]=ilu(J,setup);
                %p = -gmres(J,R,[],[],100,L,U);
                if obj.newt.linesrch.status
                    alpha = linesearch(@(alpha) obj.TimeScheme.TimeIntNLFuncDir(obj.prob,obj.sv(:,itnump1-1),obj.sv(:,itnump1-1),t,alpha,p),obj.newt.linesrch.prop);
                    p = alpha*p;
                end
                
                %if obj.printLevel > 2
                %    fprintf('---- Newton Step %4i ----- %10.7e ----- %10.7e ----- %10.7e ----- %10.7e ------\n',i_N,norm(R,2),norm(p,2),tol,tolIt);
                %end
                
                %If nan or complex encountered -> kill simulation
                if sum(isnan(p)) || ~isreal(p)
                    obj.killflag=true;
                    return;
                 end
                
                %Update solution
                obj.sv(:,itnump1) = obj.sv(:,itnump1-indexAdj) + p;
            
                %Save Newton iterates (if requested)
                if obj.saveAllIt, obj.svIter = [obj.svIter, obj.sv(:,itnump1)]; end
                            
                %Compute residual and jacobian with updated vector for next
                %iteration.  Save residual if saveNL = true in input file.
                [R,J,iprevsv] = obj.TimeScheme.TimeIntNLFunc(obj.prob,obj.sv(:,itnump1),obj.sv(:,itnump1-1),t);
                
                %Write residual to file (if requested)
                obj.writeResidual(R);
                
                indexAdj=0;
                %Stop iterating once convergence criteria is met
                if checkConverge(obj,norm(R,2),p,tol,tolIt)
                    break;
                end
            end
            obj.prob.setPrevSV(iprevsv);
            %Store the number of newton iterations at this time step
            obj.newt.iter(obj.cTimeIter) = i_N;
            
            g = obj.TimeScheme.constraintOneDomain(obj.prob,obj.sv(:,itnump1),obj.sv(:,itnump1-1),t);
            obj.Cnorm = [obj.Cnorm, norm(g)];
            %if obj.printLevel == 1.5
            %    fprintf('---- Newton Step %4i ----- %10.7e ----- %10.7e ----- %10.7e ----- %10.7e ------\n',i_N,norm(R,2),norm(p,2),tol,tolIt);
            %end
            
            %If the maxiumum newton iterations are reached, warn the user
            if ~obj.newt.quiet && (i_N==obj.newt.maxIter)
                disp(['*** Warning; Newton solver reached max number of iterations before convergence : res = ',...
                    num2str(norm(R,2)/tol),'  eps = ',num2str(obj.newt.eps), '***']);
            end
        end

        function  [] = executeModel(obj,restart)
            %This function performs the time-stepping for the FOM
            %simulation defined in this object.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - FOM object
            %
            %Outputs:
            %--------
            %There are no outputs.
            %--------------------------------------------------------------
            obj.numExecute = obj.numExecute+1;
            
            obj.killflag=false;
            nstep = obj.time.nstep; %Extract the number of time steps from the object
            obj.Cnorm = [];
            
            %Restart logic
            if nargin == 2 && restart == 1
                OldTime = obj.ontime;
                firstStep = obj.cTimeIter-1;
            else
                OldTime = 0;
                firstStep = 1;
                obj.svIter = [];
                
                %Set up state vector and store the initial condition.
                obj.sv = zeros(obj.ndof,obj.time.nstep+1);
                obj.sv(:,1) = obj.prob.ic;
                
                obj.numres=0;
                %Open file to save nonlinear terms if requsted
                if obj.saveNL > 0
                    obj.res = fopen([obj.fileRoot,'ResSnaps_basis1','.bin'],'wb');
                end
            end
            
            %Set the current parameter (from prob object) to be used in
            %NAND to know which parameter values have been used to generate
            %the current state vectors
            obj.curr_param = obj.prob.p;
            
            tFOMstart = tic;
            for i_t = firstStep:nstep %Loop over each time step to determine all state vectors
                if ((obj.printLevel>0)&&(rem(i_t,round(nstep*0.1)) == 0)) || (obj.printLevel>1)
                    fprintf('-------------------------- Time Step %4i of %4i (%2i%%) -------------------------- \n',i_t,nstep,ceil(100*i_t/nstep));
                end
                
                %Set current time iteration (so model can keep track)
                obj.cTimeIter = i_t;
                
                if obj.TimeScheme.explicit
                    obj.sv(:,i_t+1) = obj.TimeScheme.ExplicitStep(obj.prob,obj.sv(:,i_t),obj.time.T(1)+obj.TimeScheme.dt*i_t);
                    if sum(isnan(obj.sv(:,i_t+1)))>0
                        fprintf('NaN encountered on time step %i!\n',i_t);
                        obj.killflag=true;
                        return;
                    end
                    continue;
                end
                
                if ~isempty(obj.time.cfl)
                    obj.time.dt = obj.time.cfl(i_t,obj)*obj.prob.computeCFL(obj.sv(:,i_t));
                    obj.TimeScheme.setDT(obj.time.dt);
                end
                
                %Use Newton's Method to solve nonlinear system of equations
                %and store: state vector, residual and jacobian snapshots,
                %and number of newton iterations required
                obj.NewtonRaphson();
                if obj.killflag, warning('Encountered NaN at Time Step %i. Simulation Killed',obj.cTimeIter); return; end;

%                 if obj.time.steadyconverge
%                     if norm(obj.sv(:,i_t+1)-obj.sv(:,i_t)) < obj.time.steadyconverge*norm(obj.sv(:,i_t+1))
%                         obj.ontime = toc(tFOMstart) + OldTime;
%                         obj.time.nstep=i_t;
%                         obj.time.T(2) = obj.time.T(1) + i_t*obj.time.dt;
%                         obj.sv(:,i_t+2:end)=[];
% %                         obj.prob.config.setTimeObj('nstep',i_t);
% %                         obj.prob.config.setTimeObj('T',[obj.time.T]);
%                         return;
%                     end
%                 end
            end
            obj.ontime = toc(tFOMstart) + OldTime; %Record simulation time
            
            %Close file, if necessary
            if obj.saveNL > 0
                fclose(obj.res);
            end
            obj.newt.avgIter = mean(obj.newt.iter);
        end
        
        function  [] = determineNLSnapshot0p5(obj,romobj)
            %This function determines the nonlinear snapshots using method
            %0.5 and saves them to files.  Requires ROM object because
            %reduced basis is used.
            
            if ~obj.saveAllIt %For snapshot 0.5 to make sense, must collect all iterates during FOM
                fprintf('Must have saveAllIt = true in FOM to use Snapshot 0.5. Exiting...');
                return;
            end
            
            if romobj.nBases == 1
                
                obj.numres=0; %Track the number of residual vectors saved
                obj.res = fopen([obj.fileRoot,'ResSnaps_basis1.bin'],'wb');
                
                for i = 1:obj.time.nstep
                    t = obj.time.T(1) + obj.time.dt*i;
                    %Compute the orthogonal projection of the state
                    %vector on the reduced basis
                    svPrevProj = obj.prob.ic + romobj.phi*(romobj.phi'*(obj.sv(:,i)-obj.prob.ic));
                    
                    for j = 1:obj.newt.iter(i)
                        k = j + sum(obj.newt.iter(1:i-1));
                        
                        %Compute the orthogonal projection of the
                        %Newton iterate on the reduced basis.
                        svIterProj = obj.prob.ic + romobj.phi*(romobj.phi'*(obj.svIter(:,k)-obj.prob.ic));
                        
                        R = obj.TimeScheme.TimeIntNLFunc(obj.prob,svIterProj,svPrevProj,t);
                        obj.numres=obj.numres+1;
                        fwrite(obj.res,R,'double');
                    end
                end
                fclose(obj.res);
                
            else
                
                if ~(strcmpi(romobj.basisUpdate,'border_fast') || strcmpi(romobj.basisUpdate,'none'))
                    warning('Must have basisUpdate = ''border_fast'' or ''none'' in ROM to use Snapshot 0.5.  Exiting...');
                    return;
                end
                %Set numres to zero for each basis and open new files (for each
                %basis) for writing
                obj.numres = zeros(romobj.nBases,1);
                for i = 1:romobj.nBases
                    obj.res(i) = fopen([obj.fileRoot,'ResSnaps_basis',num2str(i),'.bin'],'wb');
                end
                
                wtildePrev = obj.prob.ic;
                prevBase = NaN;
                
                for i = 1:obj.time.nstep  
                    %Compute time
                    t = obj.time.T(1) + obj.time.dt*i;
                    
                    %Determine closest basis
                    cdist= ColumnwiseNorm(bsxfun(@minus,obj.sv(:,i),romobj.clusCenter),2);
                    %cdist= ColumnwiseNorm(bsxfun(@minus,wtildePrev,romobj.clusCenter),2);
                    [~,closCluster] = min(cdist);
                    closCluster=closCluster(1);
                                        
                    %Initialize ROB update components
                    if (i == 1 && strcmpi(romobj.basisUpdate,'border_fast'))
                        romobj.initializeROBcompon(closCluster);
                    end
                    
                    %Either update basis or just choose the appropriate Phi
                    if strcmpi(romobj.basisUpdate,'border_fast') && ((i > 1 && closCluster ~= prevBase)||( i == 1 && ~strcmpi(romobj.snaps.initref,'ic')))
                        %If it is AFTER the first timestep and the current basis
                        %differs from the previous basis OR if it is the first time
                        %step and the initial snapshot reference was NOT the
                        %initial condition, perform fast SVD update
                        romobj.efficientOnlineSVDupdate(closCluster,wtildePrev);
                    elseif (i > 1 && strcmpi(romobj.basisUpdate,'none') && closCluster ~= prevBase)
                        romobj.setPhi(closCluster);
                    elseif ( i == 1 )
                        romobj.setPhi(closCluster);
                        romobj.setNumSwitch(1);
                    end
                                                            
                    %Include the "initial guess" residual evaluation
                    R = obj.TimeScheme.TimeIntNLFunc(obj.prob,wtildePrev,wtildePrev,t);
                    obj.numres(closCluster)=obj.numres(closCluster)+1;
                    fwrite(obj.res(closCluster),R,'double');
                    for j = 1:obj.newt.iter(i)
                        %Determine the indices in svIter that should be written.  We
                        %want all newton iterations for the current time step.
                        k = j + sum(obj.newt.iter(1:i-1));
                        
                        wtilde     = wtildePrev + romobj.phi*(romobj.phi'*(obj.svIter(:,k)-wtildePrev));
                        
                        R = obj.TimeScheme.TimeIntNLFunc(obj.prob,wtilde,wtildePrev,t);
                        obj.numres(closCluster)=obj.numres(closCluster)+1;
                        fwrite(obj.res(closCluster),R,'double');
                    end
                    wtildePrev = wtilde;
                    prevBase = closCluster;
                end
                
                for i = 1:romobj.nBases
                    fclose(obj.res(i));
                end
            end
        end
        
        function  [] = sortNLSnapshot0byBasis(obj,romobj)
            %This function sorts the residual snapshots into the
            %appropriate file based on proximity of the corresponding state
            %to each cluster center.  ROM must be an input.
            
            if romobj.nBases == 1
                return;
            end
            
            %Read the residual files from file
            R = obj.readResJacFiles(0,1);
            
            %Set numres to zero for each basis and open new files (for each
            %basis) for writing
            obj.numres = zeros(romobj.nBases,1);
            for i = 1:romobj.nBases
                obj.res(i) = fopen([obj.fileRoot,'ResSnaps_basis',num2str(i),'.bin'],'wb');
            end
            
            for i = 1:obj.time.nstep
                %Determine cluster closest to current solution
                cdist= ColumnwiseNorm(bsxfun(@minus,obj.sv(:,i),romobj.clusCenter),2);
                [~,closCluster] = min(cdist);
                closCluster=closCluster(1);
                                
                %Determine the indices in R that should be written.  We
                %want all newton iterations for the current time step.  The
                %"plus 1" in both terms come from the fact that we save 1
                %more residual than newton iterations PER TIME STEP
                k = (1:obj.newt.iter(i)+1) + sum(1+obj.newt.iter(1:i-1));
                
                %Update numres and write the residual to the appropriate
                %file
                obj.numres(closCluster)=obj.numres(closCluster)+1;
                fwrite(obj.res(closCluster),R(:,k),'double');
            end
            
            %Close residual files for each basis
            for i = 1:romobj.nBases
                fclose(obj.res(i));
            end
            
        end
        
        function  [] = writeResidual(obj,R)
            %This function writes the residaul R to the file (assumed
            %opened) if saveNL = 1.
            if obj.saveNL == 1
                obj.numres=obj.numres+1;
                fwrite(obj.res,R,'double');
            end
        end
        
        function  [res,jac] = readResJacFiles(obj,NLSnapColl,basisNum)
            %This function reads the binary files containing the residual
            %and/or jacobian snapshots.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj      - ROM object
            %basisNum - scalar indicating the basis to read the nonlinear
            %           snapshots from.  If using global ROM, this is
            %           automatically set to zero.
            %
            %Outputs:
            %--------
            %res      - matrix containing the residual snapshots
            %jac      - []
            %--------------------------------------------------------------
            
            res=[]; jac=[];
            
            if NLSnapColl == 0 || NLSnapColl == 0.5
                fidR = fopen([obj.fileRoot,'ResSnaps_basis',num2str(basisNum),'.bin'],'r');
                if fidR == -1, disp('File couldn''t be opened'); return; end;
                
                res = reshape(fread(fidR,obj.ndof*obj.numres(basisNum),'double'),obj.ndof,obj.numres(basisNum));
                jac = [];
                
                fclose(fidR);
            end
        end
        
        function  [newobj] = hardCopy(obj,idnum)
            % Construct a new object based on a deep copy of an old object
            % of this class by copying properties over.
            newobj=FOM();
            
            props = properties(obj);
            for i = 1:length(props)
                % Use Dynamic Expressions to copy the required property.
                % For more info on usage of Dynamic Expressions, refer to
                % the section "Creating Field Names Dynamically" in:
                % web([docroot '/techdoc/matlab_prog/br04bw6-38.html#br1v5a9-1'])
                newobj.(props{i}) = obj.(props{i});
            end
            
            newobj.id = idnum;
            str=['_ID',num2str(idnum),'_'];
            startind = regexp(newobj.fileRoot,'_ID[0-9]*_','start');
            endind   = startind + length(str) - 1;
            newobj.fileRoot(startind:endind) = str;
            
            %newobj.res = fopen([newobj.fileRoot,'ResSnaps_basis1.bin'],'wb');
            newobj.numres=0;
        end
        
        function  [] = setProperty(obj,prop1,prop2,val)
            if isempty(prop2)
                obj.(prop1)=val;
            else
                obj.(prop1).(prop2)=val;
            end
        end
        
        %Optimization
        function  [R, dRdw, dRdp] = PDEconstraint(obj, w, p)
            [R,dRdw,dRdp] = obj.prob.ResSens(w, p);
        end
        
        function   [f,df] = ObjFMINCON_SAND(obj,u)
            
            w=u(1:obj.ndof,1);
            p=u(obj.ndof+1:end,1);
            
            obj.prob.updateParameters(p);
            
            [f,dfdw,dfdp] = obj.objective(w, p);
            df = [dfdw;dfdp];
        end
        
        function   [f,df] = ObjFMINCON_NAND(obj,p,sensMethod,flag)
            
            if norm(p - obj.curr_param) > 0
                obj.prob.updateParameters(p);
                obj.executeModel();
            end
            if nargin < 4 || isempty(flag)
                [df,f] = obj.reducedGrad(p,sensMethod);
            elseif strcmpi(flag,'objective')
                [~,f] = obj.reducedGrad(p,sensMethod);
                df=[];
            elseif strcmpi(flag,'gradient')
                [f,~] = obj.reducedGrad(p,sensMethod);
                df=[];
            end
        end
        
        function   [cineq,ceq,dcineq,dceq] = NLConstFMINCON_SAND(obj,u)
            
            w=u(1:obj.ndof,1);
            p=u(obj.ndof+1:end,1);
            
            obj.prob.updateParameters(p);
            
            [R,dRdw,dRdp] = obj.prob.ResSens(w,[]);
            ceq = R;
            dceq = [dRdw, dRdp]';
            
            cineq=[];
            dcineq=[];
        end
                
        function  [DfDp,f,lambda,R,dRdw,dRdp] = reducedGrad(obj,p,sensMethod)
            %This function computes the reduced gradient of the objective
            %function using either the direct or adjoint method to compute
            %the sensitivity of the state with respect to the parameter.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj    - FOM object
            %w      - state vector at the current parameter iterate
            %dcdw   - sensitivity of the residual with respect to the state
            %         at the current parameter iterate
            %dcdp   - sensitivity of the residual with respect to the
            %         parameter at the current parameter iterate
            %
            %Outputs:
            %--------
            %DfDp   - reduced gradient of the objective function at current
            %         iterate
            %f      - value of the objective function at current iterate
            %lamdba - dual variable
            %R      - residual of constraints
            %--------------------------------------------------------------
            w = obj.sv(:,end);
            
            %Compute the residual sensitivities
            [R,dRdw,dRdp] = obj.prob.ResSens(w,p);
            %Evaluate objective and sensitivities
            [f,dfdw,dfdp] = obj.objective(w,p);
            %             %Compute the residual sensitivities
            %             [R,dRdw,dRdp] = obj.prob.ResSens(obj.sv(:,end),p);
            %             %Evaluate objective and sensitivities
            %             [f,dfdw,dfdp] = obj.objective(obj.sv(:,end),p);
            switch sensMethod
                case 'direct'
                    %Using the direct method, solve for the state
                    %sensitivity (w.r.t. parameter) and compute the reduced
                    %gradient
                    DwDp = -dRdw\dRdp;
                    DfDp = dfdp + DwDp'*dfdw;
                    lambda = [];
                case 'adjoint'
                    %Using the adjoint method, solve for the dual variable
                    %and compute the reduced gradient
                    lambda = -dRdw'\dfdw;
                    DfDp = dfdp + dRdp'*lambda;
            end
        end
        
        function  [] = setOptObjective(obj,txt)
            %This function computes the reduced gradient of the objective
            %function using either the direct or adjoint method to compute
            %the sensitivity of the state with respect to the parameter.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj - FOM object
            %txt - string containing the name of the objective function
            %      which is stored in obj.prob
            %
            %Outputs:
            %--------
            %There are no outputs.
            %--------------------------------------------------------------
            
            obj.objective = eval(['@(u,z) obj.prob.',txt,'(u,z)']);
        end
        
        function  [] = setOptConstraints(obj,txt)
            %This function computes the reduced gradient of the objective
            %function using either the direct or adjoint method to compute
            %the sensitivity of the state with respect to the parameter.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj - FOM object
            %txt - string containing the name of the objective function
            %      which is stored in obj.prob
            %
            %Outputs:
            %--------
            %There are no outputs.
            %--------------------------------------------------------------
            
            obj.constraints = eval(['@(u,z) obj.prob.',txt,'(u,z)']);
        end
        
        function  [A,c] = reducedConstJac(obj,z,sensMethod)
            %This function computes the reduced Jacobian of the constraints
            %using either the direct or adjoint method to compute the
            %sensitivity of the state with respect to the parameter.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj    - FOM object
            %z      - parameter at the current parameter iterate
            %Outputs:
            %--------
            %A      - Jacobian of constraint vector
            %c      - value of constraint vector
            %--------------------------------------------------------------
            
            %INEFFICIENT - REDO THIS BY TAKING ADVANTAGE OF THE FACT THAT
            %dudw was previously computed
            
            %Compute the residual sensitivities
            [~,dRdw,dRdp] = obj.prob.ResSens(obj.sv(:,2));
            [c,dcdw,dcdp] = obj.constraints(obj.sv(:,2),z);
            
            switch sensMethod
                case 'direct'
                    %Using the direct method, solve for the state
                    %sensitivity (w.r.t. parameter) and compute the reduced
                    %gradient
                    DwDp = -dRdw\dRdp;
                    A = dcdp + dcdw*DwDp;
                case 'adjoint'
                    %Using the adjoint method, solve for the dual variable
                    %and compute the reduced gradient
                    lambda = -dRdw'\dcdw';
                    A = dcdp + lambda'*dRdp;
            end
        end
        
        function  [] = resetInitialCondition(obj,U)
            obj.sv(:,1) = U;
            obj.sv(:,2:end) = 0;
        end
        
        %Postprocessing
        function  [svF] = reconstructFullState(obj,~)
            %This function reconstructs the full state vector from the
            %FOM object
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - FOM object
            %
            %Outputs:
            %--------
            %svF     - full state vector matrix
            %--------------------------------------------------------------
            
            svF = obj.sv;
        end
        
        function  [] = emptyProperties(obj,varargin)
            %This function empties the specified properties.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj      - FOM object
            %varargin - strings indicating the properties to be cleared.
            %           This is useful if there are many instances of the
            %           class and the user only wants to keep the most
            %           relevant information to avoid memory issues.  It is
            %           usually recommended to clear variables that are
            %           purely used for computation (i.e. not
            %           postprocessing) after the computations are
            %           complete.
            %
            %Outputs:
            %--------
            %There are no outputs.
            %--------------------------------------------------------------
            
            all2clear4PP = {'res'};
            
            for i = 1:length(varargin)
                %If the user wants to clear all properties not commonly
                %needed in the postprocessing phase, call this function
                %again with the appropriate inputs.
                if strcmpi(varargin{i},'allpp')
                    obj.emptyProperties(all2clear4PP{:});
                    continue;
                end
                
                eval(['obj.',varargin{i},' = []']);
            end
        end
    end
end

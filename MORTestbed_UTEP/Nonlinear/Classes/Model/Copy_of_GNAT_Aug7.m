classdef GNAT < handle
    
    properties (SetAccess = private, GetAccess = public)
        %Quantities read from input file
        id;
        prob;
        problem; %(quasi1d...)
        probGNAT;
        ic;
        nY;
        nR;
        nJ;
        nI;
        newt;
        time;
        nSample; %number of sample nodes to use
        nGreed; %number of basis vectors to use in greed node selection
        nBases;
        nSampleProd;
        addInd;
        solver=2; % =1 original, =2 constrains
        
        objectiveNAND;
        objectiveSAND;
        objectiveSAND_ROMvar;
        curr_param;
        
        Rtype;
        
        svdUpdateData;
        SingVals;
        
        %Compute quantities
        TimeScheme;
        ontime;
        sv; %State vector (reduced) nY x (nstep+1) matrix
        A; %Online matrix A from Carlberg et. al. 2011
        B; %Online matrix B from Carlberg et. al. 2011
        Ebar;
        sampleNodes; %Vector containing the sample node numbers
        sampleInd; %Vector containing the sample indice numbers
        %(will be the same as sampleNodes if there is only 1 unknown per node)
        
        irstart; %vector of indices pointing to locations in jrow to indicate the start of a new indice
        jrow; %vector of indices indicating the indices of the full state vector where the state vector needs to be evaluated
        jdiag; %boolean vector the same size as jrow indicating whether or not that component corresponds to a diagonal entry
        phiYhat; %the partial reduced order basis (i.e. phiY(jrow,:))
        cTimeIter; %current time iteration number
        
        numExecute=0;
    end
    
    properties (Hidden=true, SetAccess = private, GetAccess = public)
        %Temporary properties
        partialU; %partial state vector at the current time step
        partialUprev; %partial state vector at the previous time step
        phiYhatInd;
        JhatTemp; %matrix of size nI x number of necessary state vectors used to store the jacobian
        JhatFluxL;
        JhatFluxR;
        reconstJhatInd; %vector that is used take the partial Jacobian from vector form into matrix form (JhatTemp),
        %i.e. JhatTemp(reconstJhatInd) = Jhat will take the vector Jhat into the appropriate matrix JhatTemp
        uniqueJROWind;
        jdiagHat;
        killflag;
        GVARtxt;
    end
    
    methods
        %Constructor
        function  [obj] = GNAT(ROMfile,romobj,id,oldobj)
            %This is the constructor of the GNAT class.  It reads the
            %appropriate entries from the input file and stores them in the
            %class instance.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %ROMfile - string indicating the filename of the rom file to
            %          use
            %romobj  - ROM object that will be associated with this GNAT
            %          instance
            %oldobj  - GNAT object.  If this is nonempty, then the object
            %          instantiated will be a deep copy of oldobj.
            %
            %Outputs:
            %--------
            %obj     - instance of the GNAT object that was constructed
            %--------------------------------------------------------------
            
            if nargin == 0
                return;
            end
            
            % Construct a new object based on a deep copy of an old object
            % of this class by copying properties over.
            if nargin == 4 && isa(oldobj,'GNAT')
                copyProperties(obj,oldobj);
                return;
            end
            obj.problem=romobj;
            %Copy the time/newt structure and store in the class instance
            %The default time structure is whatever is stored in ROM
            obj.time = romobj.time;
            obj.time.maxsteps = romobj.time.nstep;
            
            obj.newt = struct('maxIter',[],'eps',[],'iter',[],'quiet',[]);
            
            %Extract parameters from ROM file
            obj.GVARtxt = romobj.GVARtxt;
            VARtext = readInFile('VAR',ROMfile,1);
            
            %Extract the GNAT text from the rom file
            GNATtext = readInFile('GNAT',ROMfile,1);
            
            %Copy other properties from the ROM
            obj.nY = romobj.nY;
            obj.nBases = romobj.nBases; %Should always be 1
            obj.Rtype = romobj.Rtype;
            
            %Determine the GNAT properties based on the text in the ROM
            %section of the input file
            determineGNAT(obj,GNATtext,[obj.GVARtxt;VARtext],id);
            
            %Inherit time integrator from rom
            timescheme=class(romobj.TimeScheme);
            obj.TimeScheme=eval([timescheme,'(obj);']);
        end
        
        function  [] = determineGNAT(obj,GNATtext,VARtxt,id)
            %This function computes and stores properties of GNAT class from
            %the char array GNATtext
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - instance of GNAT class
            %GNATtxt - a string containing the text from the GNAT object of
            %         interest in the input file
            %id      - integer indicating the id number of the GNAT object
            %          of interest from the input file
            %
            %Outputs:
            %--------
            %No outputs.  Properties are stored in GNAT handle class.
            %--------------------------------------------------------------
            
            %Extract variables from VARtxt and evaluate them in the current
            %workspace.
            if ~isempty(VARtxt)
                for i = 1:length(VARtxt)
                    eval([VARtxt{i},';']);
                end
            end
            
            %%%%%%%%% Determine and store configuration id %%%%%%%%%
            tmp = extractInputRobust(GNATtext,'id',[]);
            if isempty(tmp)
                error('In ROM file (one of the ROM sections), the ID was not specified.  This controls the connection between different aspects of the program. Cannot continue.');
            end
            idvec = tmp(:);
            N = size(idvec,1);
            
            %Determine the entry in the input file that correspond to id
            ind = find(idvec == id);
            if isempty(ind)
                error('Can only initialize ROMs with ids defined in .rom file');
            end
            obj.id = idvec(ind);
            
            %%%%%%%%% Determine nR, nJ, nI, nSample, nGreed %%%%%%%%%
            obj.nR = extractInputROM(GNATtext,N,obj.id,'nR',[]); obj.nR = obj.nR(:);
            obj.nJ = extractInputROM(GNATtext,N,obj.id,'nJ',[]); obj.nJ = obj.nJ(:);
            obj.nI = extractInputROM(GNATtext,N,obj.id,'nI',[]); obj.nI = obj.nI(:);
            obj.nSample = extractInputROM(GNATtext,N,obj.id,'nSample',[]); obj.nSample = obj.nSample(:);
            obj.nGreed  = extractInputROM(GNATtext,N,obj.id,'nGreed',[]); obj.nGreed = obj.nGreed(:);
            
            %%%%%%%%% Determine the 'free' indices to add %%%%%%%%%
            obj.addInd = extractInputROM(GNATtext,N,obj.id,'addInd','none');
            
            %%%%%%%%% Determine time structure %%%%%%%%%
            obj.time.T     = extractInputROM(GNATtext,N,obj.id,'T',obj.time.T);
            obj.time.dt    = extractInputROM(GNATtext,N,obj.id,'dt',obj.time.dt);
            obj.time.nstep = extractInputROM(GNATtext,N,obj.id,'nstep',obj.time.nstep);
            obj.time.quiet = extractInputROM(GNATtext,N,obj.id,'timeQuiet',obj.time.quiet);
            obj.time.maxsteps = extractInputROM(GNATtext,N,obj.id,'maxsteps',obj.time.maxsteps);
            
            obj.time.steadyconverge = extractInputROM(GNATtext,N,obj.id,'steadyconverge',obj.time.steadyconverge);
            obj.time.cfl = extractInputRobust(GNATtext,'cfl',obj.time.cfl);
            
            if isempty(obj.time.dt) && isempty(obj.time.cfl)
                %if isempty(obj.time.dt)
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
            
            %             %%%%%%%%% Determine time integration %%%%%%%%%
            %             TimeSchemeString = extractInputROM(GNATtext,N,obj.id,'TimeScheme',[]);
            %             obj.TimeScheme   = determineTimeScheme(TimeSchemeString,obj);
            
            %%%%%%%%% Determine newton structure %%%%%%%%%
            obj.newt = struct('maxIter',[],'eps',[],'quiet',[],'iter',[],'linesrch',[],'converge',[],'avgIter',[]);
            obj.newt.maxIter = extractInputROM(GNATtext,N,obj.id,'maxIter',10);
            obj.newt.eps     = extractInputROM(GNATtext,N,obj.id,'eps',1e-5);
            obj.newt.quiet   = extractInputROM(GNATtext,N,obj.id,'newtQuiet',false);
            obj.newt.iter    = zeros(1,obj.time.nstep);
            %Determine linesearch status and properties
            obj.newt.linesrch.status = extractInputROM(GNATtext,N,obj.id,'linesrch',false);
            obj.newt.linesrch.prop   = extractInputROM(GNATtext,N,obj.id,'linesrchProp',[1e-4,0.9,50,5,10,0.05,0.1,1.2]);
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
        
        %External Functions
        function  [] = createGNAT(obj,probobj,romobj,phiR,phiJ,~)
            %This function performs all of the offline GNAT computations.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - instance of GNAT class
            %probobj - Problem object
            %phiR    - basis of subspace residual is constrained to lie in
            %phiJ    - basis of subspace jacobian is constrained to lie in
            %phiY    - basis of subspace model is constrained to lie in
            %          (reduced order basis)
            %
            %Outputs:
            %--------
            %No outputs.  Properties are stored in GNAT handle class.
            %--------------------------------------------------------------
            
            
            if isempty(phiJ)
                phiJ = phiR;
            end
            
            
            if ~isempty(obj.nI) && isempty(obj.nSample) && isempty(obj.nGreed)
                %Compute sample indices
                computeSampleIndices(obj,probobj,phiR,phiJ);
                %Compute additional indices required in the reduced mesh in order
                %to compute the residual and jacobian at the sample nodes.
                determineAdditionalNodes(obj,probobj);
            elseif isempty(obj.nI) && ~isempty(obj.nSample) && ~isempty(obj.nGreed)
                %Compute the sample nodes based on a greedy node selection
                %algorithm
                computeSampleNodes(obj,probobj,phiR,phiJ);
                %Compute additional nodes required in the reduced mesh in order
                %to compute the residual and jacobian at the sample nodes.
                determineAdditionalNodes(obj,probobj);
                obj.nI = length(obj.sampleInd);
            end
            
            %Extract the phiYhat from the reduced order basis
            obj.phiYhat = romobj.phi(unique(obj.jrow),:);
            
            %Check to ensure that theoretical restrictions on nR, nJ, nI,
            %nY are met.
            if (obj.nR < obj.nY) || (obj.nJ < obj.nY)
                error(['To ensure uniqueness in the reconstructed gappy quantities,',...
                    'nR >= nY and nJ >= nY are required ']);
            end
            
            if (obj.nI < obj.nR) || (obj.nI < obj.nJ)
                error(['To ensure uniqueness in the least squares gappy reconstruction,',...
                    'nI >= nR and nI >= nJ are required ']);
            end
            
            %Mask SVD update data, if border_fast
            if strcmpi(romobj.basisUpdate,'border_fast') || strcmpi(romobj.basisUpdate,'border_fast_approx')
                obj.SingVals = romobj.SingVals;
                obj.svdUpdateData.wRef = romobj.svdUpdateData.wRef(unique(obj.jrow),1);
                obj.svdUpdateData.ColSumRSVtrunc = romobj.svdUpdateData.ColSumRSVtrunc;
                obj.svdUpdateData.normOneMinusVVt = romobj.svdUpdateData.normOneMinusVVt;
            end
            
            %Compute the online matrices A and B
            computeOnlineMatrices(obj,phiR(:,1:obj.nR),phiJ(:,1:obj.nJ),romobj.phi);
            
            %Initialize the state vector matrix (reduced coordinates)
            obj.sv = zeros(obj.nY,obj.time.nstep+1); %Reduced ic = 0
            
            %Determine phiYhatInd.  Let X be a full state vector, then
            %Xpartial = X(phiYhatInd,:), where Xpartial has no repeats.
            [~,obj.phiYhatInd,J] = unique(obj.jrow); J=J(:)';
            obj.reconstJhatInd = [];
            for i = 1:length(obj.irstart)-1
                obj.reconstJhatInd = [obj.reconstJhatInd, i + obj.nI*([J(obj.irstart(i):obj.irstart(i+1)-1)]-1)];%):J(obj.irstart(i+1)-1)]-1)];
            end
            obj.uniqueJROWind=J;
            
            %Initialize (sparse) JhatTemp
            obj.JhatTemp = spalloc(obj.nI,length(unique(obj.jrow)),length(obj.reconstJhatInd));
            obj.JhatFluxL = spalloc(obj.nI,length(unique(obj.jrow)),length(obj.reconstJhatInd));
            obj.JhatFluxR = spalloc(obj.nI,length(unique(obj.jrow)),length(obj.reconstJhatInd));
            
            %Set up the vector containing the indices of the partial
            %vectors that correspond to diagonal entries
            obj.jdiagHat = obj.uniqueJROWind(logical(obj.jdiag));
            
            %Store handle to problem in instance
            obj.probGNAT = probobj.createCopy4GNAT(obj);
            %Initialize the state vector matrix (reduced coordinates)
            obj.sv = zeros(obj.nY,obj.time.nstep+1); %Reduced ic = 0
            
            %Set up the indexing in the Time stepping scheme
            obj.TimeScheme.addGNAT(obj);
        end
        
        function  [] = executeModel(obj)
            %This function performs the time-stepping for the GNAT
            %simulation defined in this object.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - GNAT object
            %probobj - Problem object
            %
            %Outputs:
            %--------
            %There are no outputs.
            %--------------------------------------------------------------
            obj.numExecute = obj.numExecute+1;
            obj.killflag = false;
            
            %Ensure obj.sv = 0 initially!
            obj.sv = zeros(obj.nY,obj.time.maxsteps+1);
            
            nstep = obj.time.maxsteps; %Extract the number of time steps from the object
            
            %Make sure partialU and partialUprev are properly set whenever this
            %function is called
            obj.setUpartial2ic;
            
            %Set the current parameter (from prob object) to be used in
            %NAND to know which parameter values have been used to generate
            %the current state vectors
            obj.curr_param = obj.probGNAT(1).p;
            
            tGNATstart = tic;
            for i_t = 1:nstep %Loop over each time step to determine all state vectors
                if ~obj.time.quiet && rem(i_t,round(nstep*0.1)) == 0
                    %Generate output so user can see progress
                    fprintf('%4.0f of %4.0f Timesteps Complete (%2.0f%%)\n',i_t,nstep,100*i_t/nstep)
                end
                
                %Set current time iteration (so model can keep track)
                obj.cTimeIter = i_t;
                
                if ~isempty(obj.time.cfl)
                    if strcmpi(class(obj.probGNAT),'quasi1dEuler')
                        obj.time.dt = obj.time.cfl(i_t,obj)*obj.probGNAT.computeCFL_GNAT(obj.probGNAT.ic + obj.phiYhat*sum(obj.sv(:,1:i_t),2));
                    else
                        obj.time.dt = obj.time.cfl(i_t,obj)*obj.probGNAT.computeCFL_GNAT(obj.probGNAT.ic(obj.uniqueJROWind) + obj.phiYhat(obj.uniqueJROWind,:)*sum(obj.sv(:,1:i_t),2));
                    end
                    %obj.time.dt = obj.time.cfl(i_t,obj.sv(:,i_t))*obj.prob.computeCFL(obj.sv(:,i_t));
                    obj.TimeScheme.setDT(obj.time.dt);
                end
                
                %Use Newton's Method to solve nonlinear system of equations
                %and store: state vector, residual and jacobian snapshots,
                %and number of newton iterations required
                %                 keyboard
                if obj.solver==1
                    NewtonRaphson(obj);
                else
                    NewtonRaphsonconstr_new_method(obj);
                end
                
                if obj.killflag, warning('Encountered NaN at Time Step %i. Simulation Killed',obj.cTimeIter); return; end;
                
                %Check for steady convergence
                %   if (i_t==1)
                %     disp('Warning: put the next  line back on in localGNAT !!');
                % end
                %if 1==0
                %if obj.time.steadyconverge
                %    if norm(obj.partialU-obj.partialUprev) < obj.time.steadyconverge*norm(obj.partialU)
                %        obj.ontime = toc(tGNATstart); %Record simulation time
                %        obj.time.nstep=i_t;
                %        obj.time.T(2) = obj.time.T(1) + i_t*obj.time.dt;
                %        obj.sv(:,i_t+2:end)=[];
                %        obj.newt.avgIter = mean(obj.newt.iter);
                %        return;
                %    end
                %end
            end
            obj.ontime = toc(tGNATstart); %Record simulation time
            
            obj.newt.avgIter = mean(obj.newt.iter);
        end
        
        function  [] = associateFullProblem(obj,probobj)
            %This function associates the full problem (i.e. not the sample
            %mesh) with the GNAT object (i.e. it is stored in the prob
            %property)
            %---------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - GNAT object
            %probobj - Problem object
            %
            %Outputs:
            %--------
            %There are no outputs.
            %---------------------------------------------------------------
            
            obj.prob = probobj;
        end %Done
        
        function  [svF] = reconstructFullState(obj,romobj)
            %This function reconstructs the full state vector from the
            %reduced state vector matrix and the reduced basis (phi)
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - GNAT object
            %romobj  - ROM object
            %
            %Outputs:
            %--------
            %svF     - full state vector matrix
            %--------------------------------------------------------------
            
            svF = zeros(size(obj.prob.ic,1),obj.time.nstep+1);
            %First column is the initial condition.
            svF(:,1) = obj.prob.ic;
            for i = 2:obj.time.nstep+1
                %The newton method that was implemented for the GNAT
                %computes the current reduced state vector relative to the
                %previous state vector.  Therefore, we reconstructed the
                %full state by looping through the time steps and adding
                %the reconstructed state difference to the previous state
                %vector.
                svF(:,i) = svF(:,i-1) + romobj.phi*obj.sv(:,i);
            end
            
        end %Done
        
        %Internal Functions
        function  [] = setUpartial2ic(obj)
            %This function sets partialU and partialUprev to the ic.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - instance of GNAT class
            %
            %Outputs:
            %--------
            %No outputs.  partialU and partialUprev are overwritten in GNAT
            %handle class.
            %--------------------------------------------------------------
            
            obj.partialU = obj.probGNAT.ic;
            obj.partialUprev = obj.partialU;
        end %Done
        
        function  [] = computeOnlineMatrices(obj,phiR,phiJ,phiY)
            %This function computes the online matrices from Carlberg et.
            %al. 2011.  This is done in the offline phase.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj  - GNAT object
            %phiR - basis for the subspace the residual is constrained to
            %        lie in
            %phiJ - basis for the subspace the jacobian is constrained to
            %       lie in
            %
            %Outputs:
            %--------
            %There are no outputs for this problem.  A and B are stored in
            %the GNAT object.
            %--------------------------------------------------------------
            
            if lower(obj.Rtype(1))=='g'
                obj.A = phiY'*(phiJ*pinv(phiJ(obj.sampleInd,:)));
                obj.B = phiY'*(phiR*pinv(phiR(obj.sampleInd,:)));
            else
                obj.A = pinv(phiJ(obj.sampleInd,:));
                obj.B = (phiJ'*phiR)*pinv(phiR(obj.sampleInd,:));
                obj.Ebar = obj.A'*obj.B;
            end
            
        end
        
        function  [] = computeSampleNodes(obj,probobj,phiR,phiJ)
            %This function determines the sample nodes and corresponding
            %indices of the full state vector for use in the GNAT method.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - instance of GNAT class
            %probobj - Problem object
            %phiR    - basis of subspace residual is constrained to lie in
            %phiJ    - basis of subspace jacobian is constrained to lie in
            %
            %Outputs:
            %--------
            %There are no outputs.  sampleNodes and sampleInd are stored in
            %the GNAT handle class.
            %--------------------------------------------------------------
            
            obj.sampleNodes = [];
            obj.sampleInd = [];
            
            %Make sure that we have not specified too few sample nodes.
            if obj.nSample*probobj.ndim < obj.nGreed
                nSampleOLD = obj.nSample;
                obj.nSample = ceil(obj.nGreed/probobj.ndim);
                fprintf('Warning: not enough sample nodes! Increasing number of sample nodes from %d to %d\n',nSampleOLD,obj.nSample);
                clear nSampleOLD;
            end
            
            QQ = 0;
            P = min(obj.nSample,obj.nGreed);
            Q = floor(obj.nGreed/P)*ones(P,1);
            ind = [1:P] < mod(obj.nGreed,P);
            Q(ind) = Q(ind) + 1;
            s = ceil(obj.nGreed/obj.nSample);
            S = floor(obj.nSample*s/obj.nGreed)*ones(P,1);
            ind = [1:P] < mod(obj.nSample,obj.nGreed);
            S(ind) = S(ind)+1;
            
            %Lines 6 and 7 from Algorithm 4 of Carlberg's thesis
            R = phiR(:,1:max(Q));
            J = phiJ(:,1:max(Q));
            
            for p = 1:P
                QQ = QQ + Q(p);
                for s = 1:S(p)
                    %Determine the set of indices not already in our sample mesh
                    validNodes = setdiff(1:length(probobj.mesh.node),obj.sampleNodes');
                    temp = zeros(length(validNodes),1);
                    for i = 1:length(validNodes)
                        for q = 1:Q(p)
                            %Convert the node numbers to all their indices
                            %in the residual
                            NodeIndInRJ = probobj.node2ind(validNodes(i));
                            %Remove DBCs
                            if sum(isnan(NodeIndInRJ)) == length(NodeIndInRJ)
                                continue;
                            elseif sum(isnan(NodeIndInRJ)) > 0
                                NodeIndInRJ(isnan(NodeIndInRJ),:) = [];
                            end
                            %Add up the contribution of the residual and
                            %jacobian for each valid node
                            temp(i) = temp(i) + norm(R(NodeIndInRJ,q),2)^2 + norm(J(NodeIndInRJ,q),2)^2;
                        end
                    end
                    %Determine the node with the max. value
                    [~,ID] = max(temp);
                    n = validNodes(ID);
                    %Add the computed nodes and indices to the sample mesh
                    obj.sampleNodes = [obj.sampleNodes;n];
                    obj.sampleInd = [obj.sampleInd;probobj.node2ind(n)];
                end
                %Lines 13-16 of Algorithm 4 of Carlberg thesis
                for q = 1:Q(p)
                    a = phiR(obj.sampleInd,1:QQ)\phiR(obj.sampleInd,QQ+q);
                    b = phiJ( obj.sampleInd,1:QQ)\phiJ(obj.sampleInd,QQ+q);
                    
                    R(:,q) = phiR(:,QQ+q) - phiR(:,1:QQ)*a;
                    J(:,q) = phiJ(:,QQ+q) - phiJ(:,1:QQ)*b;
                end
            end
            
            %Sort the sample nodes and indices
            obj.sampleInd = unique(obj.sampleInd);
            obj.sampleNodes = unique(obj.sampleNodes);
            
            %Old Algorithm
            %             P = min(obj.nSample,obj.nGreed);
            %             Q = ceil(obj.nGreed/obj.nSample);
            %             S = floor(obj.nSample/obj.nGreed) + [ones(mod(obj.nSample,obj.nGreed),1);zeros(obj.nGreed - mod(obj.nSample,obj.nGreed),1)];
            %
            %             %Lines 6 and 7 from Algorithm 4 of Carlberg's thesis
            %             R = phiR(:,1:Q);
            %             J = phiJ(:,1:Q);
            %
            %             for p = 1:P
            %                 for s = 1:S
            %                     %Determine the set of indices not already in our sample mesh
            %                     validNodes = setdiff(1:length(probobj.mesh.node),obj.sampleNodes');
            %                     temp = zeros(length(validNodes),1);
            %                     for i = 1:length(validNodes)
            %                         for q = 1:Q
            %                             %Convert the node numbers to all their indices
            %                             %in the residual
            %                             NodeIndInRJ = probobj.node2ind(validNodes(i));
            %                             %Remove DBCs
            %                             if sum(isnan(NodeIndInRJ)) == length(NodeIndInRJ)
            %                                 continue;
            %                             elseif sum(isnan(NodeIndInRJ)) > 0
            %                                 NodeIndInRJ(isnan(NodeIndInRJ),:) = [];
            %                             end
            %                             %Add up the contribution of the residual and
            %                             %jacobian for each valid node
            %                             temp(i) = temp(i) + norm(R(NodeIndInRJ,q),2)^2 + norm(J(NodeIndInRJ,q),2)^2;
            %                         end
            %                     end
            %                     %Determine the node with the max. value
            %                     [~,id] = max(temp);
            %                     n = validNodes(id);
            %                     %Add the computed nodes and indices to the sample mesh
            %                     obj.sampleNodes = [obj.sampleNodes;n];
            %                     obj.sampleInd = [obj.sampleInd;probobj.node2ind(n)];
            %                 end
            %                 %Lines 13-16 of Algorithm 4 of Carlberg thesis
            %                 for q = 1:Q
            %                     a = phiR(obj.sampleInd,1:Q*p)\phiR(obj.sampleInd,Q*p+q);
            %                     b = phiJ( obj.sampleInd,1:Q*p)\phiJ(obj.sampleInd,Q*p+q);
            %
            %                     R(:,q) = phiR(:,Q*p+q) - phiR(:,1:Q*p)*a;
            %                     J(:,q) = phiJ(:,Q*p+q) - phiJ(:,1:Q*p)*b;
            %                 end
            %             end
            %
            %             %Sort the sample nodes and indices
            %             obj.sampleInd = unique(obj.sampleInd);
            %             obj.sampleNodes = unique(obj.sampleNodes);
        end
        
        function  [] = determineAdditionalNodes(obj,probobj)
            %This function determines all additional nodes whose state
            %needs to be evaluated in order to determine the residual and
            %jacobian on the sample mesh.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - GNAT object
            %probobj - Problem object
            %
            %Outputs:
            %--------
            %There are no outputs.  jdiag, jrow, and irstart are stored in
            %the GNAT handle class.
            %--------------------------------------------------------------
            
            obj.jdiag = [];
            obj.jrow = [];
            obj.irstart = 1;
            
            for i = 1:length(obj.sampleInd)
                %Find the nonzero entries in the appropriate row of Jstruct
                temp = find(probobj.Jstruct(obj.sampleInd(i,1),:) ~= 0);
                %Store these indices in jrow
                obj.jrow = [obj.jrow,temp];
                %Determine which of these values appended onto jrow
                %corresponds to the diagonal.
                obj.jdiag = [obj.jdiag,temp == obj.sampleInd(i)];
                
                %Store the index of the start of the next node
                obj.irstart = [obj.irstart, length(obj.jrow)+1];
            end
            %Make jdiag a column vector
            obj.jdiag = obj.jdiag(:);
        end
        
        function  [] = computeSampleIndices(obj,probobj,phiR,phiJ)
            %This function determines the sample indices and corresponding
            %indices of the full state vector for use in the GNAT method.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - instance of GNAT class
            %probobj - Problem object
            %phiR    - basis of subspace residual is constrained to lie in
            %phiJ    - basis of subspace jacobian is constrained to lie in
            %
            %Outputs:
            %--------
            %There are no outputs.  sampleNodes and sampleInd are stored in
            %the GNAT handle class.
            %--------------------------------------------------------------
            
            obj.sampleNodes = [];
            obj.sampleInd = [];
            
            % Algorithm 5 of Carlberg's 2010 Paper
            R = phiR(:,1);
            J = phiJ(:,1);
            
            nbarI = 0;
            m = 1;
            
            while nbarI < obj.nI
                %Determine the set of indices not already in our sample
                %indices
                validNodes = setdiff(1:length(R),obj.sampleInd');
                %Add up the contribution of the residual and
                %jacobian for each valid index
                temp = R(validNodes,1).^2 + J(validNodes,1).^2;
                %Determine the node with the max. value
                [~,ID] = max(temp);
                n = validNodes(ID);
                %Add the computed indices to the set
                obj.sampleInd = [obj.sampleInd;n];
                %Sort the sample nodes and indices
                obj.sampleInd = unique(obj.sampleInd);
                K = checkResEvaluatedForFree(obj,probobj,n,nbarI);
                obj.sampleInd = unique([obj.sampleInd;K(:)]);
                
                nbarI = nbarI + 1 + length(K);
                m = m+1;
                pR = min(m-1,obj.nR); pJ = min(m-1,obj.nJ);
                %Lines 10-11 of Algorithm 5 of Carlberg paper 2010
                a = phiR(obj.sampleInd,1:pR)\phiR(obj.sampleInd,m);
                b = phiJ( obj.sampleInd,1:pJ)\phiJ(obj.sampleInd,m);
                
                R = phiR(:,m) - phiR(:,1:pR)*a;
                J = phiJ(:,m) - phiJ(:,1:pJ)*b;
            end
        end
        
        function  [K] = checkResEvaluatedForFree(obj,probobj,greedInd,iIndex)
            %This function computes the additional Residual and Jacobian
            %entries that can be computed (for free) with the current mask.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj        - GNAT object
            %probobj    - Problem object
            %greedInd   - scalar indicating the "current" index to add to
            %             the sample index set.
            %iIndex     - current index iteration
            %
            %Outputs:
            %--------
            %K          - vector containing the set of indices that can be
            %             added for free
            %--------------------------------------------------------------
            %Code from David Amsallem 2010, with modifications made by
            %Matthew Zahr 2011 for integration with MORTestbed.
            
            switch lower(obj.addInd)
                case 'none'
                    K= [];
                case {'samevar', 'all'}
                    K = [];
                    switch obj.addInd
                        case 'samevar'
                            Jref = find(probobj.Jstruct(greedInd,:));
                        case 'all'
                            determineAdditionalNodes(obj,probobj);
                            Jref = obj.jrow;
                    end
                    IndexComp = setdiff(1:probobj.config.ndof,obj.sampleInd);
                    for i=1:length(IndexComp);
                        k = IndexComp(i);
                        J = find(probobj.Jstruct(k,:));
                        if (isempty(setdiff(J,Jref)))
                            K = [K k];
                        end
                        if (iIndex+1+length(K)>=obj.nI) % indices limit is attained
                            if (iIndex+1+length(K)>obj.nI)
                                K(:,end) = [];
                            end
                            break;
                        end
                    end
            end
        end
        
        function  [] = NewtonRaphson(obj)
            %This function solves systems of the form:  f(x) = 0 using
            %Newton-Raphson's method for the GNAT class
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - GNAT object
            %probobj - Problem object
            %
            %Outputs:
            %--------
            %There are no outputs.  The state vector, sv, is updated in the
            %GNAT handle class
            %--------------------------------------------------------------
            
            %Determine time iteration number plus 1
            itnump1 = obj.cTimeIter + 1;
            t = obj.time.T(1) + obj.time.dt*obj.cTimeIter;
            %Set initial guess as last time step (note that sv is the
            %partial state). Initial guess for the sv = 0 for GNAT, but we
            %do need to update the partial state initial guess.
            obj.partialUprev = obj.partialU;
            %Determine the residual and jacobian based on initial guess
            [Rhat,JVhat] = obj.TimeScheme.TimeIntNLFuncGNAT(obj.probGNAT,obj.partialU,obj.partialUprev,t);
            %Update local operators
            [RAbar,E] = updateLocOperators(obj,Rhat,JVhat);
            %Determine convergence criteria to use
            if obj.newt.converge == 1
                tol = obj.newt.eps*norm(E,2);
            elseif obj.newt.converge == 2
                tol = obj.newt.eps(1);
                tolIt = obj.newt.eps(2);
            else
                tol = max(obj.newt.eps(1)*norm(E,2),obj.newt.eps(2));
                tolIt = obj.newt.eps(3);
            end
            
            for i_N = 1:obj.newt.maxIter %loop until the maximum newton iterations are reached
                %Update the state vector (reduced) with search direction
                %(unit step length or line search)
                p = -RAbar\E;
                
                if sum(isnan(p)) > 0 || sum(abs(real(p)-p)>1e-6)
                    obj.killflag=true;
                    return;
                end
                
                %Linesearch if necessary
                if obj.newt.linesrch.status
                    alpha = linesearch(@(alpha) obj.TimeScheme.TimeIntNLFuncGNATDir(obj.probGNAT,obj.partialU,obj.partialUprev,t,alpha,p),obj.newt.linesrch.prop);
                    p = alpha*p;
                end
                
                obj.sv(:,itnump1) = obj.sv(:,itnump1) + p;
                obj.partialU = obj.partialUprev + obj.phiYhat*obj.sv(:,itnump1);
                %Compute residual and jacobian with updated vector for next iteration
                [Rhat,JVhat] = obj.TimeScheme.TimeIntNLFuncGNAT(obj.probGNAT,obj.partialU,obj.partialUprev,t);
                %[Rhat,JVhat] = BackwardEulerNLFuncGNAT(obj);
                %Solve for the search direction
                [RAbar,E] = updateLocOperators(obj,Rhat,JVhat);
                
                %Stop iterating once residual is small enough
                if (norm(E)<tol)
                    %If were are using convergence criteria 2 or 3, compute the
                    %difference between two iterates and compute to absolute
                    %tolerance
                    if obj.newt.converge ~= 1
                        itDiffNorm = norm(p,2);
                        if itDiffNorm<tolIt
                            break;
                        end
                    else
                        break;
                    end
                end
            end
            
            %Store the number of newton iterations at this time step
            obj.newt.iter(obj.cTimeIter+1) = i_N;
            
            %If the maxiumum newton iterations are reached, warn the user
            if ~obj.newt.quiet && (i_N==obj.newt.maxIter)
                disp(['*** Warning; Newton solver reached max number of iterations before convergence : res = ',...
                    num2str(norm(E,2)/tol),'  eps = ',num2str(obj.newt.eps), '***']);
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %             function  [] = NewtonRaphsonconstr(obj)
        %
        %                 itnump1 = obj.cTimeIter + 1;
        %                 t = obj.time.T(1) + obj.time.dt*obj.cTimeIter;
        %                 obj.partialUprev = obj.partialU;
        %                 %Determine the residual and jacobian based on initial guess
        %                 [Rhat,JVhat] = obj.TimeScheme.TimeIntNLFuncGNAT(obj.probGNAT,obj.partialU,obj.partialUprev,t);
        %                 %Update local operators
        %                 [RAbar,E] = updateLocOperators(obj,Rhat,JVhat);
        %
        % 		keyboard
        %                 for i_N = 1:obj.newt.maxIter %loop until the maximum newton iterations are reached
        %
        %                     %obj.sv(:,itnump1)=obj.sv(:,obj.cTimeIter)+obj.phi*w_guess;
        %                     %[Res,Jres]=obj.TimeScheme.TimeIntNLFunc(obj.prob,obj.sv(:,itnump1),obj.sv(:,itnump1-1),t);
        %                     %Df=Jres*obj.phi;
        %                     U = reshape(obj.partialU,3,obj.probGNAT.nVolSamp);
        %                     startInd = 1 + obj.probGNAT.ind1;
        %                     endInd   = obj.probGNAT.nVolSamp - obj.probGNAT.indN;
        %                     [rho,u,P,c] = obj.probGNAT.conservativeToPrimitive(U(:,startInd:endInd));
        %                     e=U(3,startInd:endInd);
        %                     if obj.probGNAT.ind1
        %                         rho = [U(1,1),rho]; %Density
        %                         u   = [U(2,1),u]; %Velocity
        %                         P   = [U(3,1),P]; %Pressure
        %                         c   = [sqrt(obj.probGNAT.gamma*P(1)/rho(1)),c]; %Speed of sound
        %                         e   = [P(1)/(obj.probGNAT.gamma-1)+rho(1)*u(1)^2/2,e]; %Energy
        %                     end
        %                     if obj.probGNAT.indN
        %                         rho = [rho,U(1,end)]; %Density
        %                         u   = [u,U(2,end)]; %Velocity
        %                         P   = [P,U(3,end)]; %Pressure
        %                         c   = [c,sqrt(obj.probGNAT.gamma*P(end)/rho(end))]; %Speed of sound
        %                         e   = [e,P(end)/(obj.probGNAT.gamma-1)+rho(end)*u(end)^2/2]; %Energy
        %                     end
        %
        %                     dc_cons = [0.5*obj.probGNAT.gamma./(c.*rho).*(0.5*(obj.probGNAT.gamma-1)*u.*u - P./rho);...
        %                         -0.5*obj.probGNAT.gamma*(obj.probGNAT.gamma-1)*u./(rho.*c);...
        %                         0.5*obj.probGNAT.gamma*(obj.probGNAT.gamma-1)./(rho.*c)]';
        %
        %                     %[rho, u, P,c,e,~,dc]=obj.prob.getVariables(obj.sv(:,itnump1));
        %                     [Q,dQ]=obj.probGNAT.forceTermGNAT(u,P);
        %                     [roeF, droeF]=obj.probGNAT.roeFluxGNAT(rho,u,P,c,e,dc_cons);
        %                     if obj.probGNAT.ind1
        %                         dUdV=[1,0,0;u(1),rho(1),0;0.5*u(1)*u(1),rho(1)*u(1),1/(obj.probGNAT.gamma-1)];
        %                         droeF(:,1:3,1)=droeF(:,1:3,1)*dUdV;
        %                     end
        %                     if obj.probGNAT.indN
        %                         dUdV=[1,0,0;u(end),rho(end),0;0.5*u(end)*u(end),rho(end)*u(end),1/(obj.probGNAT.gamma-1)];
        %                         droeF(:,4:6,end)=droeF(:,4:6,end)*dUdV;
        %                     end
        %
        %                     Rright= bsxfun(@times,roeF(:,obj.probGNAT.iarrayFaceI+1),obj.probGNAT.S(obj.probGNAT.iarrayFace(1+obj.probGNAT.ind1:end-obj.probGNAT.indN)+1));
        %                     Rleft= -bsxfun(@times,roeF(:,obj.probGNAT.iarrayFaceI),obj.probGNAT.S(obj.probGNAT.iarrayFace(1+obj.probGNAT.ind1:end-obj.probGNAT.indN)));
        %
        % 		    PhiFr1=obj.problem.phiFright(1:3:end,:);
        % 		    PhiFr2=obj.problem.phiFright(2:3:end,:);
        % 		    PhiFr3=obj.problem.phiFright(3:3:end,:);
        %
        %                     PhiFl1=obj.problem.phiFleft(1:3:end,:);
        %                     PhiFl2=obj.problem.phiFleft(2:3:end,:);
        %                     PhiFl3=obj.problem.phiFleft(3:3:end,:);
        %
        %                     PhiQ1=obj.problem.phiQ(1:3:end,:);
        %                     PhiQ2=obj.problem.phiQ(2:3:end,:);
        %                     PhiQ3=obj.problem.phiQ(3:3:end,:);
        %
        %
        %                     if obj.probGNAT.ind1
        %                         kk=obj.probGNAT.iarrayFaceI;
        %                         alphaRr(:,1)=PhiFr1(obj.probGNAT.iarrayFaceI,:)\Rright(1,:)';
        %                         alphaRr(:,2)=PhiFr2(obj.probGNAT.iarrayFaceI,:)\Rright(2,:)';
        %                         alphaRr(:,3)=PhiFr3(obj.probGNAT.iarrayFaceI,:)\Rright(3,:)';
        %
        %                         alphaQ(:,1)=PhiQ1(obj.probGNAT.iarrayFaceI,:)\Q(1,:)';
        %                         alphaQ(:,2)=PhiQ2(obj.probGNAT.iarrayFaceI,:)\Q(2,:)';
        %                         alphaQ(:,3)=PhiQ3(obj.probGNAT.iarrayFaceI,:)\Q(3,:)';
        %
        %
        %                         %keyboard
        %                         g(1:3,:)=[(obj.probGNAT.SVol(2:end-1).*obj.probGNAT.dx(2:end-1))*Phi1(2:end-1,:)*w_guess;
        %                             (obj.probGNAT.SVol(2:end-1).*obj.probGNAT.dx(2:end-1))*Phi2(2:end-1,:)*w_guess;
        %                             (obj.probGNAT.SVol(2:end-1).*obj.probGNAT.dx(2:end-1))*Phi3(2:end-1,:)*w_guess ]+...
        %                             obj.time.dt*([PhiFr1(end,:)*alphaRr(:,1);
        %                             		  PhiFr2(end,:)*alphaRr(:,2);
        %                             		  PhiFr3(end,:)*alphaRr(:,3)]+ Rleft(:,1)-...
        %                             [(obj.probvGNAT.SVol(2:end-1).*obj.probGNAT.dx(2:end-1))*PhiQ1(2:end-1,:)*alphaQ(:,1);
        %                             (obj.probGNAT.SVol(2:end-1).*obj.probGNAT.dx(2:end-1))*PhiQ2(2:end-1,:)*alphaQ(:,2);
        %                             (obj.probGNAT.SVol(2:end-1).*obj.probGNAT.dx(2:end-1))*PhiQ3(2:end-1,:)*alphaQ(:,3)]);
        %
        %                     end
        %
        %
        %                     if obj.probGNAT.indN
        %                         kk=obj.probGNAT.iarrayFaceI;
        %                         alphaRl(:,1)=PhiFl1(obj.probGNAT.iarrayFaceI,:)\Rleft(1,:)';
        %                         alphaRl(:,2)=PhiFl2(obj.probGNAT.iarrayFaceI,:)\Rleft(2,:)';
        %                         alphaRl(:,3)=PhiFl3(obj.probGNAT.iarrayFaceI,:)\Rleft(3,:)';
        %
        %                         alphaQ(:,1)=PhiQ1(obj.probGNAT.iarrayFaceI,:)\Q(1,:)';
        %                         alphaQ(:,2)=PhiQ2(obj.probGNAT.iarrayFaceI,:)\Q(2,:)';
        %                         alphaQ(:,3)=PhiQ3(obj.probGNAT.iarrayFaceI,:)\Q(3,:)';
        %
        %
        %                         %keyboard
        %                         g(1:3,:)=[(obj.probGNAT.SVol(2:end-1).*obj.probGNAT.dx(2:end-1))*Phi1(2:end-1,:)*w_guess;
        %                             (obj.probGNAT.SVol(2:end-1).*obj.probGNAT.dx(2:end-1))*Phi2(2:end-1,:)*w_guess;
        %                             (obj.probGNAT.SVol(2:end-1).*obj.probGNAT.dx(2:end-1))*Phi3(2:end-1,:)*w_guess ]+...
        %                             obj.time.dt*([PhiFl1(1,:)*alphaRl(:,1);
        %                                           PhiFl2(1,:)*alphaRl(:,2);
        %                                           PhiFl3(1,:)*alphaRl(:,3)]+ Rright(:,1)-...
        %                             [(obj.probvGNAT.SVol(2:end-1).*obj.probGNAT.dx(2:end-1))*PhiQ1(2:end-1,:)*alphaQ(:,1);
        %                             (obj.probGNAT.SVol(2:end-1).*obj.probGNAT.dx(2:end-1))*PhiQ2(2:end-1,:)*alphaQ(:,2);
        %                             (obj.probGNAT.SVol(2:end-1).*obj.probGNAT.dx(2:end-1))*PhiQ3(2:end-1,:)*alphaQ(:,3)]);
        %
        %                     end
        %                     p = -RAbar\E;
        %
        %                     obj.sv(:,itnump1) = obj.sv(:,itnump1) + p;
        %                     obj.partialU = obj.partialUprev + obj.phiYhat*obj.sv(:,itnump1);
        %                     %Compute residual and jacobian with updated vector for next iteration
        %                     [Rhat,JVhat] = obj.TimeScheme.TimeIntNLFuncGNAT(obj.probGNAT,obj.partialU,obj.partialUprev,t);
        %                     %[Rhat,JVhat] = BackwardEulerNLFuncGNAT(obj);
        %                     %Solve for the search direction
        %                     [RAbar,E] = updateLocOperators(obj,Rhat,JVhat);
        %
        %                     %Stop iterating once residual is small enough
        %                     if (norm(E)<tol)
        %
        %                     end
        %                 end
        %                 %Store the number of newton iterations at this time step
        %                 obj.newt.iter(obj.cTimeIter+1) = i_N;
        %
        %                 %If the maxiumum newton iterations are reached, warn the user
        %                 if ~obj.newt.quiet && (i_N==obj.newt.maxIter)
        %                     disp(['*** Warning; Newton solver reached max number of iterations before convergence : res = ',...
        %                         num2str(norm(E,2)/tol),'  eps = ',num2str(obj.newt.eps), '***']);
        %                 end
        %             end
        
        
        function  [] = NewtonRaphsonconstr(obj)
            
            itnump1 = obj.cTimeIter + 1;
            t = obj.time.T(1) + obj.time.dt*obj.cTimeIter;
            obj.partialUprev = obj.partialU;
            %Determine the residual and jacobian based on initial guess
            [Rhat,JVhat] = obj.TimeScheme.TimeIntNLFuncGNAT(obj.probGNAT,obj.partialU,obj.partialUprev,t);
%             keyboard
            Df=JVhat;
            w_guess=obj.phiYhat\(obj.partialU-obj.partialUprev);
            
%              keyboard
            for i_N = 1:obj.newt.maxIter %loop until the maximum newton iterations are reached
   
                U = reshape(obj.partialU,3,obj.probGNAT.nVolSamp);
                startInd = 1 + obj.probGNAT.ind1;
                endInd   = obj.probGNAT.nVolSamp - obj.probGNAT.indN;
                [rho,u,P,c] = obj.probGNAT.conservativeToPrimitive(U(:,startInd:endInd));
                e=U(3,startInd:endInd);
                if obj.probGNAT.ind1
                    rho = [U(1,1),rho]; %Density
                    u   = [U(2,1),u]; %Velocity
                    P   = [U(3,1),P]; %Pressure
                    c   = [sqrt(obj.probGNAT.gamma*P(1)/rho(1)),c]; %Speed of sound
                    e   = [P(1)/(obj.probGNAT.gamma-1)+rho(1)*u(1)^2/2,e]; %Energy
                end
                if obj.probGNAT.indN
                    rho = [rho,U(1,end)]; %Density
                    u   = [u,U(2,end)]; %Velocity
                    P   = [P,U(3,end)]; %Pressure
                    c   = [c,sqrt(obj.probGNAT.gamma*P(end)/rho(end))]; %Speed of sound
                    e   = [e,P(end)/(obj.probGNAT.gamma-1)+rho(end)*u(end)^2/2]; %Energy
                end
                
                dc_cons = [0.5*obj.probGNAT.gamma./(c.*rho).*(0.5*(obj.probGNAT.gamma-1)*u.*u - P./rho);...
                    -0.5*obj.probGNAT.gamma*(obj.probGNAT.gamma-1)*u./(rho.*c);...
                    0.5*obj.probGNAT.gamma*(obj.probGNAT.gamma-1)./(rho.*c)]';
                
                %[rho, u, P,c,e,~,dc]=obj.prob.getVariables(obj.sv(:,itnump1));
                [Q,dQ]=obj.probGNAT.forceTermGNAT(u,P);
                %keyboard
                [roeF, droeF]=obj.probGNAT.roeFluxGNAT(rho,u,P,c,e,dc_cons);
                if obj.probGNAT.ind1
                    dUdV=[1,0,0;u(1),rho(1),0;0.5*u(1)*u(1),rho(1)*u(1),1/(obj.probGNAT.gamma-1)];
                    droeF(:,1:3,1)=droeF(:,1:3,1)*dUdV;
                end
                if obj.probGNAT.indN
                    dUdV=[1,0,0;u(end),rho(end),0;0.5*u(end)*u(end),rho(end)*u(end),1/(obj.probGNAT.gamma-1)];
                    droeF(:,4:6,end)=droeF(:,4:6,end)*dUdV;
                end
                %keyboard 
                Rright= bsxfun(@times,roeF(:,obj.probGNAT.iarrayFaceI+1),obj.probGNAT.S(obj.probGNAT.iarrayFace(1+obj.probGNAT.ind1:end-obj.probGNAT.indN)+1));
                Rleft =-bsxfun(@times,roeF(:,obj.probGNAT.iarrayFaceI),obj.probGNAT.S(obj.probGNAT.iarrayFace(1+obj.probGNAT.ind1:end-obj.probGNAT.indN)));
               
% indFR=[];
%for ii=1:298
%            for jj=1:9
%                if Right(:,jj)==Fresh(:,ii)
%                    indFr=[indFr,ii];
%                end
%            end
%        end 
                J2R=zeros(3*9*(obj.probGNAT.nVolMask-obj.probGNAT.ind1-obj.probGNAT.indN),1);
                J2L=zeros(3*9*(obj.probGNAT.nVolMask-obj.probGNAT.ind1-obj.probGNAT.indN),1);
                
%                  keyboard
                for k = 1:obj.probGNAT.nVolMask-obj.probGNAT.ind1-obj.probGNAT.indN
                   
%                    JnewL(3*(k-1)+1:3*k,1:6)=[-obj.probGNAT.S(obj.probGNAT.iarrayFace(k+obj.probGNAT.ind1))*droeF(:,1:3,obj.probGNAT.iarrayFaceI(k)),...
%                        - obj.probGNAT.S(obj.probGNAT.iarrayFace(k+obj.probGNAT.ind1))*droeF(:,4:6,obj.probGNAT.iarrayFaceI(k))];
%                    JnewR(3*(k-1)+1:3*k,1:6)=[ obj.probGNAT.S(obj.probGNAT.iarrayFace(k+obj.probGNAT.ind1)+1)*droeF(:,1:3,obj.probGNAT.iarrayFaceI(k)+1),...
%                        obj.probGNAT.S(obj.probGNAT.iarrayFace(k+obj.probGNAT.ind1)+1)*droeF(:,4:6,obj.probGNAT.iarrayFaceI(k)+1)];
                    
%                     tempL = [-obj.probGNAT.S(obj.probGNAT.iarrayFace(k+obj.probGNAT.ind1))*droeF(:,1:3,obj.probGNAT.iarrayFaceI(k)),...
%                         - obj.probGNAT.S(obj.probGNAT.iarrayFace(k+obj.probGNAT.ind1))*droeF(:,4:6,obj.probGNAT.iarrayFaceI(k)),...
%                         zeros(3)]';
    
                    tempL = -obj.probGNAT.S(obj.probGNAT.iarrayFace(k+obj.probGNAT.ind1))*[droeF(:,:,obj.probGNAT.iarrayFaceI(k)), zeros(3)]';
                    J2L(27*(k-1)+1:27*k,:)=tempL(:);
                    
                    tempR =  obj.probGNAT.S(obj.probGNAT.iarrayFace(k+obj.probGNAT.ind1)+1)*[zeros(3),droeF(:,:,obj.probGNAT.iarrayFaceI(k)+1)]';
                    J2R(27*(k-1)+1:27*k,:)=tempR(:);
                    
                end
                
                
                Phi1=obj.problem.phi(1:3:end,:); %basis for the first conserved quantity rho
                Phi2=obj.problem.phi(2:3:end,:); %basis for the second conserved quantity rho*u
                Phi3=obj.problem.phi(3:3:end,:); %basis for the thirs conserved quantity e
                
                PhiFr1=obj.problem.phiFright(1:3:end,:);
                PhiFr2=obj.problem.phiFright(2:3:end,:);
                PhiFr3=obj.problem.phiFright(3:3:end,:);
%                 
%                 PhiFl1=obj.problem.phiFleft(1:3:end,:);
%                 PhiFl2=obj.problem.phiFleft(2:3:end,:);
%                 PhiFl3=obj.problem.phiFleft(3:3:end,:);
                
                PhiQ1=obj.problem.phiQ(1:3:end,:);
                PhiQ2=obj.problem.phiQ(2:3:end,:);
                PhiQ3=obj.problem.phiQ(3:3:end,:);
                %dQ has only the derivative of Q2, since others are 0
%                  keyboard
                derQ=PhiQ2*pinv(PhiQ2(obj.sampleNodes,:))*squeeze(dQ(1,:,:))'; %nVol x 3 matrix, 3 corresponds to derivative wrt rho, rho*u and e                
                derivQ=zeros(1, length(w_guess));
                
                for i=2:obj.problem.prob.nVol-1
                    derivQ=derivQ+derQ(i,:)*[Phi1(i,:);Phi2(i,:);Phi3(i,:)]*(obj.problem.prob.SVol(i).*obj.problem.prob.dx(i));
                end
                
                if obj.probGNAT.ind1
%                    keyboard
                    g(1:3,:)=[(obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi1(2:end-1,:)*w_guess;
                        (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi2(2:end-1,:)*w_guess;
                        (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi3(2:end-1,:)*w_guess ]+...
                        obj.time.dt*([PhiFr1(end-1,:)*pinv(PhiFr1(setdiff(obj.sampleNodes,obj.probGNAT.ind1)-1,:))*Rright(1,:)';
                        PhiFr2(end-1,:)*pinv(PhiFr2(setdiff(obj.sampleNodes,obj.probGNAT.ind1)-1,:))*Rright(2,:)'; 
                        PhiFr3(end-1,:)*pinv(PhiFr3(setdiff(obj.sampleNodes,obj.probGNAT.ind1)-1,:))*Rright(3,:)']+ Rleft(:,1)-...
                        [(obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*PhiQ1(2:end-1,:)*pinv(PhiQ1(obj.sampleNodes,:))*Q(1,:)';
                        (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*PhiQ2(2:end-1,:)*pinv(PhiQ2(obj.sampleNodes,:))*Q(2,:)';
                        (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*PhiQ3(2:end-1,:)*pinv(PhiQ3(obj.sampleNodes,:))*Q(3,:)']);
                    J2L = [zeros(18,1);J2L];
                    J2R = [zeros(18,1);J2R];
                    
%                      keyboard
                    obj.JhatFluxL(obj.reconstJhatInd) = J2L; % 30xlength(unique(obj.jrow))
                    obj.JhatFluxR(obj.reconstJhatInd) = J2R;
                    
                    %                     index=[];
                    %                     for jj=1:length(obj.probGNAT.iarrayFaceI)
                    %                         index=[index, 3*(obj.probGNAT.iarrayFaceI(jj)-1)+1:3*obj.probGNAT.iarrayFaceI(jj)];
                    %                     end
                    %
                    %                     keyboard
                    %                     JRfull=obj.problem.phiFright*pinv(obj.problem.phiFright(index,:))*JnewR;
                    %                     JLfull=obj.problem.phiFleft*pinv(obj.problem.phiFleft(index,:))*JnewL;
                     keyboard
%                     Dg=[(obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi1(2:end-1,:);
%                         (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi2(2:end-1,:);
%                         (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi3(2:end-1,:) ]+...
%                         obj.time.dt*(JLfull(1:3,1:3)*obj.problem.phi(1:3,:) +JLfull(1:3,4:6)*obj.problem.phi(4:6,:)+ ...
%                         JRfull(end-2:end,1:3)*obj.problem.phi(end-5:end-3,:)+ JRfull(end-2:end,4:6)*obj.problem.phi(end-2:end,:)- [zeros(1,5);derivQ;zeros(1,5)]);
%                  keyboard
                    Dg=[(obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi1(2:end-1,:);
                        (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi2(2:end-1,:);
                        (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi3(2:end-1,:) ]+...
                        obj.time.dt*(obj.JhatFluxL(4:6,:)*obj.problem.phi(unique(obj.jrow),:) + ...
                               obj.JhatFluxR(end-2:end,:)*obj.problem.phi(unique(obj.jrow),:)- [zeros(1,length(w_guess));derivQ;zeros(1,length(w_guess))]);
                
                
                end
                
                
%                 if obj.probGNAT.indN
%                     
%                     %keyboard
%                     g(1:3,:)=[(obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi1(2:end-1,:)*w_guess;
%                         (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi2(2:end-1,:)*w_guess;
%                         (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi3(2:end-1,:)*w_guess ]+...
%                         obj.time.dt*([PhiFl1(1,:)*pinv(PhiFl1(setdiff(obj.sampleNodes, obj.probGNAT.indN),:))*Rleft(1,:)';
%                         PhiFl2(1,:)*pinv(PhiFl2(setdiff(obj.sampleNodes, obj.probGNAT.indN),:))*Rleft(2,:)';
%                         PhiFl3(1,:)*pinv(PhiFl3(setdiff(obj.sampleNodes, obj.probGNAT.indN),:))*Rleft(3,:)']+ Rright(:,1)-...
%                          (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*PhiQ2(2:end-1,:)*pinv(PhiQ2(obj.problem.prob.sampleNodes,:))*Q(2,:)';
%                         (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*PhiQ3(2:end-1,:)*pinv(PhiQ3(obj.problem.prob.sampleNodes,:))*Q(3,:)']);
%                     
%                 end
%                 keyboard
                H=Df'*Df; h=Df'*Rhat;
                P=H\Dg';
                x=H\h;
                S=-inv(Dg*P);
                Qzim=P*S;
                del_w=Qzim*g-(x+Qzim*Dg*x);
                w_guess=w_guess+del_w;
%keyboard
                if norm(del_w,2)<10^(-6)
                    break;
                end
                obj.sv(:,itnump1)=w_guess;
%                 obj.sv(:,itnump1) = obj.sv(:,itnump1) + w_guess;
                obj.partialU = obj.partialUprev + obj.phiYhat*obj.sv(:,itnump1);
                %Compute residual and jacobian with updated vector for next iteration
                [Rhat,JVhat] = obj.TimeScheme.TimeIntNLFuncGNAT(obj.probGNAT,obj.partialU,obj.partialUprev,t);
                
            end
             disp(['norm of the constraint   ', num2str(norm(g))])
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function  [] = NewtonRaphsonconstr_new(obj)
            
            itnump1 = obj.cTimeIter + 1;
            t = obj.time.T(1) + obj.time.dt*obj.cTimeIter;
            obj.partialUprev = obj.partialU;
            %Determine the residual and jacobian based on initial guess
            [Rhat,JVhat] = obj.TimeScheme.TimeIntNLFuncGNAT(obj.probGNAT,obj.partialU,obj.partialUprev,t);
%             keyboard
            Df=JVhat;
            w_guess=obj.phiYhat\(obj.partialU-obj.partialUprev);
            
%               keyboard
            for i_N = 1:obj.newt.maxIter %loop until the maximum newton iterations are reached
   
                U = reshape(obj.partialU,3,obj.probGNAT.nVolSamp);
                startInd = 1 + obj.probGNAT.ind1;
                endInd   = obj.probGNAT.nVolSamp - obj.probGNAT.indN;
                [rho,u,P,c] = obj.probGNAT.conservativeToPrimitive(U(:,startInd:endInd));
                e=U(3,startInd:endInd);
                if obj.probGNAT.ind1
                    rho = [U(1,1),rho]; %Density
                    u   = [U(2,1),u]; %Velocity
                    P   = [U(3,1),P]; %Pressure
                    c   = [sqrt(obj.probGNAT.gamma*P(1)/rho(1)),c]; %Speed of sound
                    e   = [P(1)/(obj.probGNAT.gamma-1)+rho(1)*u(1)^2/2,e]; %Energy
                end
                if obj.probGNAT.indN
                    rho = [rho,U(1,end)]; %Density
                    u   = [u,U(2,end)]; %Velocity
                    P   = [P,U(3,end)]; %Pressure
                    c   = [c,sqrt(obj.probGNAT.gamma*P(end)/rho(end))]; %Speed of sound
                    e   = [e,P(end)/(obj.probGNAT.gamma-1)+rho(end)*u(end)^2/2]; %Energy
                end
                
                dc_cons = [0.5*obj.probGNAT.gamma./(c.*rho).*(0.5*(obj.probGNAT.gamma-1)*u.*u - P./rho);...
                    -0.5*obj.probGNAT.gamma*(obj.probGNAT.gamma-1)*u./(rho.*c);...
                    0.5*obj.probGNAT.gamma*(obj.probGNAT.gamma-1)./(rho.*c)]';
                
                %[rho, u, P,c,e,~,dc]=obj.prob.getVariables(obj.sv(:,itnump1));
                [Q,dQ]=obj.probGNAT.forceTermGNAT(u,P);
                %keyboard
                [roeF, droeF]=obj.probGNAT.roeFluxGNAT(rho,u,P,c,e,dc_cons);
                if obj.probGNAT.ind1
                    dUdV=[1,0,0;u(1),rho(1),0;0.5*u(1)*u(1),rho(1)*u(1),1/(obj.probGNAT.gamma-1)];
                    droeF(:,1:3,1)=droeF(:,1:3,1)*dUdV;
                end
                if obj.probGNAT.indN
                    dUdV=[1,0,0;u(end),rho(end),0;0.5*u(end)*u(end),rho(end)*u(end),1/(obj.probGNAT.gamma-1)];
                    droeF(:,4:6,end)=droeF(:,4:6,end)*dUdV;
                end
                
                Rright= bsxfun(@times,roeF(:,obj.probGNAT.iarrayFaceI+1),obj.probGNAT.S(obj.probGNAT.iarrayFace(1+obj.probGNAT.ind1:end-obj.probGNAT.indN)+1));
                Rleft =-bsxfun(@times,roeF(:,obj.probGNAT.iarrayFaceI),obj.probGNAT.S(obj.probGNAT.iarrayFace(1+obj.probGNAT.ind1:end-obj.probGNAT.indN)));
               
                for k = 1:obj.probGNAT.nVolMask-obj.probGNAT.ind1-obj.probGNAT.indN
                
                    J2Lnew(3*k-2:3*k,:) =-obj.probGNAT.S(obj.probGNAT.iarrayFace(k+obj.probGNAT.ind1)) * droeF(:,:,obj.probGNAT.iarrayFaceI(k));
                    J2Rnew(3*k-2:3*k,:) = obj.probGNAT.S(obj.probGNAT.iarrayFace(k+obj.probGNAT.ind1)+1)*droeF(:,:,obj.probGNAT.iarrayFaceI(k)+1);
                         
                end
                
%                 keyboard 
                Phi1=obj.problem.phi(1:3:end,:); %basis for the first conserved quantity rho
                Phi2=obj.problem.phi(2:3:end,:); %basis for the second conserved quantity rho*u
                Phi3=obj.problem.phi(3:3:end,:); %basis for the thirs conserved quantity e
                PhiQ1=obj.problem.phiQ(1:3:end,:);
                PhiQ2=obj.problem.phiQ(2:3:end,:);
                PhiQ3=obj.problem.phiQ(3:3:end,:);
                %dQ has only the derivative of Q2, since others are 0
%               keyboard
                derQ=PhiQ2*pinv(PhiQ2(obj.sampleNodes,:))*squeeze(dQ(1,:,:))'; %nVol x 3 matrix, 3 corresponds to derivative wrt rho, rho*u and e                
                derivQ=zeros(1, length(w_guess));
                
                for i=2:obj.problem.prob.nVol-1
                    derivQ=derivQ+derQ(i,:)*[Phi1(i,:);Phi2(i,:);Phi3(i,:)]*(obj.problem.prob.SVol(i).*obj.problem.prob.dx(i));
                end
                
%                 keyboard
%                 
%                 load romroeF
%                 reshF=reshape(romroeF(:,1),3, 299);
%                 Index=[];
%                 for jj=1:299
%                     for kk=1:size(roeF,2)
%                         if abs(roeF(:,kk)-reshF(:,jj))<10^(-10)
%                             Index=[Index,jj];
%                         end
%                     end
%                 end
                
                if obj.probGNAT.ind1

                    fluxLeft = obj.problem.phiFleft *pinv(obj.problem.phiFleft(setdiff(obj.sampleInd,[1,2,3])-3,:))*Rleft(:);
                    fluxRight= obj.problem.phiFright*pinv(obj.problem.phiFright(setdiff(obj.sampleInd,[1,2,3])-3,:))*Rright(:);
%                     J2L = [zeros(18,1);J2L];
%                     J2R = [zeros(18,1);J2R];
%                     keyboard
                    JRnew =zeros(3*(obj.probGNAT.nVolMask-obj.probGNAT.ind1), 3*obj.problem.prob.nVol);
                    JLnew =zeros(3*(obj.probGNAT.nVolMask-obj.probGNAT.ind1), 3*obj.problem.prob.nVol);
                    for i=1:obj.probGNAT.nVolMask-obj.probGNAT.ind1
                        ii=obj.sampleNodes(i+1);
                        JRnew(3*i-2:3*i,3*ii-2:3*(ii+1))=J2Rnew(3*i-2:3*i,:);
                        JLnew(3*i-2:3*i,3*ii-5:3*ii)=J2Lnew(3*i-2:3*i,:);
                    end
                    
%                     keyboard
%                     obj.JhatFluxL(obj.reconstJhatInd) = J2L; % 30xlength(unique(obj.jrow))
%                     obj.JhatFluxR(obj.reconstJhatInd) = J2R;

%                     JRnew(end-2:end,:)=[]; JRnew(1:3,:)=[];
%                     JLnew(end-2:end,:)=[]; JLnew(1:3,:)=[];
                    
                    keyboard
                    g(1:3,:)=[(obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi1(2:end-1,:)*w_guess;
                        (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi2(2:end-1,:)*w_guess;
                        (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi3(2:end-1,:)*w_guess ]+...
                        obj.time.dt*(fluxLeft(1:3) + fluxRight(end-2:end)-...
                        [(obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*PhiQ1(2:end-1,:)*pinv(PhiQ1(obj.sampleNodes,:))*Q(1,:)';
                        (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*PhiQ2(2:end-1,:)*pinv(PhiQ2(obj.sampleNodes,:))*Q(2,:)';
                        (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*PhiQ3(2:end-1,:)*pinv(PhiQ3(obj.sampleNodes,:))*Q(3,:)']);
                    
                    Jleft =obj.problem.phiFleft*pinv(obj.problem.phiFleft(setdiff(obj.sampleInd,[1,2,3])-3,:))*JLnew;
                    Jright=obj.problem.phiFright*pinv(obj.problem.phiFright(setdiff(obj.sampleInd,[1,2,3])-3,:))*JRnew;
                    
                    Dg=[(obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi1(2:end-1,:);
                        (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi2(2:end-1,:);
                        (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi3(2:end-1,:) ]+...
                        obj.time.dt*(Jleft(1:3,1:6)*obj.problem.phi(1:6,:) + ...
                               Jright(end-2:end,end-5:end)*obj.problem.phi(end-5:end,:)- ...
                           [zeros(1,length(w_guess));derivQ;zeros(1,length(w_guess))]);
                    
                    
%                     keyboard
%                     start = 1;
%                     finish = obj.problem.prob.nVol-1;
%                     left = -obj.problem.prob.S(start+1)*[obj.problem.phiRoeF1(start,:)*pinv(obj.problem.phiRoeF1(Index,:))*roeF(1,:)';
%                                     obj.problem.phiRoeF2(start,:)*pinv(obj.problem.phiRoeF2(Index,:))*roeF(2,:)';
%                                     obj.problem.phiRoeF3(start,:)*pinv(obj.problem.phiRoeF3(Index,:))*roeF(3,:)'];
%                     right= obj.problem.prob.S(finish+1)*[obj.problem.phiRoeF1(finish,:)*pinv(obj.problem.phiRoeF1(Index,:))*roeF(1,:)';
%                                     obj.problem.phiRoeF2(finish,:)*pinv(obj.problem.phiRoeF2(Index,:))*roeF(2,:)';
%                                     obj.problem.phiRoeF3(finish,:)*pinv(obj.problem.phiRoeF3(Index,:))*roeF(3,:)'];
%                     JL= -obj.problem.prob.S(start+1)*[obj.problem.phiRoeF1(start,:)*pinv(obj.problem.phiRoeF1(Index,:))*droeF(1,:, :)';
%                                     obj.problem.phiRoeF2(start,:)*pinv(obj.problem.phiRoeF2(Index,:))*droeF(2,:,start)';
%                                     obj.problem.phiRoeF3(start,:)*pinv(obj.problem.phiRoeF3(Index,:))*droeF(3,:,start)'];
%                     JR= obj.problem.prob.S(finish+1)*[obj.problem.phiRoeF1(finish,:)*pinv(obj.problem.phiRoeF1(Index,:))*droeF(1,:, finish)';
%                                     obj.problem.phiRoeF2(finish,:)*pinv(obj.problem.phiRoeF2(Index,:))*droeF(2,:,finish)';
%                                     obj.problem.phiRoeF3(finish,:)*pinv(obj.problem.phiRoeF3(Index,:))*droeF(3,:,finish)'];            
%                     
%                     
% %                    keyboard
%                     g(1:3,:)=[(obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi1(2:end-1,:)*w_guess;
%                         (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi2(2:end-1,:)*w_guess;
%                         (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi3(2:end-1,:)*w_guess ]+...
%                         obj.time.dt*(left + right-...
%                         [(obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*PhiQ1(2:end-1,:)*pinv(PhiQ1(obj.sampleNodes,:))*Q(1,:)';
%                         (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*PhiQ2(2:end-1,:)*pinv(PhiQ2(obj.sampleNodes,:))*Q(2,:)';
%                         (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*PhiQ3(2:end-1,:)*pinv(PhiQ3(obj.sampleNodes,:))*Q(3,:)']);
%                     
                    
                    %                     index=[];
                    %                     for jj=1:length(obj.probGNAT.iarrayFaceI)
                    %                         index=[index, 3*(obj.probGNAT.iarrayFaceI(jj)-1)+1:3*obj.probGNAT.iarrayFaceI(jj)];
                    %                     end
                    %
                    %                     keyboard
                    %                     JRfull=obj.problem.phiFright*pinv(obj.problem.phiFright(index,:))*JnewR;
                    %                     JLfull=obj.problem.phiFleft*pinv(obj.problem.phiFleft(index,:))*JnewL;
%                      keyboard
%                     Dg=[(obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi1(2:end-1,:);
%                         (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi2(2:end-1,:);
%                         (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi3(2:end-1,:) ]+...
%                         obj.time.dt*(JLfull(1:3,1:3)*obj.problem.phi(1:3,:) +JLfull(1:3,4:6)*obj.problem.phi(4:6,:)+ ...
%                         JRfull(end-2:end,1:3)*obj.problem.phi(end-5:end-3,:)+ JRfull(end-2:end,4:6)*obj.problem.phi(end-2:end,:)- [zeros(1,5);derivQ;zeros(1,5)]);
%                  keyboard
%                     Dg=[(obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi1(2:end-1,:);
%                         (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi2(2:end-1,:);
%                         (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi3(2:end-1,:) ]+...
%                         obj.time.dt*(obj.JhatFluxL(4:6,:)*obj.problem.phi(unique(obj.jrow),:) + ...
%                                obj.JhatFluxR(end-2:end,:)*obj.problem.phi(unique(obj.jrow),:)- [zeros(1,length(w_guess));derivQ;zeros(1,length(w_guess))]);
                
%                     JL= obj.problem.phi(4:6,:)* pinv(obj.phiYhat)*obj.JhatFluxL
%                     Dg=[(obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi1(2:end-1,:);
%                         (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi2(2:end-1,:);
%                         (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi3(2:end-1,:) ];
% %                         obj.time.dt*(obj.problem.prob.S(2)*

                
                end
                
                
%                 if obj.probGNAT.indN
%                     
%                     %keyboard
%                     g(1:3,:)=[(obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi1(2:end-1,:)*w_guess;
%                         (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi2(2:end-1,:)*w_guess;
%                         (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi3(2:end-1,:)*w_guess ]+...
%                         obj.time.dt*([PhiFl1(1,:)*pinv(PhiFl1(setdiff(obj.sampleNodes, obj.probGNAT.indN),:))*Rleft(1,:)';
%                         PhiFl2(1,:)*pinv(PhiFl2(setdiff(obj.sampleNodes, obj.probGNAT.indN),:))*Rleft(2,:)';
%                         PhiFl3(1,:)*pinv(PhiFl3(setdiff(obj.sampleNodes, obj.probGNAT.indN),:))*Rleft(3,:)']+ Rright(:,1)-...
%                          (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*PhiQ2(2:end-1,:)*pinv(PhiQ2(obj.problem.prob.sampleNodes,:))*Q(2,:)';
%                         (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*PhiQ3(2:end-1,:)*pinv(PhiQ3(obj.problem.prob.sampleNodes,:))*Q(3,:)']);
%                     
%                 end
%                 keyboard
                H=Df'*Df; h=Df'*Rhat;
                P=H\Dg';
                x=H\h;
                S=-inv(Dg*P);
                Qzim=P*S;
                del_w=Qzim*g-(x+Qzim*Dg*x);
                w_guess=w_guess+del_w;
%keyboard
                if norm(del_w,2)<10^(-6)
                    break;
                end
                obj.sv(:,itnump1)=w_guess;
%                 obj.sv(:,itnump1) = obj.sv(:,itnump1) + w_guess;
                obj.partialU = obj.partialUprev + obj.phiYhat*obj.sv(:,itnump1);
                %Compute residual and jacobian with updated vector for next iteration
                [Rhat,JVhat] = obj.TimeScheme.TimeIntNLFuncGNAT(obj.probGNAT,obj.partialU,obj.partialUprev,t);
                Df=JVhat;
            end
             disp(['norm of the constraint   ', num2str(norm(g))])
        end
        
        


function  [] = NewtonRaphsonconstr_new_method(obj)
            
            itnump1 = obj.cTimeIter + 1;
            t = obj.time.T(1) + obj.time.dt*obj.cTimeIter;
            obj.partialUprev = obj.partialU;
            %Determine the residual and jacobian based on initial guess
            [Rhat,JVhat] = obj.TimeScheme.TimeIntNLFuncGNAT(obj.probGNAT,obj.partialU,obj.partialUprev,t);
%             keyboard
            Df=JVhat;
            w_guess=obj.phiYhat\(obj.partialU-obj.partialUprev);
            
%               keyboard
            for i_N = 1:obj.newt.maxIter %loop until the maximum newton iterations are reached
   
                U = reshape(obj.partialU,3,obj.probGNAT.nVolSamp);
                startInd = 1 + obj.probGNAT.ind1;
                endInd   = obj.probGNAT.nVolSamp - obj.probGNAT.indN;
                [rho,u,P,c] = obj.probGNAT.conservativeToPrimitive(U(:,startInd:endInd));
                e=U(3,startInd:endInd);
                if obj.probGNAT.ind1
                    rho = [U(1,1),rho]; %Density
                    u   = [U(2,1),u]; %Velocity
                    P   = [U(3,1),P]; %Pressure
                    c   = [sqrt(obj.probGNAT.gamma*P(1)/rho(1)),c]; %Speed of sound
                    e   = [P(1)/(obj.probGNAT.gamma-1)+rho(1)*u(1)^2/2,e]; %Energy
                end
                if obj.probGNAT.indN
                    rho = [rho,U(1,end)]; %Density
                    u   = [u,U(2,end)]; %Velocity
                    P   = [P,U(3,end)]; %Pressure
                    c   = [c,sqrt(obj.probGNAT.gamma*P(end)/rho(end))]; %Speed of sound
                    e   = [e,P(end)/(obj.probGNAT.gamma-1)+rho(end)*u(end)^2/2]; %Energy
                end
                
                dc_cons = [0.5*obj.probGNAT.gamma./(c.*rho).*(0.5*(obj.probGNAT.gamma-1)*u.*u - P./rho);...
                    -0.5*obj.probGNAT.gamma*(obj.probGNAT.gamma-1)*u./(rho.*c);...
                    0.5*obj.probGNAT.gamma*(obj.probGNAT.gamma-1)./(rho.*c)]';
                
                %[rho, u, P,c,e,~,dc]=obj.prob.getVariables(obj.sv(:,itnump1));
                [Q,dQ]=obj.probGNAT.forceTermGNAT(u,P);
                %keyboard
                [roeF, droeF]=obj.probGNAT.roeFluxGNAT(rho,u,P,c,e,dc_cons);
                if obj.probGNAT.ind1
                    dUdV=[1,0,0;u(1),rho(1),0;0.5*u(1)*u(1),rho(1)*u(1),1/(obj.probGNAT.gamma-1)];
                    droeF(:,1:3,1)=droeF(:,1:3,1)*dUdV;
                end
                if obj.probGNAT.indN
                    dUdV=[1,0,0;u(end),rho(end),0;0.5*u(end)*u(end),rho(end)*u(end),1/(obj.probGNAT.gamma-1)];
                    droeF(:,4:6,end)=droeF(:,4:6,end)*dUdV;
                end
                
%                 Rright= bsxfun(@times,roeF(:,obj.probGNAT.iarrayFaceI+1),obj.probGNAT.S(obj.probGNAT.iarrayFace(1+obj.probGNAT.ind1:end-obj.probGNAT.indN)+1));
%                 Rleft =-bsxfun(@times,roeF(:,obj.probGNAT.iarrayFaceI),obj.probGNAT.S(obj.probGNAT.iarrayFace(1+obj.probGNAT.ind1:end-obj.probGNAT.indN)));
 
                z=obj.sampleNodes(2:end)-1;
                Dr1L = zeros(obj.probGNAT.nVolMask-1, 900);
                Dr1R = zeros(obj.probGNAT.nVolMask-1, 900);
                Dr2L = zeros(obj.probGNAT.nVolMask-1, 900);
                Dr2R = zeros(obj.probGNAT.nVolMask-1, 900);
                Dr3L = zeros(obj.probGNAT.nVolMask-1, 900);
                Dr3R = zeros(obj.probGNAT.nVolMask-1, 900);
                for k = 1:obj.probGNAT.nVolMask-obj.probGNAT.ind1-obj.probGNAT.indN
                   ii=z(k);
                        
                    Dr1L(k,3*ii-2:3*(ii+1)) =- squeeze(droeF(1,:,obj.probGNAT.iarrayFaceI(k)))';
                    Dr1R(k,3*(ii+1)-2:3*(ii+2)) = squeeze(droeF(1,:,obj.probGNAT.iarrayFaceI(k)+1))';
                    Dr2L(k,3*ii-2:3*(ii+1)) =- squeeze(droeF(2,:,obj.probGNAT.iarrayFaceI(k)))';
                    Dr2R(k,3*(ii+1)-2:3*(ii+2)) = squeeze(droeF(2,:,obj.probGNAT.iarrayFaceI(k)+1))';
                    Dr3L(k,3*ii-2:3*(ii+1)) =- squeeze(droeF(3,:,obj.probGNAT.iarrayFaceI(k)))';
                    Dr3R(k,3*(ii+1)-2:3*(ii+2)) = squeeze(droeF(3,:,obj.probGNAT.iarrayFaceI(k)+1))';
                         
                end
                
%                 keyboard 
                Phi1=obj.problem.phi(1:3:end,:); %basis for the first conserved quantity rho
                Phi2=obj.problem.phi(2:3:end,:); %basis for the second conserved quantity rho*u
                Phi3=obj.problem.phi(3:3:end,:); %basis for the thirs conserved quantity e
                PhiQ1=obj.problem.phiQ(1:3:end,:);
                PhiQ2=obj.problem.phiQ(2:3:end,:);
                PhiQ3=obj.problem.phiQ(3:3:end,:);
                %dQ has only the derivative of Q2, since others are 0
%               keyboard
                derQ=PhiQ2*pinv(PhiQ2(obj.sampleNodes,:))*squeeze(dQ(1,:,:))'; %nVol x 3 matrix, 3 corresponds to derivative wrt rho, rho*u and e                
                derivQ=zeros(1, length(w_guess));
                
                for i=2:obj.problem.prob.nVol-1
                    derivQ=derivQ+derQ(i,:)*[Phi1(i,:);Phi2(i,:);Phi3(i,:)]*(obj.problem.prob.SVol(i).*obj.problem.prob.dx(i));
                end
                
%                 keyboard
%                 
%                 load romroeF
%                 reshF=reshape(romroeF(:,1),3, 299);
%                 Index=[];
%                 for jj=1:299
%                     for kk=1:size(roeF,2)
%                         if abs(roeF(:,kk)-reshF(:,jj))<10^(-10)
%                             Index=[Index,jj];
%                         end
%                     end
%                 end
                
                if obj.probGNAT.ind1
                                         
%                     Dg=[(obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi1(2:end-1,:);
%                         (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi2(2:end-1,:);
%                         (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi3(2:end-1,:) ]+...
%                         obj.time.dt*(Jleft(1:3,1:6)*obj.problem.phi(1:6,:) + ...
%                                Jright(end-2:end,end-5:end)*obj.problem.phi(end-5:end,:)- ...
%                            [zeros(1,length(w_guess));derivQ;zeros(1,length(w_guess))]);
                    
                    
%                     keyboard
                    start = 1;
                    finish = obj.problem.prob.nVol-1;
                    left = -obj.problem.prob.S(start+1)*[obj.problem.phiRoeF1(start,:)*pinv(obj.problem.phiRoeF1(obj.sampleNodes(2:end)-1,:))*roeF(1,obj.probGNAT.iarrayFaceI)';
                                    obj.problem.phiRoeF2(start,:)*pinv(obj.problem.phiRoeF2(obj.sampleNodes(2:end)-1,:))*roeF(2,obj.probGNAT.iarrayFaceI)';
                                    obj.problem.phiRoeF3(start,:)*pinv(obj.problem.phiRoeF3(obj.sampleNodes(2:end)-1,:))*roeF(3,obj.probGNAT.iarrayFaceI)'];
                    right= obj.problem.prob.S(finish+1)*[obj.problem.phiRoeF1(finish,:)*pinv(obj.problem.phiRoeF1(obj.sampleNodes(2:end),:))*roeF(1,obj.probGNAT.iarrayFaceI+1)';
                                    obj.problem.phiRoeF2(finish,:)*pinv(obj.problem.phiRoeF2(obj.sampleNodes(2:end),:))*roeF(2,obj.probGNAT.iarrayFaceI+1)';
                                    obj.problem.phiRoeF3(finish,:)*pinv(obj.problem.phiRoeF3(obj.sampleNodes(2:end),:))*roeF(3,obj.probGNAT.iarrayFaceI+1)'];
                    
                    JL= -obj.problem.prob.S(start+1)*[obj.problem.phiRoeF1(start,:)*pinv(obj.problem.phiRoeF1(obj.sampleNodes(2:end)-1,:))*Dr1L;
                                    obj.problem.phiRoeF2(start,:)*pinv(obj.problem.phiRoeF2(obj.sampleNodes(2:end)-1,:))*Dr2L;
                                    obj.problem.phiRoeF3(start,:)*pinv(obj.problem.phiRoeF3(obj.sampleNodes(2:end)-1,:))*Dr3L];
                    JR= obj.problem.prob.S(finish+1)*[obj.problem.phiRoeF1(finish,:)*pinv(obj.problem.phiRoeF1(obj.sampleNodes(2:end),:))*Dr1R;
                                    obj.problem.phiRoeF2(finish,:)*pinv(obj.problem.phiRoeF2(obj.sampleNodes(2:end),:))*Dr2R;
                                    obj.problem.phiRoeF3(finish,:)*pinv(obj.problem.phiRoeF3(obj.sampleNodes(2:end),:))*Dr3R];            
                    
                    
%                    keyboard
                    g(1:3,:)=[(obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi1(2:end-1,:)*w_guess;
                        (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi2(2:end-1,:)*w_guess;
                        (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi3(2:end-1,:)*w_guess ]+...
                        obj.time.dt*(left + right-...
                        [(obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*PhiQ1(2:end-1,:)*pinv(PhiQ1(obj.sampleNodes,:))*Q(1,:)';
                        (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*PhiQ2(2:end-1,:)*pinv(PhiQ2(obj.sampleNodes,:))*Q(2,:)';
                        (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*PhiQ3(2:end-1,:)*pinv(PhiQ3(obj.sampleNodes,:))*Q(3,:)']);
                    
                    Dg=[(obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi1(2:end-1,:);
                        (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi2(2:end-1,:);
                        (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi3(2:end-1,:) ]+...
                        obj.time.dt*(JL*obj.problem.phi + JR*obj.problem.phi- [zeros(1,length(w_guess));derivQ;zeros(1,length(w_guess))]);
                    
                    %                     index=[];
                    %                     for jj=1:length(obj.probGNAT.iarrayFaceI)
                    %                         index=[index, 3*(obj.probGNAT.iarrayFaceI(jj)-1)+1:3*obj.probGNAT.iarrayFaceI(jj)];
                    %                     end
                    %
                    %                     keyboard
                    %                     JRfull=obj.problem.phiFright*pinv(obj.problem.phiFright(index,:))*JnewR;
                    %                     JLfull=obj.problem.phiFleft*pinv(obj.problem.phiFleft(index,:))*JnewL;
%                      keyboard
%                     Dg=[(obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi1(2:end-1,:);
%                         (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi2(2:end-1,:);
%                         (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi3(2:end-1,:) ]+...
%                         obj.time.dt*(JLfull(1:3,1:3)*obj.problem.phi(1:3,:) +JLfull(1:3,4:6)*obj.problem.phi(4:6,:)+ ...
%                         JRfull(end-2:end,1:3)*obj.problem.phi(end-5:end-3,:)+ JRfull(end-2:end,4:6)*obj.problem.phi(end-2:end,:)- [zeros(1,5);derivQ;zeros(1,5)]);
%                  keyboard
%                     Dg=[(obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi1(2:end-1,:);
%                         (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi2(2:end-1,:);
%                         (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi3(2:end-1,:) ]+...
%                         obj.time.dt*(obj.JhatFluxL(4:6,:)*obj.problem.phi(unique(obj.jrow),:) + ...
%                                obj.JhatFluxR(end-2:end,:)*obj.problem.phi(unique(obj.jrow),:)- [zeros(1,length(w_guess));derivQ;zeros(1,length(w_guess))]);
                
%                     JL= obj.problem.phi(4:6,:)* pinv(obj.phiYhat)*obj.JhatFluxL
%                     Dg=[(obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi1(2:end-1,:);
%                         (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi2(2:end-1,:);
%                         (obj.problem.prob.SVol(2:end-1).*obj.problem.prob.dx(2:end-1))*Phi3(2:end-1,:) ];
% %                         obj.time.dt*(obj.problem.prob.S(2)*

                
                end
                
%                 keyboard
                H=Df'*Df; h=Df'*Rhat;
                P=H\Dg';
                x=H\h;
                S=-inv(Dg*P);
                Qzim=P*S;
                del_w=Qzim*g-(x+Qzim*Dg*x);
                w_guess=w_guess+del_w;
%keyboard
                if norm(del_w,2)<10^(-6)
                    break;
                end
                obj.sv(:,itnump1)=w_guess;
%                 obj.sv(:,itnump1) = obj.sv(:,itnump1) + w_guess;
                obj.partialU = obj.partialUprev + obj.phiYhat*obj.sv(:,itnump1);
                %Compute residual and jacobian with updated vector for next iteration
                [Rhat,JVhat] = obj.TimeScheme.TimeIntNLFuncGNAT(obj.probGNAT,obj.partialU,obj.partialUprev,t);
                Df=JVhat;
            end
             disp(['norm of the constraint   ', num2str(norm(g))])
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        function  [RAbar,E] = updateLocOperators(obj,Rhat,JVhat)
            %This function updates the local operators (reduced size)
            %during the online phase.  From Algorithm 3 of Carlberg et. al.
            %2010.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj   - GNAT object
            %Rhat  - partial residual
            %JVhat - partial jacobian multiplied by the partial phi
            %
            %Outputs:
            %--------
            %RAbar - RAbar from Algorithm 3 of Carlberg et. al. 2010
            %E     - Q'*Bbar from Algorithm 3 of Carlberg et. al 2010
            %--------------------------------------------------------------
            
            [Q,RAbar] = qr(obj.A*JVhat,0);
            E = Q'*(obj.B*Rhat);
            
        end
        
        function  [JVhat] = computeJPHIhat(obj,Jhat)
            %This function multiplies the partial Jacobian by the partial
            %phi matrix.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj   - GNAT object
            %Jhat  - rows of Jacobian matrix corresponding to sample
            %        indices stored in a vector
            %
            %Outputs:
            %--------
            %JVhat - product of partial Jacobian with partial phi
            %--------------------------------------------------------------
            
            %Reshape Jhat into matrix JHatTemp
            %keyboard
            obj.JhatTemp(obj.reconstJhatInd) = Jhat;
            JVhat = obj.JhatTemp*obj.phiYhat;
            
        end
        
        function  [] = modifyICandTranslateBasisRk1(obj,newIC_hat)
            
            if obj.nBases > 1
                error('IC modification and basis translation only works for global Basis FOR NOW');
            end
            
            if norm(newIC_hat-obj.svdUpdateData.wRef) < 1e-4
                fprintf('Skipping basis update...\n');
                return;
            end
            
            obj.probGNAT.setProperty('ic',newIC_hat);
            obj.phiYhat=ThinSVD_AddVec2Cols_Fast(obj.phiYhat,obj.SingVals(1:obj.nY),...
                obj.svdUpdateData.ColSumRSVtrunc,...
                obj.svdUpdateData.normOneMinusVVt,...
                obj.svdUpdateData.wRef-newIC_hat);
            obj.phiYhat=obj.phiYhat(:,1:end-1);
            
            obj.svdUpdateData.wRef=newIC_hat;
            if length(obj.SingVals)==obj.nY, err=0; else err = 100*obj.SingVals(obj.nY+1)/obj.SingVals(1); end;
            fprintf('Relative 2-norm error in update = %f%%\n',err);
            
        end
        
        %Optimization
        function   [f,df] = ObjFMINCON_SAND(obj,u)
            
            w=u(1:obj.nY,1);
            p=u(obj.nY+1:end,1);
            
            obj.probGNAT.updateParametersGNAT(p);
            
            [f,dfdw,dfdp] = obj.objectiveSAND(w, p);
            
            
            
            df = [obj.phiYhat'*dfdw;dfdp];
            %disp('Current Value of p = ');
            %disp(p);
        end
        
        function   [f,df] = ObjFMINCON_SAND_ROMVAR(obj,u)
            
            w=u(1:obj.nY,1);
            p=u(obj.nY+1:end,1);
            
            obj.probGNAT.updateParametersGNAT(p);
            
            [f,dfdw,dfdp] = obj.objectiveSAND_ROMvar(w, p);
            
            
            % disp('To change!');
            %df = [obj.phiYhat'*dfdw;dfdp];
            df = [dfdw;dfdp];
            %disp('Current Value of p = ');
            %disp(p);
        end
        
        function   [f,df] = ObjFMINCON_NAND(obj,p,sensMethod)
            if norm(p - obj.curr_param) > 0
                obj.probGNAT.updateParametersGNAT(p);
                obj.executeModel;
            end
            
            [df,f] = obj.reducedGrad(p,sensMethod);
            if (~isreal(f))
                f = NaN;
            end
            if (~isreal(df))
                df = zeros(size(df));
            end
        end
        
        function   [cineq,ceq,dcineq,dceq] = NLConstFMINCON_SAND(obj,u)
            
            w=u(1:obj.nY,1);
            p=u(obj.nY+1:end,1);
            obj.probGNAT.updateParametersGNAT(p);
            
            if strcmpi(obj.Rtype,'Galerkin')
                [RHat,dRdwHat,dRdpHat] = obj.probGNAT.ResSensGNAT(obj.probGNAT.ic(obj.uniqueJROWind) + obj.phiYhat(obj.uniqueJROWind,:)*w,[]);
                obj.JhatTemp(obj.reconstJhatInd) = dRdwHat;
                
                ceq = obj.B*RHat;
                dceq = [obj.B*obj.JhatTemp*obj.phiYhat, obj.B*dRdpHat]';
                
                cineq=[];
                dcineq=[];
            else
                [RHat,dRdwHat,dRdpHat] = obj.probGNAT.ResSensGNAT(obj.probGNAT.ic(obj.uniqueJROWind) + obj.phiYhat(obj.uniqueJROWind,:)*w,[]);
                obj.JhatTemp(obj.reconstJhatInd) = dRdwHat;
                
                ceq = obj.phiYhat'*obj.JhatTemp'*obj.Ebar*RHat;
                dRdwPhi = obj.JhatTemp*obj.phiYhat;
                dceq = [(dRdwPhi'*obj.Ebar*dRdwPhi), dRdwPhi'*obj.Ebar*dRdpHat]';
                
                cineq=[];
                dcineq=[];
            end
        end
        
        function  [DfDp,f,lambda] = reducedGrad(obj,z,sensMethod)
            %This function computes the reduced gradient of the objective
            %function using either the direct or adjoint method to compute
            %the sensitivity of the state with respect to the parameter.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj    - ROM object
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
            %--------------------------------------------------------------
            
            if obj.killflag
                f  = nan;
                DfDp = nan(size(z));
                lambda=[];
                return;
            end
            
            what = obj.probGNAT.ic + obj.phiYhat*sum(obj.sv,2);
            
            %Compute the residual sensitivities
            [~,dRhdw,dRhdp,~,~] = obj.probGNAT.ResSensGNAT(what(obj.uniqueJROWind));
            %Evaluate objective and sensitivities
            [f,dfdw,dfdp] = obj.objectiveNAND(what,z);
            
            obj.JhatTemp(obj.reconstJhatInd) = dRhdw;
            
            switch lower(obj.Rtype)
                case 'petrov-galerkin'
                    
                    dRrdy = obj.phiYhat'*(obj.JhatTemp'*(obj.Ebar*(obj.JhatTemp*obj.phiYhat)));
                    dRrdp = obj.phiYhat'*(obj.JhatTemp'*(obj.Ebar*dRhdp));
                    
                case 'galerkin'
                    
                    dRrdy = obj.B*(obj.JhatTemp*obj.phiYhat);
                    dRrdp = obj.B*dRhdp;
            end
            
            switch sensMethod
                case 'direct'
                    %Using the direct method, solve for the state
                    %sensitivity (w.r.t. parameter) and compute the reduced
                    %gradient
                    DwDp = -dRrdy\dRrdp;
                    DfDp = dfdp + DwDp'*(obj.phi'*dfdw);
                    %DwDp = -dRdw\dRdp;
                    %DfDp = dfdp + DwDp'*dfdw;
                    lambda = [];
                case 'adjoint'
                    %Using the adjoint method, solve for the dual variable
                    %and compute the reduced gradient
                    lambda = -(dRrdy'\(obj.phiYhat'*dfdw));
                    DfDp = dfdp + dRrdp'*lambda;
                    %lambda = -dRdw'\dfdw;
                    %DfDp = dfdp + dRdp'*lambda;
            end
        end
        
        function  [] = setOptObjective(obj,txt)
            %This function computes the reduced gradient of the objective
            %function using either the direct or adjoint method to compute
            %the sensitivity of the state with respect to the parameter.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj - ROM object
            %txt - string containing the name of the objective function
            %      which is stored in obj.prob
            %
            %Outputs:
            %--------
            %There are no outputs.
            %--------------------------------------------------------------
            %disp('To change!');
            obj.objectiveSAND_ROMvar = eval(['@(u,z) obj.probGNAT.',txt,'(u,z)']);
            obj.objectiveSAND = eval(['@(u,z) obj.probGNAT.',txt,'(obj.probGNAT.ic + obj.phiYhat*u,z)']);
            obj.objectiveNAND = eval(['@(u,z) obj.probGNAT.',txt,'(u,z)']);
        end
        
        %Unused?
        function  [nRnew,nJnew,nInew,bestSpeed,err] = improveGNATparams(obj,fomobj,romobj,nRdata,nJdata,nIdata,maxErr,normType,type,stopFlag)
            %This function determines a better selection (in terms of
            %speed) of nR, nJ, nI from the current selection (inside some
            %3D rectangle specified by the bounds) such that the maximum
            %relative error is below some acceptable level specified by
            %maxErr.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj      - GNAT object
            %fomobj   - FOM object
            %romobj   - ROM object
            %nRdata   - 1 x 3 vector indicating the lower and upper bounds
            %           of the nR search space, and the frequency at which
            %           to check points, i.e. [20, 50, 10] means that nR
            %           will be tested at 20, 30, 40, 50
            %nJdata   - 1 x 3 vector indicating the lower and upper bounds
            %           of the nJ search space, and the frequency at which
            %           to check points, i.e. [20, 50, 10] means that nJ
            %           will be tested at 20, 30, 40, 50
            %nIdata   - 1 x 3 vector indicating the lower and upper bounds
            %           of the nI search space, and the frequency at which
            %           to check points, i.e. [20, 50, 10] means that nI
            %           will be tested at 20, 30, 40, 50
            %maxErr   - scalar indicating maximum acceptable max. relative
            %           error.
            %normType - scalar or 'inf' indicating the norm type to use in
            %           quantifying the error.
            %type     - string indicating type of improvement algorithm to
            %           perform (only 'global' supported currently).
            %stopFlag - boolean indicating whether or not to stop searching
            %           once a better nR, nJ, nI are found.
            %
            %Outputs:
            %--------
            %nRnew     - new value of nR
            %nJnew     - new value of nJ
            %nInew     - new value of nI
            %bestSpeed - speedup from FOM using the new parameters
            %err       - maximum relative error using new parameters
            %--------------------------------------------------------------
            
            %Compute phiR and phiJ
            [phiR,~] = pod(romobj.res);
            [phiJ,~] = pod(romobj.jac);
            
            %Run the simulation with the original nR, nJ, nI combination
            fprintf('Re-running original parameter combination so speedup parameters can be compared fairly\n');
            
            %Make sure user doesn't specify invalid property combinations.
            if obj.nR < obj.nY
                fprintf('The input nR was not valid, changing from %u to %u\n',obj.nR,obj.nY);
                obj.nR = obj.nY;
            end
            if obj.nJ < obj.nY
                fprintf('The input nJ was not valid, changing from %u to %u\n',obj.nJ,obj.nY);
                obj.nJ = obj.nY;
            end
            if obj.nR > obj.nI || obj.nJ > obj.nI
                fprintf('The input nI was not valid, changing from %u to %u\n',obj.nI,max(obj.nR,obj.nJ));
                obj.nI = max(obj.nR,obj.nJ);
            end
            
            newGnat = GNAT([],[],[],obj);
            newGnat.createGNAT(fomobj.prob,phiR,phiJ,romobj.phi);
            newGnat.executeModel;
            
            %Set vectors that contain all of the points to be checked in
            %each dimension.
            z1 = nRdata(1):nRdata(3):nRdata(2);
            z2 = nJdata(1):nJdata(3):nJdata(2);
            z3 = nIdata(1):nIdata(3):nIdata(2);
            
            %Determine the speedup of the initial nR, nJ, nI combination
            %and set a bestSpeed variable.
            speed0 = fomobj.ontime/newGnat.ontime;
            bestSpeed = speed0;
            SV0  = newGnat.reconstructFullState(romobj);
            err0 = max(computeError(SV0,fomobj.sv,normType));
            
            %Record initial nR, nJ, nI values
            nR0 = obj.nR;
            nJ0 = obj.nJ;
            nI0 = obj.nI;
            
            %Determine ROM error to check if the specified maxErr has a
            %chance to be realized.
            errROM = max(computeError(romobj.sv,fomobj.sv,normType));
            if errROM > maxErr
                disp(['Feasible set empty since the ROM error is greater than the maximum error specified.  Increase maxErr to ',num2str(errROM),' + e, for some e > 0']);
                %This will ensure that the maximum error between the ROM
                %and GNAT is less than 2%, but I have quantified the error
                %with reference to the FOM, which is why the conversion
                %below is used.
                temp = max(ColumnwiseNorm(romobj.sv,normType)./ColumnwiseNorm(fomobj.sv,normType));
                maxErr = 0.02*temp + errROM;
                disp(['The maximum error has been changed to ',num2str(maxErr),'(the error between the GNAT and the FOM) which will give you a GNAT model within 2% of the ROM object']);
            end;
            
            %Count the total number of configurations of nR, nJ, nI that
            %will be check.
            cnt = 0;
            for i = z1
                for j = z2
                    for k = z3
                        if (i > k) || (j > k) || (i < newGnat.nY) || (j < newGnat.nY)
                            continue;
                        end
                        cnt = cnt + 1;
                    end
                end
            end
            
            %Initialize the outputs in case a better combination is not
            %found
            nRnew = obj.nR;
            nJnew = obj.nJ;
            nInew = obj.nI;
            err = maxErr;
            
            %Initialize iteration count
            iter = 0;
            switch type
                case 'global'
                    for j = z2 %There is relative insensitivity to nJ and will thus be the outer loop
                        for k = z3 %nI is important for speed and will be the intermediate loop
                            for i = z1 %nR is usually the most important for error control,
                                %so this will be the inside loop
                                %If nR > nI or nJ > nI or nR < nY or nJ <
                                %nY skip the iteration
                                if i > k || j > k || i < newGnat.nY || j < newGnat.nY
                                    continue;
                                end
                                iter = iter + 1;
                                fprintf('GNAT Test %u of %u\n',iter,cnt);
                                
                                %Clear newGnat object and reinitialize.
                                %Set new nR, nJ, nI
                                clear newGnat;
                                newGnat = GNAT([],[],[],obj);
                                newGnat.nR = i;
                                newGnat.nJ = j;
                                newGnat.nI = k;
                                newGnat.nSample = [];
                                newGnat.nGreed  = [];
                                
                                %Perform offline and online computations
                                newGnat.createGNAT(fomobj.prob,phiR,phiJ,romobj.phi);
                                newGnat.executeModel;
                                
                                %Extract the state vector, maximum relative
                                %error, and speedup
                                newSV  = newGnat.reconstructFullState(romobj);
                                newErr = max(computeError(newSV,fomobj.sv,normType));
                                newSpeed = fomobj.ontime/newGnat.ontime;
                                
                                if newErr <= maxErr
                                    %If we have found a feasible point
                                    %(i.e. MRerr < maxErr)...
                                    if (err0 > maxErr) && (newSpeed <= bestSpeed) && (nRnew == nR0 && nJnew == nJ0 && nInew == nI0)
                                        %If the original error is above the
                                        %tolerance, the speedup is not
                                        %better than the original (b/c the
                                        %original used parameters that were
                                        %too small) and this is the first
                                        %feasible point encounter (this is
                                        %the last condition), perform
                                        %the code below.
                                        
                                        %Store new 'best' values
                                        nRnew = i; nJnew = j; nInew = k;
                                        bestSpeed = newSpeed;
                                        err = newErr;
                                        
                                        %Update the input GNAT object to have
                                        %these 'best' values
                                        obj.nR = nRnew;
                                        obj.nJ = nJnew;
                                        obj.nI = nInew;
                                        
                                        fprintf('new Speedup = %.5f vs. original Speedup = %.5f\n',bestSpeed,speed0);
                                        continue;
                                        
                                    elseif newSpeed > bestSpeed
                                        %Also perform the code below if the
                                        %speedup of this parameter
                                        %combination is faster, which means
                                        %we have found a better combo.
                                    else
                                        %Otherwise, move to the next loop
                                        %iteration.
                                        fprintf('Did not find a better GNAT parameter selection\n');
                                        continue;
                                    end
                                    
                                    %Store new 'best' values
                                    nRnew = i; nJnew = j; nInew = k;
                                    bestSpeed = newSpeed;
                                    err = newErr;
                                    
                                    %Update the input GNAT object to have
                                    %these 'best' values
                                    obj.nR = nRnew;
                                    obj.nJ = nJnew;
                                    obj.nI = nInew;
                                    
                                    %Print progress to the user
                                    fprintf('Found a Better GNAT parameter selection!\n');
                                    fprintf('(nR, nJ, nI) = (%u,%u,%u)\n',i,j,k);
                                    fprintf('new Speedup = %.5f vs. original Speedup = %.5f\n',bestSpeed,speed0);
                                    fprintf('new Error = %.5f vs. maximum acceptable error = %.5E \n',err,maxErr);
                                    fprintf('NOTE: When this function call is completed, the GNAT object you passed in will have nR, nJ, nI set to the values indicated above\n');
                                    if stopFlag
                                        %If the user wants the function to
                                        %exit when it finds any better
                                        %combination that the original,
                                        %return.
                                        return;
                                    end
                                else
                                    fprintf('Did not find a better GNAT parameter selection\n');
                                end
                            end
                        end
                    end
            end
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
            
            all2clear4PP = {'A','B','irstart','jrow','jdiag','phiYhat',...
                'partialU','partialUprev','JhatTemp',...
                'reconstJhatInd','uniqueJROWind','jdiagHat'};
            
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
        
        function  [] = adjustProp_nRnJnI(obj,nRnew,nJnew,nInew)
            obj.nR = nRnew;
            obj.nJ = nJnew;
            obj.nI = nInew;
        end
        
        function  [] = copyProperties(obj,oldobj)
            %This function copies the properties from oldobj to obj, where both oldobj
            %and obj are instances of a Problem class.
            %--------------------------------------------------------------------------
            %Inputs:
            %-------
            %obj    - Problem class
            %oldobj - Problem class whose properties are to be copied into obj
            %
            %Outputs:
            %--------
            %There are no outptus.
            %--------------------------------------------------------------------------
            
            props = properties(oldobj);
            for i = 1:length(props)
                % Use Dynamic Expressions to copy the required property.
                % For more info on usage of Dynamic Expressions, refer to
                % the section "Creating Field Names Dynamically" in:
                % web([docroot '/techdoc/matlab_prog/br04bw6-38.html#br1v5a9-1'])
                obj.(props{i}) = oldobj.(props{i});
            end
            
            
        end
    end
end

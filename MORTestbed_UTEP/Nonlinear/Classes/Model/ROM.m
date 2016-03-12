classdef ROM < handle
    
    properties (SetAccess = private, GetAccess = public)
        %Input properties
        fileRoot = [pwd,SlashBS];
        prob; %Problem object
        Rtype = 'G'; %ROM type: 'G' or 'PG'
        precompFlag = false;
        
        id; %This is the id of the Rtype ROM module
        cfgid; %This is the id of the configuration that
        snaps; %structure with fields: coll, int, num, dist, distparams, distseed
        time;
        newt; %structure with fields: maxIter, eps, nIter (average # of newton steps)
        ndof;
        nY;
        nYrelEnergy=0.9997;
        nYMin;
        nYMax;
        nYflag;
        saveNL = 0;
        nBases = 2;
        addElemTol = 0.1;
        basisSelect = 'closest';
        clustering = 'k_means';
        TimeScheme;
        lastTimeScheme;
        projerr;
        mu;
        sigma;
        FullBasis;
        %trunc=10;
        %ncell=3;
        
        gZim;
        gOrig;
        constr;
        
        curr_param;
        augment=false;
        
        %Online SVD update properties
        basisUpdate;
        snapMATRIX;
        svdUpdateData;
        ROBcompon;
        numswitch;
        prevSwitchState;
        wr_switch;
        
        ClusterIndex;
        SingVals;
        
        objectiveNAND;
        objectiveSAND;
        objectiveSAND_ROMvar;
        fom_samples=[];
        
        %Computed Properties
        ontime;
        sv;
        res = [];
        jac = [];
        phi0; %Cell array of the original bases (before svd updating)
        phi; %The current online reduced basis (double matrix)
        phiFright;
        
        
        
        phiFleft;
        phiRoeF1;
        phiRoeF2;
        phiRoeF3;
        phiQ;
        phiFrnew;
        phiFlnew;
        LocBasisHist;
        clusCenter;
        preCompDistQuant;
        UrLoc;
        cTimeIter; %current time iteration number
        resBound = [];
        errBound = [];
        
        numExecute=0;
        newtonSolver ; % 1 = original,  2= Zimmerman's method;  3=  more constraints
    end
    
    properties(SetAccess=public,GetAccess=public)
        ncell;
        trunc;
    end
    
    properties (Hidden = true, SetAccess = private, GetAccess = public)
        %Global variable text
        GVARtxt = [];
        killflag = false;
        
        numres=0;
        printLevel=1;
    end
    
    methods
        %Constructor
        function [obj] = ROM(ROMfile,romtype,romid,probobj, method, cellnum, basNum)
            %This is the constructor of the ROM class.  It reads the
            %appropriate entries from the input file and stores them in the
            %class instance.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %ROMfile - string indicating the filename of the rom file to
            %          use
            %romtype - string indicating the type of projection to use ('g'
            %          for Galerkin or 'pg' for Petrov-Galerkin)
            %romid   - integer containing the id number of the romtype ROM
            %          object to use
            %probobj - instance of a problem object containing
            %          problem-specific information for the analysis
            %
            %Outputs:
            %--------
            %obj     - instance of the ROM object that was constructed
            %--------------------------------------------------------------
            
            if nargin == 0
                return;
            end
            
            % Construct a new object based on a deep copy of an old object
            % of this class by copying properties over.
            if nargin == 2
                oldobj = romtype;
                if isa(oldobj,'ROM')
                    props = properties(oldobj);
                    for i = 1:length(props)
                        % Use Dynamic Expressions to copy the required property.
                        % For more info on usage of Dynamic Expressions, refer to
                        % the section "Creating Field Names Dynamically" in:
                        % web([docroot '/techdoc/matlab_prog/br04bw6-38.html#br1v5a9-1'])
                        obj.(props{i}) = oldobj.(props{i});
                    end
                    return;
                end
            end
            
            %Store handle to problem in instance
            obj.prob = probobj;
            obj.ncell = cellnum;
            obj.newtonSolver = method;
            obj.trunc = basNum;
           % keyboard
	    %Copy the time structure and store in the class instance
            %The default time structure is whatever is stored in CONFIG
            obj.time = probobj.config.time;
            obj.newt = struct('maxIter',[],'eps',[],'iter',[],'quiet',[]);
            %Initialize snapshot structure
            obj.snaps = struct('coll',[],'int',[],'num',[],'distseed',[],'dist',[],'distparams',[]);
            
            %Extract parameters from ROM file
            obj.GVARtxt = probobj.config.GVARtxt;
            VARtext = readInFile('VAR',ROMfile,1);
            %Extract the ROM text from the rom file (depending on the value
            %of romtype).  Notice, the third argument is a 1 because I only
            %allow 1 block of each type in a rom file.
            if strcmpi(romtype,'g')
                obj.Rtype = 'Galerkin';
                ROMtext    = readInFile('GROM',ROMfile,1);
            elseif strcmpi(romtype,'pg')
                obj.Rtype = 'Petrov-Galerkin';
                ROMtext    = readInFile('PGROM',ROMfile,1);
            elseif strcmpi(romtype,'pg-reg')
                obj.Rtype = 'Petrov-Galerkin-Regularized';
                ROMtext    = readInFile('PGROM',ROMfile,1);
            end
            
            %Determine the ROM properties based on the text in the ROM
            %section of the input file
            determineROM(obj,ROMtext,[obj.GVARtxt;VARtext],romid);
            
            %Extract and store the configuration id to which this ROM
            %instance is associated with.
            obj.cfgid = probobj.config.id;
            
            %             %Set up state vector and store the initial condition.
            %             obj.sv = zeros(probobj.config.ndof,obj.time.nstep+1);
            %             obj.sv(:,1) = probobj.ic;
            
            %Set up number of dofs
            obj.ndof = probobj.config.ndof;
            
            %Set res and jac properties
            if obj.nBases == 1
                obj.res = [];
                obj.jac = [];
            else
                obj.res = cell(obj.nBases,1);
                obj.jac = cell(obj.nBases,1);
            end
        end
        
        function  [] = determineROM(obj,ROMtext,VARtxt,id)
            %This function computes and stores properties of ROM class from
            %the char array ROMtext
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - instance of ROM class
            %ROMtxt  - a string containing the text from the ROM object of
            %          interest in the input file
            %id      - integer indicating the id number of the ROM object
            %          of interest from the input file
            %
            %Outputs:
            %--------
            %No outputs.  Properties are stored in ROM handle class.
            %--------------------------------------------------------------
            
            %Extract variables from VARtxt and evaluate them in the current
            %workspace.
            if ~isempty(VARtxt)
                for i = 1:length(VARtxt)
                    eval([VARtxt{i},';']);
                end
            end
            
            %All of the below properties may be vectors (or cell arrays),
            %in which case, we ensure that they are of the appropriate size
            %and extract the current one:
            
            %%%%%%%%% Determine and store configuration id %%%%%%%%%
            tmp = extractInputRobust(ROMtext,'id',[]);
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
            
            %%%%%%%%% Determine fileRoot %%%%%%%%%
            obj.fileRoot = [extractInputROM(ROMtext,N,obj.id,'fileRoot',obj.fileRoot),'_ID',num2str(id),'_'];
            
            %%%%%%%%% Determine number of Bases %%%%%%%%%
            obj.nBases = extractInputROM(ROMtext,N,obj.id,'nBases',obj.nBases);
            
            obj.numres=zeros(obj.nBases,1);
            obj.printLevel = extractInputROM(ROMtext,N,obj.id,'printLevel',obj.printLevel);
            
            %%%%%%%%% Determine nY, nYrelEnergy, nYMinMax, nYflag %%%%%%%%%
            obj.nY = extractInputROM(ROMtext,N,obj.id,'nY',[]); obj.nY = obj.nY(:);
            if isempty(obj.nY)
                %If obj.nY is empty, we will use energy to determine the size of the ROB
                obj.nYflag = 2;
                
                obj.nYrelEnergy = extractInputROM(ROMtext,N,obj.id,'nYrelEnergy',obj.nYrelEnergy);
                obj.nYMin       = extractInputROM(ROMtext,N,obj.id,'nYMin',[]);
                obj.nYMax       = extractInputROM(ROMtext,N,obj.id,'nYMax',[]);
                %Make size(nYMin,1) == obj.nBases if nYMin is a scalar
                if size(obj.nYMin,1) == 1
                    obj.nYMin = obj.nYMin*ones(obj.nBases,1);
                end
                %Error if nYMin and nBases not compatible
                if size(obj.nYMin,1) ~= obj.nBases
                    error(['Size of nYMin and value of nBases not compatible for id ',num2str(obj.id)]);
                end
                
                %Make size(nYMax,1) == obj.nBases if nYMax is a scalar
                if size(obj.nYMax,1) == 1
                    obj.nYMax = obj.nYMax*ones(obj.nBases,1);
                end
                %Error if nYMax and nBases not compatible
                if size(obj.nYMax,1) ~= obj.nBases
                    error(['Size of nYMax and value of nBases not compatible for id ',num2str(obj.id)]);
                end
            else
                %If obj.nY is not empty, we will use this to determine the size of the ROB (not energy)
                obj.nYflag = 1;
                %Make size(nY,1) == obj.nBases if nY is a scalar
                if size(obj.nY,1) == 1
                    obj.nY = obj.nY*ones(obj.nBases,1);
                end
                %Error if nY and nBases not compatible
                if size(obj.nY,1) ~= obj.nBases
                    error(['Size of nY and value of nBases not compatible for id ',num2str(obj.id)]);
                end
            end
            
            % Regularization
            obj.mu    = extractInputROM(ROMtext,N,obj.id,'mu',obj.mu);
            obj.sigma = extractInputROM(ROMtext,N,obj.id,'sigma',obj.sigma);
            
            %%%%%%%%% Determine type of nonlinear snaps to save %%%%%%%%%
            obj.saveNL = extractInputROM(ROMtext,N,obj.id,'saveNL',obj.saveNL);
            
            %%%%%%%%% Determine type of nonlinear snaps to save %%%%%%%%%
            if (obj.nBases > 1)
                obj.clustering = extractInputROM(ROMtext,N,obj.id,'clustering',obj.clustering);
            else
                obj.clustering = [];
            end
            %%%%%%%%% Determine sharing tolerance %%%%%%%%%
            obj.addElemTol = extractInputROM(ROMtext,N,obj.id,'addElemTol',obj.addElemTol);
            
            %%%%%%%%% Determine method of selecting local basis %%%%%%%%%
            obj.basisSelect = extractInputROM(ROMtext,N,obj.id,'basisSelect',obj.basisSelect);
            
            %%%%%%%%% Determine time structure %%%%%%%%%
            %obj.time = struct('T',[],'dt',[],'nstep',[],'quiet',[]);
            obj.time.T     = extractInputROM(ROMtext,N,obj.id,'T',obj.time.T);
            obj.time.dt    = extractInputROM(ROMtext,N,obj.id,'dt',obj.time.dt);
            obj.time.nstep = extractInputROM(ROMtext,N,obj.id,'nstep',obj.time.nstep);
            obj.time.quiet = extractInputROM(ROMtext,N,obj.id,'timeQuiet',obj.time.quiet);
            obj.time.steadyconverge = extractInputROM(ROMtext,N,obj.id,'steadyconverge',obj.time.steadyconverge);
            obj.time.cfl = extractInputRobust(ROMtext,'cfl',obj.time.cfl);
            obj.time.genalpha = extractInputRobust(ROMtext,'genalpha',[]);
            
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
            TimeSchemeString = extractInputROM(ROMtext,N,obj.id,'TimeScheme',[]);
            obj.TimeScheme     = determineTimeScheme(TimeSchemeString,obj);
            lastTimeSchemeString = extractInputROM(ROMtext,N,obj.id,'lastTimeScheme',[]);
            if ~isempty(lastTimeSchemeString), obj.lastTimeScheme = determineTimeScheme(lastTimeSchemeString,obj); end;
            
            %%%%%%%%% Determine newton structure %%%%%%%%%
            obj.newt = struct('maxIter',[],'eps',[],'quiet',[],'iter',[],'linesrch',[],'converge',[],'avgIter',[]);
            obj.newt.maxIter = extractInputROM(ROMtext,N,obj.id,'maxIter',10);
            obj.newt.eps     = extractInputROM(ROMtext,N,obj.id,'eps',1e-5);
            obj.newt.quiet   = extractInputROM(ROMtext,N,obj.id,'newtQuiet',false);
            obj.newt.iter    = zeros(1,obj.time.nstep);
            %Determine linesearch status and properties
            obj.newt.linesrch.status = extractInputROM(ROMtext,N,obj.id,'linesrch',false);
            obj.newt.linesrch.prop   = extractInputROM(ROMtext,N,obj.id,'linesrchProp',[1e-4,0.9,50,5,10,0.05,0.1,1.2]);
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
            
            %%%%%%%%% Initialize SVD update structure and determine buffer %%%%%%%%%
            obj.basisUpdate = extractInputROM(ROMtext,N,obj.id,'basisUpdate','border_fast');
            
            if strcmpi(obj.basisUpdate,'border')
                obj.svdUpdateData = repmat(struct('bufferSize',[],'wRef',[],'ColSumRSVtrunc',[],'normOneMinusVVt',[],'bufferedLSingVecs',[]),obj.nBases,1);
            elseif strcmpi(obj.basisUpdate,'border_fast')
                obj.svdUpdateData = repmat(struct('bufferSize',[],'wRef',[],'ColSumRSVtrunc',[],'normOneMinusVVt',[],'bufferedLSingVecs',[],...
                    'a',[],'b',[],'c',[],'d',[],'e',[],'F',[]),obj.nBases,1);
            elseif strcmpi(obj.basisUpdate,'border_exact')
                obj.svdUpdateData = repmat(struct('snapMAT',[],'prevRef',[]),obj.nBases,1);
            end
            
            if strcmpi(obj.basisUpdate,'border') || strcmpi(obj.basisUpdate,'border_fast') || strcmpi(obj.basisUpdate,'border_fast_approx')
                tmp = extractInputROM(ROMtext,N,obj.id,'svdUpdateBuffer',0); tmp = tmp(:);
                %Adjust size/check for compatibility
                if size(tmp,1) == 1
                    tmp = tmp*ones(obj.nBases,1);
                elseif size(tmp,1) ~= obj.nBases
                    error(['Size of svdUpdateBuffer not compatible with nBases in id ',num2str(obj.id)]);
                end
                %Load values into svdUpdateData structure
                for i = 1:obj.nBases
                    obj.svdUpdateData(i).bufferSize = tmp(i);
                end
            else
                for i = 1:obj.nBases
                    obj.svdUpdateData(i).bufferSize = 0;
                end
            end
            %%%%%%%%% Determine snaps structure %%%%%%%%%
            obj.snaps = struct('initref',[],'int',[],'num',[],'distseed',[],'dist',[],'distparams',[]);
            obj.snaps.initref    = extractInputROM(ROMtext,N,obj.id,'initsnapref','ic');
            obj.snaps.int        = extractInputROM(ROMtext,N,obj.id,'snapInt',[]);
            obj.snaps.num        = extractInputROM(ROMtext,N,obj.id,'nsnap','all');
            snapDist             = extractInputROM(ROMtext,N,obj.id,'snapDist',{[],[],[]});
            if isempty(snapDist)
                obj.snaps.distseed = [];
                obj.snaps.dist = [];
                obj.snaps.distparams = [];
            else
                obj.snaps.distseed   = snapDist{1};
                obj.snaps.dist       = snapDist{2};
                obj.snaps.distparams = snapDist{3};
            end
        end
        
        %Internal Functions
        function  [snapMat,ind] = extractRawSnapshots(obj,fomobj,flag,indices)
            %This function determines the state vector snapshot matrix from
            %FOM object input while adhereing to the user requests
            %indicated in the snaps structure
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - ROM object
            %fomobj  - FOM object to extract state vectors from
            %trainIC - ndof x 1 vector defining the training initial
            %          condition.  This is to be used if multiple
            %          simulations comprise the snapshots for a Local POD
            %          simulation.  If this is empty, the initial condition
            %          of the fomobj is used.
            %flag    - boolean = 1 if you wish to extract snapshots
            %          according to the CFG file and = 0 if you wish to
            %          extract all snapshots (in the 0 case, static
            %          simulaitons will draw from either fom.svIter or
            %          fom.sv (depending on availability) and dynamic
            %          simulations will always return fom.sv(:,2:end)
            %
            %Outputs:
            %--------
            %snapMat - snapshot matrix
            %--------------------------------------------------------------
            
            if size(fomobj.sv,2) == 2
                ind=2;
                snapMat = fomobj.sv(:,2);
                return;
            end
            
            %Handle the steady state case.
            if strcmpi(class(fomobj.TimeScheme),'SteadyTime')
                if fomobj.saveAllIt
                    fprintf('Since using "SteadyTime" integration and saveAllIt = true, snapshots taken from svIter\n');
                    snapMat = fomobj.svIter;
                    ind = 1:size(fomobj.svIter,2);
                else
                    snapMat = fomobj.sv(:,2:end);
                    ind = 2:size(fomobj.sv,2);
                end
                return;
            end
            
            %Determine the all possible snapshots from the state vector
            %matrix without considering the snapshot distribution,
            %interval, etc.
            if ~flag
                if nargin == 4
                    snapMat = fomobj.sv(:,indices);
                    ind=indices;
                else
                    snapMat = fomobj.sv(:,2:end);
                    ind = 2:size(fomobj.sv,2);
                end
                return;
            end
            
            Fsv = fomobj.sv(:,2:end);
            %Determine the time vector from the fomobj time structure
            %and remove the first time (corresponds to ic)
            t = linspace(fomobj.time.T(1),fomobj.time.T(2),fomobj.time.nstep+1);
            t(:,1) = [];
            
            %Determine the indicies of t that are valid for use as
            %snapshots (i.e. in the snapshot interval) and the maximum
            %number of snapshots in the interval specified.
            if ~isempty(obj.snaps.int) && size(obj.snaps.int,2) == 2
                indValid = find((t >= obj.snaps.int(1,1)) & (t <= obj.snaps.int(1,2)));
                nsnapMAX = numel(indValid);
            else
                %if obj.snaps.int is empty, the entire time domain is used
                indValid = 1:fomobj.time.nstep;
                nsnapMAX = size(fomobj.sv,2)-1;
            end
            
            %Determine the actual number of snapshots that will be used
            if isempty(obj.snaps.num) || (ischar(obj.snaps.num) && strcmpi(obj.snaps.num,'all'))
                %If obj.snaps.num is a char array whose value is 'all',
                %use all snapshots in the indicated interval
                numsnaps = nsnapMAX;
            elseif ischar(obj.snaps.num) && ~strcmpi(obj.snaps.num,'all')
                %If obj.snaps.num is a char array whose value is not 'all',
                %return an error
                error('Only valid keyword in the nsnap property is all');
            else
                %Otherwise, set the number of snapshots to the number
                %indicated in the input file (as long as obj.snaps.num
                %< nstep)
                if obj.snaps.num > fomobj.time.nstep
                    fprintf('Warning:  Too many snapshots specified.  Using nsnapMAX = %d snapshots',nsnapMAX);
                    numsnaps = nsnapMAX;
                elseif obj.snaps.num < 1 %percentage!
                    numsnaps = ceil(obj.snaps.num*nsnapMAX);
                else
                    numsnaps = obj.snaps.num;
                end
            end
            
            %Determine the indices (in Fsv) of the snapshots that will be
            %used.
            if numsnaps == nsnapMAX
                %If all snapshots in the interval are to be used, simply
                %return those indices.
                ind = indValid;
            else
                %Otherwise, compute them from the specified
                %distribution.
                if ~isempty(obj.snaps.dist) && ~isempty(obj.snaps.distseed) && ~isempty(obj.snaps.distparams)
                    ind = computeSnapshotFromDist(numsnaps,obj.snaps.distparams,...
                        obj.snaps.dist,obj.snaps.distseed,indValid(1),indValid(end));
                else
                    ind = computeSnapshotFromDist(numsnaps,indValid(end),'uniform',1,[],[]);
                end
            end
            %Return the appropriate snapshots
            snapMat = Fsv(:,ind);
            ind = ind+1; %Add 1 because we want to take into account the fact that we skipped over the initial condition.
        end %Done
        
        function  [] = NewtonRaphson(obj)
            %This function solves systems of the form:  f(x) = 0 using
            %Newton-Raphson's method
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - ROM object
            %
            %Outputs:
            %--------
            %No outputs. The state vector at the current time
            %(indicated by the cTimeIter property) is stored in the FOM
            %handle class.
            %--------------------------------------------------------------
            
            %Determine time iteration number plus 1
            itnump1 = obj.cTimeIter + 1;
            
            %%%%%%UPDATE BASIS ON THE FLY%%%%%%
            %             if ( itnump1-1 == 1 && ~strcmpi(obj.snaps.initref,'ic'))
            %                 disp('Updating state ROB at initial time');
            %                 %If  the initial snapshot reference was NOT the
            %                 %initial condition, perform thin SVD update
            %                 n = obj.nY+obj.svdUpdateData.bufferSize;
            %
            %                 vec = obj.svdUpdateData.wRef - obj.sv(:,itnump1-1);
            %                 U = ThinSVD_AddVec2Cols_Fast([obj.phi,obj.svdUpdateData.bufferedLSingVecs],...
            %                             obj.SingVals(1:n,1),obj.svdUpdateData.ColSumRSVtrunc,...
            %                             obj.svdUpdateData.normOneMinusVVt,vec);
            %                 obj.phi = U(:,1:obj.nY);
            %             end
            
            %Determine the residual and jacobian based on initial guess
            t = obj.time.T(1) + obj.time.dt*obj.cTimeIter;
            
            [R,J] = obj.TimeScheme.TimeIntNLFunc(obj.prob,obj.sv(:,itnump1-1),obj.sv(:,itnump1-1),t);
            %Determine convergence criteria to use
            %             Hres=(J*obj.phi).'*(J*obj.phi);
            res0 = norm(R,2);
            [tol,tolIt] = determineConvergeCriterion(obj,res0);
            
            indexAdj=1;
            if obj.printLevel > 1
                fprintf('---- Newton Step    # --------- ||R|| ------------ ||W''*R|| ----------- ||du|| ------------- tol -------------- tol_du ------ snap_coll ----------\n');
            end
            
            
            if obj.cTimeIter==1
                obj.gOrig=[];
            end
            
            
            for i_N = 1:obj.newt.maxIter
                %Solve for the search direction using Newton or Gauss-Newton
                if lower(obj.Rtype(1)) == 'g'
                    newR = obj.phi'*R;
                    du = -((obj.phi'*J*obj.phi)\newR);
                    conv = norm(newR,2);
                    %conv(i_N) = norm(newR,2);
                elseif strcmpi(obj.Rtype,'Petrov-Galerkin')
                    newJ = J*obj.phi;
                    if obj.augment
                        if i_N == 1
                            S_k = zeros(obj.nY);
                            du = -(newJ\R);
                        else
                            [Ro,Jo] = obj.TimeScheme.TimeIntNLFunc(obj.prob,obj.sv(:,itnump1)-p,obj.sv(:,itnump1-1),t);
                            S_k = obj.GN_HessianMod(S_k,du,newJ,R,Jo*obj.phi,Ro);
                            du = -((newJ'*newJ + S_k)\(newJ'*R));
                        end
                    else
                        du = -(newJ\R);
                    end
                    
                    conv = norm(newJ'*R,2);
                elseif strcmpi(obj.Rtype,'Petrov-Galerkin-Regularized')
                    %                     mu    = 0.1;
                    %                     sigma = 1;
                    
                    newJ = [ J*obj.phi; obj.mu*eye(obj.nY)] ;
                    newR = [ obj.sigma*R ; zeros(obj.nY,1) ];
                    du = -(newJ\newR);
                    conv = norm(newJ'*newR,2);
                    %conv(i_N) = norm(newJ'*R,2);
                end
                
                if i_N == 1, [tol,tolIt] = determineConvergeCriterion(obj,conv); end;
                
                %If nan or complex encountered -> kill simulation
                if sum(isnan(du)) || ~isreal(du)
                    obj.killflag=true;
                    return;
                end
                
                %Linesearch if necessary
                if obj.newt.linesrch.status
                    %alpha = linesearch(@(alpha) obj.TimeScheme.TimeIntNLFuncDir(obj.prob,obj.sv(:,itnump1-indexAdj),obj.sv(:,itnump1-1),t,alpha,du,obj),obj.newt.linesrch.prop);
                    alpha = linesrchBackNewton(@(alpha) obj.TimeScheme.TimeIntNLFuncDir(obj.prob,obj.sv(:,itnump1-indexAdj),obj.sv(:,itnump1-1),t,alpha,du,obj));
                    du = alpha*du;
                end
                %Reconstructed search direction
                p = obj.phi*du;
                
                obj.sv(:,itnump1)=obj.sv(:,itnump1-indexAdj)+p;
                
                if obj.printLevel > 2
                    fprintf('---- Newton Step %4i ----- %10.7e ----- %10.7e ----- %10.7e ----- %10.7e ----- %10.7e ------- %1i --------------\n',i_N,norm(R,2),conv,norm(du,2),tol,tolIt,obj.saveNL);
                end
                
                %WHY DO WE SKIP THE FIRST ONE?
                if i_N > 1 || obj.newt.maxIter == 1
                    %Write nonlinear snapshots
                    obj.writeNonlinearSnapshot(R,J,p,1);
                end
                
                indexAdj=0;
                if checkConverge(obj,conv,du,tol,tolIt)
                    break;
                end
                
                %Compute residual and jacobian with updated vector for next iteration
                [R,J] = obj.TimeScheme.TimeIntNLFunc(obj.prob,obj.sv(:,itnump1),obj.sv(:,itnump1-1),t);
                %                 Hres=(J*obj.phi).'*(J*obj.phi);
            end
            
            w_guess=obj.phi'*(obj.sv(:,itnump1)-obj.sv(:,obj.cTimeIter));
            g(1,:)=obj.constr1(w_guess);
            g(2,:)=obj.constr2(w_guess);
            g(3,:)=obj.constr3(w_guess);
            
            obj.gOrig=[obj.gOrig g];
            Corig=obj.gOrig;
            save constrOrig Corig
            
            
            %             if itnump1==2
            %                 y_NR=obj.sv(:,itnump1);
            %                 save NR y_NR;
            %                 g1=obj.constr1(w_guess)
            %                 g2=obj.constr2(w_guess)
            %                 g3=obj.constr3(w_guess)
            %             end
            %Write the last nonlinear snapshot.  We note that we are using
            %J^(k+1)*0 if using Snapshot 1.5 or 2 (we don't have a
            %search direction p^(k+1) because we converged!).  It shouldn't
            %matter much because ||p|| should be small since we converged,
            %i.e. the next step shouldn't take us away from the solution!
            %obj.writeNonlinearSnapshot(R,J,zeros(size(p)),1);
            
            %Store the number of newton iterations at this time step
            obj.newt.iter(obj.cTimeIter) = i_N;
            
            if obj.printLevel == 1.5
                fprintf('---- Newton Step %4i ----- %10.7e ----- %10.7e ----- %10.7e ----- %10.7e ----- %10.7e ------- %1i --------------\n',i_N,norm(R,2),conv,norm(p,2),tol,tolIt,obj.saveNL);
            end
            
            %If the maxiumum newton iterations are reached, warn the user
            if ~obj.newt.quiet && (i_N==obj.newt.maxIter)
                disp(['*** Warning; Newton solver reached max number of iterations before convergence : res = ',...
                    num2str(norm(R,2)/tol),'  eps = ',num2str(obj.newt.eps), '***']);
            end
        end %Done
        
        function  [] = NewtonRaphson_trunc(obj)
            %This function solves systems of the form:  f(x) = 0 using
            
            %Newton-Raphson's method
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - ROM object
            %
            %Outputs:
            %--------
            %No outputs. The state vector at the current time
            %(indicated by the cTimeIter property) is stored in the FOM
            %handle class.
            %--------------------------------------------------------------
            
            %Determine time iteration number plus 1
            itnump1 = obj.cTimeIter + 1;
            
            %%%%%%UPDATE BASIS ON THE FLY%%%%%%
            %             if ( itnump1-1 == 1 && ~strcmpi(obj.snaps.initref,'ic'))
            %                 disp('Updating state ROB at initial time');
            %                 %If  the initial snapshot reference was NOT the
            %                 %initial condition, perform thin SVD update
            %                 n = obj.nY+obj.svdUpdateData.bufferSize;
            %
            %                 vec = obj.svdUpdateData.wRef - obj.sv(:,itnump1-1);
            %                 U = ThinSVD_AddVec2Cols_Fast([obj.phi,obj.svdUpdateData.bufferedLSingVecs],...
            %                             obj.SingVals(1:n,1),obj.svdUpdateData.ColSumRSVtrunc,...
            %                             obj.svdUpdateData.normOneMinusVVt,vec);
            %                 obj.phi = U(:,1:obj.nY);
            %             end
            
            %Determine the residual and jacobian based on initial guess
            t = obj.time.T(1) + obj.time.dt*obj.cTimeIter;
            Phi=obj.phi(:, 1:obj.trunc);
            
            [R,J] = obj.TimeScheme.TimeIntNLFunc(obj.prob,obj.sv(:,itnump1-1),obj.sv(:,itnump1-1),t);
            %Determine convergence criteria to use
            % Hres=(J*obj.phi).'*(J*obj.phi);
            res0 = norm(R,2);
            [tol,tolIt] = determineConvergeCriterion(obj,res0);
            
            indexAdj=1;
            if obj.printLevel > 1
                fprintf('---- Newton Step    # --------- ||R|| ------------ ||W''*R|| ----------- ||du|| ------------- tol -------------- tol_du ------ snap_coll ----------\n');
            end
            
            if obj.cTimeIter==1
                obj.gOrig=[];
                size(Phi)
            end
            
            for i_N = 1:obj.newt.maxIter
                %Solve for the search direction using Newton or Gauss-Newton
                if lower(obj.Rtype(1)) == 'g'
                    newR = Phi'*R;
                    du = -((Phi'*J*Phi)\newR);
                    conv = norm(newR,2);
                    %conv(i_N) = norm(newR,2);
                    if obj.cTimeIter==1 &i_N==1
                        disp(['Galerkin'])
                    end
                elseif strcmpi(obj.Rtype,'Petrov-Galerkin')
                    
                    if obj.cTimeIter==1 &i_N==1
                        disp(['Petrov-Galerkin'])
                    end
                    newJ = J*Phi;
                    if obj.augment
                        if i_N == 1
                            S_k = zeros(obj.nY);
                            du = -(newJ\R);
                        else
                            [Ro,Jo] = obj.TimeScheme.TimeIntNLFunc(obj.prob,obj.sv(:,itnump1)-p,obj.sv(:,itnump1-1),t);
                            S_k = obj.GN_HessianMod(S_k,du,newJ,R,Jo*Phi,Ro);
                            du = -((newJ'*newJ + S_k)\(newJ'*R));
                        end
                    else
                        du = -(newJ\R);
                    end
                    
                    conv = norm(newJ'*R,2);
                elseif strcmpi(obj.Rtype,'Petrov-Galerkin-Regularized')
                    % mu    = 0.1;
                    % sigma = 1;
                    
                    newJ = [ J*Phi; obj.mu*eye(obj.nY)] ;
                    newR = [ obj.sigma*R ; zeros(obj.nY,1) ];
                    du = -(newJ\newR);
                    conv = norm(newJ'*newR,2);
                    %conv(i_N) = norm(newJ'*R,2);
                end
                
                if i_N == 1, [tol,tolIt] = determineConvergeCriterion(obj,conv); end;
                
                %If nan or complex encountered -> kill simulation
                if sum(isnan(du)) || ~isreal(du)
                    obj.killflag=true;
                    return;
                end
                
                %Linesearch if necessary
                if obj.newt.linesrch.status
                    %alpha = linesearch(@(alpha) obj.TimeScheme.TimeIntNLFuncDir(obj.prob,obj.sv(:,itnump1-indexAdj),obj.sv(:,itnump1-1),t,alpha,du,obj),obj.newt.linesrch.prop);
                    alpha = linesrchBackNewton(@(alpha) obj.TimeScheme.TimeIntNLFuncDir(obj.prob,obj.sv(:,itnump1-indexAdj),obj.sv(:,itnump1-1),t,alpha,du,obj));
                    du = alpha*du;
                end
                %Reconstructed search direction
                p = Phi*du;
                
                obj.sv(:,itnump1)=obj.sv(:,itnump1-indexAdj)+p;
                
                if obj.printLevel > 2
                    fprintf('---- Newton Step %4i ----- %10.7e ----- %10.7e ----- %10.7e ----- %10.7e ----- %10.7e ------- %1i --------------\n',i_N,norm(R,2),conv,norm(du,2),tol,tolIt,obj.saveNL);
                end
                
                %WHY DO WE SKIP THE FIRST ONE?
                if i_N > 1 || obj.newt.maxIter == 1
                    %Write nonlinear snapshots
                    obj.writeNonlinearSnapshot(R,J,p,1);
                end
                
                indexAdj=0;
                if checkConverge(obj,conv,du,tol,tolIt)
                    break;
                end
                
                %Compute residual and jacobian with updated vector for next iteration
                [R,J] = obj.TimeScheme.TimeIntNLFunc(obj.prob,obj.sv(:,itnump1),obj.sv(:,itnump1-1),t);
                %                 Hres=(J*Phi).'*(J*obj.phi);
            end
            
            w_guess=Phi'*(obj.sv(:,itnump1)-obj.sv(:,obj.cTimeIter));
            
            %Write the last nonlinear snapshot.  We note that we are using
            %J^(k+1)*0 if using Snapshot 1.5 or 2 (we don't have a
            %search direction p^(k+1) because we converged!).  It shouldn't
            %matter much because ||p|| should be small since we converged,
            %i.e. the next step shouldn't take us away from the solution!
            %obj.writeNonlinearSnapshot(R,J,zeros(size(p)),1);
            
            %Store the number of newton iterations at this time step
            obj.newt.iter(obj.cTimeIter) = i_N;
            
            if obj.printLevel == 1.5
                fprintf('---- Newton Step %4i ----- %10.7e ----- %10.7e ----- %10.7e ----- %10.7e ----- %10.7e ------- %1i --------------\n',i_N,norm(R,2),conv,norm(p,2),tol,tolIt,obj.saveNL);
                
            end
            
            %If the maxiumum newton iterations are reached, warn the user
            if ~obj.newt.quiet && (i_N==obj.newt.maxIter)
                disp(['*** Warning; Newton solver reached max number of iterations before convergence : res = ',...
                    num2str(norm(R,2)/tol),'  eps = ',num2str(obj.newt.eps), '***']);
            end
        end %Done
        
        
        
        function  [] = improvedLSPG(obj)
            
            itnump1 = obj.cTimeIter + 1;
            %             keyboard
            Phi=obj.phi(:,1:obj.trunc);
            
            t = obj.time.T(1) + obj.time.dt*obj.cTimeIter;
            obj.sv(:,itnump1)=obj.sv(:,obj.cTimeIter);
            
            w_guess = Phi'* (obj.sv(:,itnump1)-obj.sv(:,obj.cTimeIter)); %initial guess for generalized coordinates for Newton iteration
            
            if obj.printLevel > 1
                fprintf('---- Newton Step    # --------- ||R|| ------------ ||du|| ------------- \n');
            end
            
            if obj.cTimeIter==1
                obj.constr=[];
            end
            
            for k=1:20 %obj.newt.maxIter
                
                if obj.cTimeIter==1 & k==1
                    size(Phi)
                    disp(['improved number of cells  ', num2str(obj.ncell)])
                end
                obj.sv(:,itnump1)=obj.sv(:,obj.cTimeIter)+Phi*w_guess;
                [Res,Jres]=obj.TimeScheme.TimeIntNLFunc(obj.prob,obj.sv(:,itnump1),obj.sv(:,itnump1-1),t);
                Df=Jres*Phi;
                
                [g,Dg]=obj.constraintsOneDomain_old(w_guess);
                %                 [g_check,Dg_check]=obj.constraintsForGNAT(w_guess, obj.sv(:,itnump1), obj.time.dt);
                %                 if  norm(g-g_check)>1e-5 || norm(Dg-Dg_check)>1e-5
                %                     disp('constraints are different')
                %                     keyboard
                %                 end
                %
                %                 params.DifferenceStep = 1e-6;
                %                 params.DifferenceType='centered';
                %                 for ii=1:3
                %                     out=gradientcheck(@(w) obj.testConst_old(w, ii), w_guess, params);%poblano toolbox
                %                     if out.RelError>1e-4, keyboard, end
                %                 end
                
                H=Df'*Df; h=Df'*Res;
                P=H\Dg';
                x=H\h;
                Sm=-inv(Dg*P);
                Qzim=P*Sm;
                del_w=Qzim*g-(x+Qzim*Dg*x);
                
                w_guess=w_guess+del_w;
                
                if norm(del_w,2)<10^(-5)
                    break;
                end
            end
            [g,~]=obj.constraintsOneDomain_old(w_guess);
            obj.constr=[obj.constr,g];
            %             constr=obj.constr;
            %             save Zimconstr9999 constr
            
            disp(['norm of the constraint   ', num2str(norm(g))])
            disp(['number of iterations in the improvedLSPG  ',num2str(k)])
        end
        
        
        function [g, Dg]=constraintsOneDomain_old(obj, w_increment)
            
            [rFlux, drFlux]=obj.RightFlux_old(w_increment);
            [lFlux, dlFlux]=obj.LeftFlux_old(w_increment);
            [sumQ,sumdQ]= obj.sumForceTerms_old(w_increment);
            
            Phi1=obj.phi(1:3:end,:); %basis for the first conserved quantity rho
            Phi2=obj.phi(2:3:end,:); %basis for the second conserved quantity rho*u
            Phi3=obj.phi(3:3:end,:); %basis for the thirs conserved quantity e
            
            
            g=[(obj.prob.SVol(2:end-1).*obj.prob.dx(2:end-1))*Phi1(2:end-1,:)*w_increment;
                (obj.prob.SVol(2:end-1).*obj.prob.dx(2:end-1))*Phi2(2:end-1,:)*w_increment;
                (obj.prob.SVol(2:end-1).*obj.prob.dx(2:end-1))*Phi3(2:end-1,:)*w_increment]+...
                obj.time.dt*(rFlux+lFlux-sumQ);
            
            Dg= [(obj.prob.SVol(2:end-1).*obj.prob.dx(2:end-1))*Phi1(2:end-1,:);
                (obj.prob.SVol(2:end-1).*obj.prob.dx(2:end-1))*Phi2(2:end-1,:);
                (obj.prob.SVol(2:end-1).*obj.prob.dx(2:end-1))*Phi3(2:end-1,:)]+...
                obj.time.dt*(dlFlux+ drFlux- sumdQ);
            
        end
        
        function   [r,dr]=RightFlux_old(obj, w_increment)
            Phi=obj.phi(:,1:obj.trunc);          
            SV=obj.sv(:,obj.cTimeIter)+Phi(:,1:obj.trunc)*w_increment;
            [rho, u, P,c,e,~,dc]=obj.prob.getVariables(SV);
            [roeF, droeF]=obj.prob.roeFlux(rho,u,P,c,e,dc);
            dUdV=[1,0,0;u(1),rho(1),0;0.5*u(1)*u(1),rho(1)*u(1),1/(obj.prob.gamma-1)];
            droeF(:,1:3,1)=droeF(:,1:3,1)*dUdV;
            dUdV=[1,0,0;u(end),rho(end),0;0.5*u(end)*u(end),rho(end)*u(end),1/(obj.prob.gamma-1)];
            droeF(:,4:6,end)=droeF(:,4:6,end)*dUdV;
            
            % governEqn
            Rright= bsxfun(@times,roeF(:,2:end),obj.prob.S(3:end-1));
            r=Rright(:,end);
            finish=obj.prob.nVol-1;
            derJR=obj.prob.S(finish+1)*droeF(:,:,finish);
            dr=derJR*Phi(3*(finish-1)+1:3*(finish+1),:);
            %             if obj.cTimeIter==1
            %                 keyboard
            %                 dr_true=dr; save dr_true dr_true;
            %             end
        end
        
        
        function   [l,dl]=LeftFlux_old(obj, w_increment)
            
            Phi = obj.phi(:,1:obj.trunc);
            SV=obj.sv(:,obj.cTimeIter) + Phi * w_increment;
            
            [rho, u, P,c,e,~,dc]=obj.prob.getVariables(SV);
            [roeF, droeF]=obj.prob.roeFlux(rho,u,P,c,e,dc);
            dUdV=[1,0,0;u(1),rho(1),0;0.5*u(1)*u(1),rho(1)*u(1),1/(obj.prob.gamma-1)];
            droeF(:,1:3,1)=droeF(:,1:3,1)*dUdV;
            dUdV=[1,0,0;u(end),rho(end),0;0.5*u(end)*u(end),rho(end)*u(end),1/(obj.prob.gamma-1)];
            droeF(:,4:6,end)=droeF(:,4:6,end)*dUdV;
            
            Rleft= -bsxfun(@times,roeF(:,1:end-1),obj.prob.S(2:end-2));
            l=Rleft(:,1);
            start=1;
            derJL=-obj.prob.S(start+1)*droeF(:,:,start);
            dl=derJL*Phi(3*(start-1)+1:3*(start+1),:);
            
        end
        
        function [sq, sdq]=sumForceTerms_old(obj, w_increment)
            Phi=obj.phi(:,1:obj.trunc);
            SV=obj.sv(:,obj.cTimeIter)+Phi*w_increment;
            [~, u, P,~,~,~,~]=obj.prob.getVariables(SV);
            [Q,dQ]=obj.prob.forceTerm(u,P);
            %keyboard
            sq=Q(:,2:end-1)*(obj.prob.SVol(2:end-1).*obj.prob.dx(2:end-1))';
            sdq=zeros(3,obj.trunc);
            for i=2:obj.prob.nVol-1
                %            derivQ=derivQ+squeeze(dQ(:,:,i))*[Phi1(i,:);Phi2(i,:);Phi3(i,:)]*(obj.prob.SVol(i).*obj.prob.dx(i));
                sdq=sdq+squeeze(dQ(:,:,i))*Phi(3*i-2:3*i,:)*(obj.prob.SVol(i).*obj.prob.dx(i));
            end

        end
        
        function [ci,dci]=testConst_old(obj, w_increment,ind)
            [g,Dg]=obj.constraintsOneDomain_old(w_increment);
            ci=g(ind);
            dci=Dg(ind,:)';
        end
        
        
        
        function [g, Dg]=constraintsForGNAT(obj, w_increment,sv,delt)
            
            % keyboard
            [rFlux, drFlux]=obj.RightFlux(sv);
            [lFlux, dlFlux]=obj.LeftFlux(sv);
            [sumQ,sumdQ]= obj.sumForceTerms(sv);
            
            Phi1=obj.phi(1:3:end,:); %basis for the first conserved quantity rho
            Phi2=obj.phi(2:3:end,:); %basis for the second conserved quantity rho*u
            Phi3=obj.phi(3:3:end,:); %basis for the thirs conserved quantity e
            
            
            g=[(obj.prob.SVol(2:end-1).*obj.prob.dx(2:end-1))*Phi1(2:end-1,:)*w_increment;
                (obj.prob.SVol(2:end-1).*obj.prob.dx(2:end-1))*Phi2(2:end-1,:)*w_increment;
                (obj.prob.SVol(2:end-1).*obj.prob.dx(2:end-1))*Phi3(2:end-1,:)*w_increment]+...
                delt*(rFlux+lFlux-sumQ);
            
            Dg= [(obj.prob.SVol(2:end-1).*obj.prob.dx(2:end-1))*Phi1(2:end-1,:);
                (obj.prob.SVol(2:end-1).*obj.prob.dx(2:end-1))*Phi2(2:end-1,:);
                (obj.prob.SVol(2:end-1).*obj.prob.dx(2:end-1))*Phi3(2:end-1,:)]+...
                delt*(dlFlux+ drFlux- sumdQ);
            
        end
        
        function   [r,dr]=RightFlux(obj, SV)
            
            [rho, u, P,c,e,~,dc]=obj.prob.getVariables(SV);
            [roeF, droeF]=obj.prob.roeFlux(rho,u,P,c,e,dc);
            dUdV=[1,0,0;u(1),rho(1),0;0.5*u(1)*u(1),rho(1)*u(1),1/(obj.prob.gamma-1)];
            droeF(:,1:3,1)=droeF(:,1:3,1)*dUdV;
            dUdV=[1,0,0;u(end),rho(end),0;0.5*u(end)*u(end),rho(end)*u(end),1/(obj.prob.gamma-1)];
            droeF(:,4:6,end)=droeF(:,4:6,end)*dUdV;
            
            % code from fucntion governEqn
            Rright= bsxfun(@times,roeF(:,2:end),obj.prob.S(3:end-1));
            r=Rright(:,end);
            finish=obj.prob.nVol-1;
            derJR=obj.prob.S(finish+1)*droeF(:,:,finish);
            dr=derJR*obj.phi(3*(finish-1)+1:3*(finish+1),:);
        end
        
        function   [l,dl]=LeftFlux(obj, SV)
            
            [rho, u, P,c,e,~,dc]=obj.prob.getVariables(SV);
            [roeF, droeF]=obj.prob.roeFlux(rho,u,P,c,e,dc);
            dUdV=[1,0,0;u(1),rho(1),0;0.5*u(1)*u(1),rho(1)*u(1),1/(obj.prob.gamma-1)];
            droeF(:,1:3,1)=droeF(:,1:3,1)*dUdV;
            dUdV=[1,0,0;u(end),rho(end),0;0.5*u(end)*u(end),rho(end)*u(end),1/(obj.prob.gamma-1)];
            droeF(:,4:6,end)=droeF(:,4:6,end)*dUdV;
            
            Rleft= -bsxfun(@times,roeF(:,1:end-1),obj.prob.S(2:end-2));
            l=Rleft(:,1);
            start=1;
            derJL=-obj.prob.S(start+1)*droeF(:,:,start);
            dl=derJL*obj.phi(3*(start-1)+1:3*(start+1),:);
            
        end
        
        function [sq, sdq]=sumForceTerms(obj,SV)
            
            [~, u, P,~,~,~,~]=obj.prob.getVariables(SV);
            [Q,dQ]=obj.prob.forceTerm(u,P);
            sq=Q(:,2:end-1)*(obj.prob.SVol(2:end-1).*obj.prob.dx(2:end-1))';
            sdq=zeros(3,obj.trunc);
            for i=2:obj.prob.nVol-1
                %            derivQ=derivQ+squeeze(dQ(:,:,i))*[Phi1(i,:);Phi2(i,:);Phi3(i,:)]*(obj.prob.SVol(i).*obj.prob.dx(i));
                sdq=sdq+squeeze(dQ(:,:,i))*obj.phi(3*i-2:3*i,:)*(obj.prob.SVol(i).*obj.prob.dx(i));
            end
        end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 function r = RightFluxAll(obj, SV)
      
      [rho, u, P, c, e, ~, dc] = obj.prob.getVariables(SV);
      [roeF, ~ ] = obj.prob.roeFlux(rho,u,P,c,e,dc);
%       dUdV=[1,0,0;u(1),rho(1),0;0.5*u(1)*u(1),rho(1)*u(1),1/(obj.prob.gamma-1)];
%       droeF(:,1:3,1)=droeF(:,1:3,1)*dUdV;
%       dUdV=[1,0,0;u(end),rho(end),0;0.5*u(end)*u(end),rho(end)*u(end),1/(obj.prob.gamma-1)];
%       droeF(:,4:6,end)=droeF(:,4:6,end)*dUdV;

      % code from fucntion governEqn      
      r = reshape(bsxfun(@times,roeF(:,2:end),obj.prob.S(3:end-1)),...
          size(roeF,1)*size(roeF(:, 2:end),2),1);
     
  end
  
  function   l = LeftFluxAll(obj, SV)
      
      [rho, u, P, c, e, ~, dc] = obj.prob.getVariables(SV);
      [roeF, ~] = obj.prob.roeFlux(rho,u,P,c,e,dc);
%       dUdV=[1,0,0;u(1),rho(1),0;0.5*u(1)*u(1),rho(1)*u(1),1/(obj.prob.gamma-1)];
%       droeF(:,1:3,1)=droeF(:,1:3,1)*dUdV;
%       dUdV=[1,0,0;u(end),rho(end),0;0.5*u(end)*u(end),rho(end)*u(end),1/(obj.prob.gamma-1)];
%       droeF(:,4:6,end)=droeF(:,4:6,end)*dUdV;
      
      l = reshape(-bsxfun(@times, roeF(:,1:end-1), obj.prob.S(2:end-2)),...
          size(roeF, 1) * size(roeF(:, 2:end), 2), 1);
      
      
  end
  
  function [force_q,dforce_q] = ForceTermsAll(obj, SV)
      [~, u, P, ~,~,~,~]  = obj.prob.getVariables(SV);
      [force_q, dforce_q] = obj.prob.forceTerm(u,P);
  end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        
        function [ci,dci]=testConstGNAT(obj, w_increment,sv, ind)
            [g,Dg]=obj.constraintsOneDomain(w_increment,sv);
            ci=g(ind);
            dci=Dg(ind,:)';
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function chooseSolver(obj)
            
            % 	obj.trunc=min(find(cumsum(obj.SingVals)/sum(obj.SingVals)>0.9999));
            
            if obj.ncell==1
                if obj.newtonSolver==1
                    obj.NewtonRaphson_trunc();
                else
                    obj.improvedLSPG()
                end
            else
                if 3*obj.ncell<obj.trunc
                    obj.LSPG_subdomains();
                elseif 3*obj.ncell==obj.trunc
                    obj.solveConstraints();
                elseif 3*obj.ncell>obj.trunc
                    obj.minimizeConstraints();
                end
            end
        end
        
        
        
        
        function  [] = LSPG_subdomains(obj)
            %%use if 3*obj.ncell < obj.trunc
            itnump1 = obj.cTimeIter + 1;
            Phi=obj.phi(:,1:obj.trunc);
            t = obj.time.T(1) + obj.time.dt*obj.cTimeIter;
            obj.sv(:,itnump1)=obj.sv(:,obj.cTimeIter);
            
            w_guess = Phi'* (obj.sv(:,itnump1)-obj.sv(:,obj.cTimeIter));
            if obj.printLevel > 1
                fprintf('---- Newton Step    # --------- ||R|| ------------ ||du|| ------------- \n');
            end
            
            for k=1:obj.newt.maxIter
                
                if obj.cTimeIter==1 & k==1
                    size(Phi)
                    disp(['LSPG, constraints on several domaines; num of domains ', num2str(obj.ncell)])
                end
                obj.sv(:,itnump1)=obj.sv(:,obj.cTimeIter)+Phi*w_guess;
                [Res,Jres]=obj.TimeScheme.TimeIntNLFunc(obj.prob,obj.sv(:,itnump1),obj.sv(:,itnump1-1),t);
                Df=Jres*Phi;
                %keyboard
                [g,Dg]=obj.myconstraints(w_guess);
                for i=1:length(g)
                    out=gradientcheck( @(w) obj.testmyconstraints(w,i), w_guess);
                    if out.RelError>1e-4,keyboard, end
                end
                
                H=Df'*Df; h=Df'*Res;
                P=H\Dg';
                x=H\h;
                S=-inv(Dg*P);
                Qzim=P*S;
                del_w=Qzim*g-(x+Qzim*Dg*x);
                
                
                w_guess=w_guess+del_w;
                
                if norm(del_w,2)<10^(-6)
                    break;
                end
            end
                    
            [g,~]=obj.myconstraints(w_guess);
            disp(['norm of the constraints ', num2str(norm(g))])
            %disp(['number of iterations in truncated case ',num2str(k)])
        end
        
%   function [cnstr, Dcnstr]=myconstraints(obj, w_increment)
%             
%             dn=floor(obj.prob.nVol/(obj.ncell));
%             points=1:dn:obj.prob.nVol+1;
%             if points(end)~=obj.prob.nVol+1 && points(end)~=obj.prob.nVol
%                 points(end+1)=obj.prob.nVol+1;
%             else
%                 points(end)= obj.prob.nVol+1;
%             end
%             
%             Phi=obj.phi(:,1:obj.trunc);
%             Phi1=Phi(1:3:end,:); %basis for the first conserved quantity rho
%             Phi2=Phi(2:3:end,:); %basis for the second conserved quantity rho*u
%             Phi3=Phi(3:3:end,:); %basis for the thirs conserved quantity e
%             
%             
%             SV=obj.sv(:,obj.cTimeIter)+Phi*w_increment;
%             [rho, u, P,c,e,~,dc]=obj.prob.getVariables(SV);
%             [Q,dQ]=obj.prob.forceTerm(u,P);
%             [roeF, droeF]=obj.prob.roeFlux(rho,u,P,c,e,dc);
%             
%             %from function governEqn 
%             dUdV=[1,0,0;u(1),rho(1),0;0.5*u(1)*u(1),rho(1)*u(1),1/(obj.prob.gamma-1)];
%             droeF(:,1:3,1)=droeF(:,1:3,1)*dUdV;
%             dUdV=[1,0,0;u(end),rho(end),0;0.5*u(end)*u(end),rho(end)*u(end),1/(obj.prob.gamma-1)];
%             droeF(:,4:6,end)=droeF(:,4:6,end)*dUdV;
%             
%             start=1;
%             finish=points(2)-1;
%             right=   roeF(:,finish)*obj.prob.S(finish+1);
%             left = - roeF(:,start)*obj.prob.S(start+1);
%             %keyboard
%             current_points = 2:points(2)-1;
%             cnstr(1:3,:)=[(obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi1(current_points,:)*w_increment;
%                 (obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi2(current_points,:)*w_increment;
%                 (obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi3(current_points,:)*w_increment ]+...
%                 obj.time.dt*(right+left-Q(:,current_points)*(obj.prob.SVol(current_points).*obj.prob.dx(current_points))');
%             
%             derivQ1=zeros(3,obj.trunc);
%             for i=2:points(2)-1 %1:points(2)-1
%                 derivQ1=derivQ1+squeeze(dQ(:,:,i))*[Phi1(i,:);Phi2(i,:);Phi3(i,:)]*(obj.prob.SVol(i).*obj.prob.dx(i));
%             end
%             
%             derJL= -obj.prob.S(start+1)*droeF(:,:,start); % 3x6 , first dimension  is w.r.t. (rho, rho*u, e); second dimension is w.r.t. cells j and j+1
%             derJR=  obj.prob.S(finish+1)*droeF(:,:,finish);
%             
%             Dcnstr(1:3,:)= [(obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi1(current_points,:);
%                 (obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi2(current_points,:);
%                 (obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi3(current_points,:)]+...
%                 obj.time.dt*(derJL*Phi(3*(start-1)+1:3*(start+1),:)+ ...
%                 derJR*Phi(3*(finish-1)+1:3*(finish+1),:)- derivQ1);
%             
%             
%             if obj.ncell>2
%                 for kk=2:obj.ncell-1
%                     start=points(kk)-1;
%                     finish=points(kk+1)-1;
%                     right=   roeF(:,finish)*obj.prob.S(finish+1);
%                     left = - roeF(:,start)*obj.prob.S(start+1);
%                     current_points=start+1:finish;
%                     cnstr(3*(kk-1)+1:3*kk,:)=[(obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi1(current_points,:)*w_increment;
%                         (obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi2(current_points,:)*w_increment;
%                         (obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi3(current_points,:)*w_increment]+...
%                         obj.time.dt*(right+left-Q(:,current_points)*(obj.prob.SVol(current_points).*obj.prob.dx(current_points))');
%                     %obj.time.dt*(right(:,points(kk+1)-2)+left(:,points(kk)-1)-Q(:,points(kk):points(kk+1)-1)*(obj.prob.SVol(points(kk):points(kk+1)-1).*obj.prob.dx(points(kk):points(kk+1)-1))');
%                     derivQ=zeros(3,obj.trunc);
%                     for ii=points(kk):points(kk+1)-1
%                         derivQ=derivQ+squeeze(dQ(:,:,ii))*[Phi1(ii,:);Phi2(ii,:);Phi3(ii,:)]*(obj.prob.SVol(ii).*obj.prob.dx(ii));
%                     end
%                     
%                     derJL= -obj.prob.S(start+1)*droeF(:,:,start); % 3x6 , first dimension  is w.r.t. (rho, rho*u, e); second dimension is w.r.t. cells j and j+1
%                     derJR= obj.prob.S(finish+1)*droeF(:,:,finish);
%                     
%                     Dcnstr(3*(kk-1)+1:3*kk,:)= [(obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi1(current_points,:);
%                         (obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi2(current_points,:);
%                         (obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi3(current_points,:)]+...
%                         obj.time.dt*(derJL*Phi(3*(start-1)+1:3*(start+1),:)+ ...
%                         derJR*Phi(3*(finish-1)+1:3*(finish+1),:)- derivQ);
%                     
%                 end
%             end
%             
%             
%             %keyboard
%             start=points(end-1)-1;
%             finish=points(end)-2;
%             right=  roeF(:,finish)*obj.prob.S(finish+1);
%             left = -roeF(:,start)*obj.prob.S(start+1);
%             current_points=start+1:finish;
%             cnstr((obj.ncell-1)*3+1:obj.ncell*3,:)=[(obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi1(current_points,:)*w_increment;
%                 (obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi2(current_points,:)*w_increment;
%                 (obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi3(current_points,:)*w_increment ]+...
%                 obj.time.dt*(right+left-Q(:,current_points)*(obj.prob.SVol(current_points).*obj.prob.dx(current_points))');            
%             
%             derivQ3=zeros(3,obj.trunc);
%             for ii= points(end-1):obj.prob.nVol-1 %points(end-1):obj.prob.nVol
%                 derivQ3=derivQ3+squeeze(dQ(:,:,ii))*[Phi1(ii,:);Phi2(ii,:);Phi3(ii,:)]*(obj.prob.SVol(ii).*obj.prob.dx(ii));
%             end
%             
%            
%             derJL=-obj.prob.S(start+1)* droeF(:,:,start); % 3x6 , first dimension  is w.r.t. (rho, rho*u, e); second dimension is w.r.t. cells j and j+1
%             derJR=obj.prob.S(finish+1)* droeF(:,:,finish);
%             
%             Dcnstr((obj.ncell-1)*3+1:obj.ncell*3,:)= [(obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi1(current_points,:);
%                 (obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi2(current_points,:);
%                 (obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi3(current_points,:)]+...
%                 obj.time.dt*(derJL*Phi(3*(start-1)+1:3*(start+1),:)+ ...
%                 derJR*Phi(3*(finish-1)+1:3*(finish+1),:)- derivQ3);
%                        
%         end
        
  %%%%
  function [cnstr, Dcnstr]=myconstraints(obj, w_increment)
                %keyboard
                dn=floor(obj.prob.nVol/(obj.ncell));
                points=1:dn:obj.prob.nVol+1;
                if points(end)~=obj.prob.nVol+1 && points(end)~=obj.prob.nVol
                        points(end+1)=obj.prob.nVol+1;
                else
                        points(end)= obj.prob.nVol+1;
                end

                Phi=obj.phi(:,1:obj.trunc);
                Phi1=Phi(1:3:end,:); %basis for the first conserved quantity rho
                Phi2=Phi(2:3:end,:); %basis for the second conserved quantity rho*u
                Phi3=Phi(3:3:end,:); %basis for the thirs conserved quantity e


                SV=obj.sv(:,obj.cTimeIter)+Phi*w_increment;
                [rho, u, P,c,e,~,dc]=obj.prob.getVariables(SV);
                [Q,dQ]=obj.prob.forceTerm(u,P);
                [roeF, droeF]=obj.prob.roeFlux(rho,u,P,c,e,dc);

                %from governEqn 
                dUdV=[1,0,0;u(1),rho(1),0;0.5*u(1)*u(1),rho(1)*u(1),1/(obj.prob.gamma-1)];
                droeF(:,1:3,1)=droeF(:,1:3,1)*dUdV;
                dUdV=[1,0,0;u(end),rho(end),0;0.5*u(end)*u(end),rho(end)*u(end),1/(obj.prob.gamma-1)];
                droeF(:,4:6,end)=droeF(:,4:6,end)*dUdV;
                
                %                 for kl= 2:obj.prob.nVol-1
%                     i=kl-1;
%                     J2L(3*(i-1)+1:3*i,1:6) = [-obj.prob.S(kl)*droeF(:,1:3,kl-1), -obj.prob.S(kl)*droeF(:,4:6,kl-1)];
%                     J2R(3*(i-1)+1:3*i,1:6) = [obj.prob.S(kl+1)*droeF(:,1:3,kl), obj.prob.S(kl+1)*droeF(:,4:6,kl)];
%                 end

%                 right = bsxfun(@times,roeF(:,2:end),obj.prob.S(3:end-1));
%                 left=  -bsxfun(@times,roeF(:,1:end-1),obj.prob.S(2:end-2));
                start=1;
                finish=points(2)-1;
                right=   roeF(:,finish)*obj.prob.S(finish+1);
                left = - roeF(:,start)*obj.prob.S(start+1);
%                 keyboard
                current_points = 2:points(2)-1;
                cnstr(1:3,:)=[(obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi1(current_points,:)*w_increment;
                    (obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi2(current_points,:)*w_increment;
                    (obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi3(current_points,:)*w_increment ]+...
                    obj.time.dt*(right+left-Q(:,current_points)*(obj.prob.SVol(current_points).*obj.prob.dx(current_points))');

                derivQ1=zeros(3,obj.trunc);
                for i=2:points(2)-1 %1:points(2)-1
                    derivQ1=derivQ1+squeeze(dQ(:,:,i))*[Phi1(i,:);Phi2(i,:);Phi3(i,:)]*(obj.prob.SVol(i).*obj.prob.dx(i));
                end

%                 %keyboard
%                 Dcnstr(1:3,:) = [(obj.prob.SVol(2:points(2)-1).*obj.prob.dx(2:points(2)-1))*Phi1(2:points(2)-1,:);
%                     (obj.prob.SVol(2:points(2)-1).*obj.prob.dx(2:points(2)-1))*Phi2(2:points(2)-1,:);
%                     (obj.prob.SVol(2:points(2)-1).*obj.prob.dx(2:points(2)-1))*Phi3(2:points(2)-1,:)]+...
%                     obj.time.dt*(J2L(1:3,1:3)*Phi(1:3,:) +J2L(1:3,4:6)*Phi(4:6,:)+ ...
%                      J2R(3*(points(2)-3)+1:3*(points(2)-2),1:3)*Phi(3*(points(2)-2)+1:3*(points(2)-1),:)+J2R(3*(points(2)-3)+1:3*(points(2)-2),4:6)*Phi(3*(points(2)-1)+1:3*points(2),:)- derivQ1);
%                 % J2R(3*(points(2)-3)+1:3*(points(2)-2),1:3)*Phi(3*(points(2)-2)+1:3*(points(2)-1),:)+J2R(3*(points(2)-3)+1:3*(points(2)-2),4:6)*Phi(3*(points(2)-1)+1:3*points(2),:)- derivQ1);

 derJL= -obj.prob.S(start+1)*droeF(:,:,start); % 3x6 , first dimension  is w.r.t. (rho, rho*u, e); second dimension is w.r.t. cells j and j+1
                derJR=  obj.prob.S(finish+1)*droeF(:,:,finish);

                Dcnstr(1:3,:)= [(obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi1(current_points,:);
                    (obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi2(current_points,:);
                    (obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi3(current_points,:)]+...
                    obj.time.dt*(derJL*Phi(3*(start-1)+1:3*(start+1),:)+ ...
                    derJR*Phi(3*(finish-1)+1:3*(finish+1),:)- derivQ1);


                if obj.ncell>2
                    for kk=2:obj.ncell-1
                        start=points(kk)-1;
                        finish=points(kk+1)-1;
                        right=   roeF(:,finish)*obj.prob.S(finish+1);
                        left = - roeF(:,start)*obj.prob.S(start+1);
                        current_points=start+1:finish;
                        cnstr(3*(kk-1)+1:3*kk,:)=[(obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi1(current_points,:)*w_increment;
                            (obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi2(current_points,:)*w_increment;
                            (obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi3(current_points,:)*w_increment]+...
                            obj.time.dt*(right+left-Q(:,current_points)*(obj.prob.SVol(current_points).*obj.prob.dx(current_points))');
                        %obj.time.dt*(right(:,points(kk+1)-2)+left(:,points(kk)-1)-Q(:,points(kk):points(kk+1)-1)*(obj.prob.SVol(points(kk):points(kk+1)-1).*obj.prob.dx(points(kk):points(kk+1)-1))');
                        derivQ=zeros(3,obj.trunc);
                        for ii=points(kk):points(kk+1)-1
                            derivQ=derivQ+squeeze(dQ(:,:,ii))*[Phi1(ii,:);Phi2(ii,:);Phi3(ii,:)]*(obj.prob.SVol(ii).*obj.prob.dx(ii));
                        end
                        %keyboard
                        %                         Dcnstr(3*(kk-1)+1:3*kk,:)=[(obj.prob.SVol(points(kk):points(kk+1)-1).*obj.prob.dx(points(kk):points(kk+1)-1))*Phi1(points(kk):points(kk+1)-1,:);
                        %                             (obj.prob.SVol(points(kk):points(kk+1)-1).*obj.prob.dx(points(kk):points(kk+1)-1))*Phi2(points(kk):points(kk+1)-1,:);
                        %                             (obj.prob.SVol(points(kk):points(kk+1)-1).*obj.prob.dx(points(kk):points(kk+1)-1))*Phi3(points(kk):points(kk+1)-1,:)]+...
                        %                             obj.time.dt*(J2L(3*(points(kk)-2)+1:3*(points(kk)-1),1:3)*Phi(3*(points(kk)-2)+1:3*(points(kk)-1),:) +J2L(3*(points(kk)-2)+1:3*(points(kk)-1),4:6)*Phi(3*(points(kk)-1)+1:3*(points(kk)),:)+ ...
                        %                             J2R(3*(points(kk+1)-3)+1:3*(points(kk+1)-2),1:3)*Phi(3*(points(kk+1)-2)+1:3*(points(kk+1)-1),:)+ J2R(3*(points(kk+1)-3)+1:3*(points(kk+1)-2),4:6)*Phi(3*(points(kk+1)-1)+1:3*(points(kk+1)),:)- derivQ);

                        derJL= -obj.prob.S(start+1)*droeF(:,:,start); % 3x6 , first dimension  is w.r.t. (rho, rho*u, e); second dimension is w.r.t. cells j and j+1 
                        derJR= obj.prob.S(finish+1)*droeF(:,:,finish);

                        Dcnstr(3*(kk-1)+1:3*kk,:)= [(obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi1(current_points,:);
                            (obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi2(current_points,:);
                            (obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi3(current_points,:)]+...
                            obj.time.dt*(derJL*Phi(3*(start-1)+1:3*(start+1),:)+ ...
                            derJR*Phi(3*(finish-1)+1:3*(finish+1),:)- derivQ);

                    end
                end


%                 keyboard
                start=points(end-1)-1;
                finish=points(end)-2;
                right=  roeF(:,finish)*obj.prob.S(finish+1);
                left = -roeF(:,start)*obj.prob.S(start+1);
                current_points=start+1:finish;
                cnstr((obj.ncell-1)*3+1:obj.ncell*3,:)=[(obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi1(current_points,:)*w_increment;
                    (obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi2(current_points,:)*w_increment;
                    (obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi3(current_points,:)*w_increment ]+...
                    obj.time.dt*(right+left-Q(:,current_points)*(obj.prob.SVol(current_points).*obj.prob.dx(current_points))');
%                     obj.time.dt*(right(:,end)+left(:,points(end-1)-1)-Q(:,points(end-1):end)*(obj.prob.SVol(points(end-1):end).*obj.prob.dx(points(end-1):end))');
                derivQ3=zeros(3,obj.trunc);
                for ii= points(end-1):obj.prob.nVol-1 %points(end-1):obj.prob.nVol
                    derivQ3=derivQ3+squeeze(dQ(:,:,ii))*[Phi1(ii,:);Phi2(ii,:);Phi3(ii,:)]*(obj.prob.SVol(ii).*obj.prob.dx(ii));
                end

                %keyboard
%                 Dcnstr((obj.ncell-1)*3+1:obj.ncell*3,:) = [(obj.prob.SVol(points(end-1):end-1).*obj.prob.dx(points(end-1):end-1))*Phi1(points(end-1):end-1,:);
%                     (obj.prob.SVol(points(end-1):end-1).*obj.prob.dx(points(end-1):end-1))*Phi2(points(end-1):end-1,:);
%                     (obj.prob.SVol(points(end-1):end-1).*obj.prob.dx(points(end-1):end-1))*Phi3(points(end-1):end-1,:)]+...
%                     obj.time.dt*(J2L(3*(points(end-1)-2)+1:3*(points(end-1)-1),1:3)*Phi(3*(points(end-1)-2)+1:3*(points(end-1)-1),:)+J2L(3*(points(end-1)-2)+1:3*(points(end-1)-1),4:6)*Phi(3*(points(end-1)-1)+1:3*(points(end-1)),:)+ ...
%                     J2R(end-2:end,1:3)*Phi(end-5:end-3,:)+ J2R(end-2:end,4:6)*Phi(end-2:end,:)- derivQ3);
                derJL=-obj.prob.S(start+1)* droeF(:,:,start); % 3x6 , first dimension  is w.r.t. (rho, rho*u, e); second dimension is w.r.t. cells j and j+1
%                 finish=size(roeF,2);
                derJR=obj.prob.S(finish+1)* droeF(:,:,finish);

                Dcnstr((obj.ncell-1)*3+1:obj.ncell*3,:)= [(obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi1(current_points,:);
                    (obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi2(current_points,:);
                    (obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi3(current_points,:)]+...
                     obj.time.dt*(derJL*Phi(3*(start-1)+1:3*(start+1),:)+ ...
                                  derJR*Phi(3*(finish-1)+1:3*(finish+1),:)- derivQ3);



    end

  
  
  
  %%%%
  function [c, dc]=testmyconstraints(obj, w_increment,i)
            [Cnst,Dcnstr]=obj.myconstraints(w_increment);
            c=Cnst(i); dc=Dcnstr(i,:)';
        end
       
        function  [] = solveConstraints(obj)
            %use if 3*obj.ncell==obj.trunc
            
            %keyboard
            itnump1 = obj.cTimeIter + 1;
            Phi=obj.phi(:,1:obj.trunc);
            obj.sv(:,itnump1)=obj.sv(:,obj.cTimeIter);
            
            w_guess = Phi'* (obj.sv(:,itnump1)-obj.sv(:,obj.cTimeIter));
            if obj.printLevel > 1
                fprintf('---- Newton Step    # --------- ||R|| ------------ ||du|| ------------- \n');
            end
            %keyboard
            if obj.cTimeIter==1 
                size(Phi)
                disp('num constraint=num basis vectors ')
                disp(['several domaines ', num2str(obj.ncell)])
            end
            
            for k=1:20
                
                [g,Dg]=obj.myconstraints(w_guess);
                
                for i=1:length(g)
                    out=gradientcheck( @(w) obj.testmyconstraints(w,i), w_guess);
                    if out.RelError>1e-4,keyboard, end
                end
                
                if norm(g,2)<10^(-6)
                    break
                end
                del_w=-Dg\g;
                w_guess=w_guess+del_w;
                
            end
            obj.sv(:,itnump1)=obj.sv(:,obj.cTimeIter)+Phi*w_guess;
            [g,~]=obj.myconstraints(w_guess);
            %if norm(g)>10^(-4), keyboard, end
            disp(['norm of the constraints ', num2str(norm(g))])
            disp(['number of iterations  ',num2str(k)])
        end
        
     
        function  [] = minimizeConstraints(obj)
            % %use if 3*obj.ncell > obj.trunc
            itnump1 = obj.cTimeIter + 1;
            Phi=obj.phi(:,1:obj.trunc);
            obj.sv(:,itnump1)=obj.sv(:,obj.cTimeIter);
            
            w_guess = Phi'* (obj.sv(:,itnump1)-obj.sv(:,obj.cTimeIter));
            if obj.printLevel > 1
                fprintf('---- Newton Step    # --------- ||R|| ------------ ||du|| ------------- \n');
            end
            %keyboard
            if obj.cTimeIter==1 
                size(Phi)
                disp('num constraint=num basis vectors ')
                disp(['several domaines ', num2str(obj.ncell)])
            end
            
            for k=1:50
                
                [g,Dg]=obj.myconstraints(w_guess);
                if norm(g,2)<10^(-6)
                    break
                end
                
                J=Dg'*g; Hess=Dg'*Dg;
                del_w=-Hess\J;
                w_guess=w_guess+del_w;
                if norm(del_w,2)<10^(-6)
                    break
                end
            end
            
            [g,~]=obj.myconstraints(w_guess);
            %if norm(g)>10^(-4), keyboard, end
            disp(['norm of the constraints ', num2str(norm(g))])
            disp(['number of iterations  ', num2str(k)])
            obj.sv(:,itnump1)=obj.sv(:,obj.cTimeIter)+Phi*w_guess;
            
        end
               
%%%%%%%%%%%%%%%%%%%%%%% Multiple domains for GNAT  %%%%%%%%%%%
function [cnstr, Dcnstr]=constrMultiDomain(obj,w_increment,SV,delt)
    
    dn=floor(obj.prob.nVol/(obj.ncell));
    points=1:dn:obj.prob.nVol+1;
    if points(end)~=obj.prob.nVol+1 && points(end)~=obj.prob.nVol
        points(end+1)=obj.prob.nVol+1;
    else
        points(end)= obj.prob.nVol+1;
    end
    
    Phi=obj.phi(:,1:obj.trunc);
    Phi1=Phi(1:3:end,:); %basis for the first conserved quantity rho
    Phi2=Phi(2:3:end,:); %basis for the second conserved quantity rho*u
    Phi3=Phi(3:3:end,:); %basis for the thirs conserved quantity e
    
    
    [rho, u, P,c,e,~,dc]=obj.prob.getVariables(SV);
    [Q,dQ]=obj.prob.forceTerm(u,P);
    [roeF, droeF]=obj.prob.roeFlux(rho,u,P,c,e,dc);
    
    %from governEqn
    dUdV=[1,0,0;u(1),rho(1),0;0.5*u(1)*u(1),rho(1)*u(1),1/(obj.prob.gamma-1)];
    droeF(:,1:3,1)=droeF(:,1:3,1)*dUdV;
    dUdV=[1,0,0;u(end),rho(end),0;0.5*u(end)*u(end),rho(end)*u(end),1/(obj.prob.gamma-1)];
    droeF(:,4:6,end)=droeF(:,4:6,end)*dUdV;
    
    start=1;
    finish=points(2)-1;
    right=   roeF(:,finish)*obj.prob.S(finish+1);
    left = - roeF(:,start)*obj.prob.S(start+1);
    
    current_points = 2:points(2)-1;
    cnstr(1:3,:)=[(obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi1(current_points,:)*w_increment;
        (obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi2(current_points,:)*w_increment;
        (obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi3(current_points,:)*w_increment ]+...
        delt*(right+left-Q(:,current_points)*(obj.prob.SVol(current_points).*obj.prob.dx(current_points))');
    
    derivQ1=zeros(3,obj.trunc);
    for i=2:points(2)-1 %1:points(2)-1
        derivQ1=derivQ1+squeeze(dQ(:,:,i))*[Phi1(i,:);Phi2(i,:);Phi3(i,:)]*(obj.prob.SVol(i).*obj.prob.dx(i));
    end
    
    derJL= -obj.prob.S(start+1)*droeF(:,:,start); % 3x6 , first dimension  is w.r.t. (rho, rho*u, e); second dimension is w.r.t. cells j and j+1
    derJR=  obj.prob.S(finish+1)*droeF(:,:,finish);
    
    Dcnstr(1:3,:)= [(obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi1(current_points,:);
        (obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi2(current_points,:);
        (obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi3(current_points,:)]+...
        delt*(derJL*Phi(3*(start-1)+1:3*(start+1),:)+ ...
        derJR*Phi(3*(finish-1)+1:3*(finish+1),:)- derivQ1);
    
    
    if obj.ncell>2
        for kk=2:obj.ncell-1
            start=points(kk)-1;
            finish=points(kk+1)-1;
            right=   roeF(:,finish)*obj.prob.S(finish+1);
            left = - roeF(:,start)*obj.prob.S(start+1);
            current_points=start+1:finish;
            cnstr(3*(kk-1)+1:3*kk,:)=[(obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi1(current_points,:)*w_increment;
                (obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi2(current_points,:)*w_increment;
                (obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi3(current_points,:)*w_increment]+...
                delt*(right+left-Q(:,current_points)*(obj.prob.SVol(current_points).*obj.prob.dx(current_points))');
            derivQ=zeros(3,obj.trunc);
            for ii=points(kk):points(kk+1)-1
                derivQ=derivQ+squeeze(dQ(:,:,ii))*[Phi1(ii,:);Phi2(ii,:);Phi3(ii,:)]*(obj.prob.SVol(ii).*obj.prob.dx(ii));
            end
            
            derJL= -obj.prob.S(start+1)*droeF(:,:,start); % 3x6 , first dimension  is w.r.t. (rho, rho*u, e); second dimension is w.r.t. cells j and j+1
            derJR= obj.prob.S(finish+1)*droeF(:,:,finish);
            
            Dcnstr(3*(kk-1)+1:3*kk,:)= [(obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi1(current_points,:);
                (obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi2(current_points,:);
                (obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi3(current_points,:)]+...
                delt*(derJL*Phi(3*(start-1)+1:3*(start+1),:)+ ...
                derJR*Phi(3*(finish-1)+1:3*(finish+1),:)- derivQ);
            
        end
    end
    
    
    start=points(end-1)-1;
    finish=points(end)-2;
    right=  roeF(:,finish)*obj.prob.S(finish+1);
    left = -roeF(:,start)*obj.prob.S(start+1);
    current_points=start+1:finish;
    cnstr((obj.ncell-1)*3+1:obj.ncell*3,:)=[(obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi1(current_points,:)*w_increment;
        (obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi2(current_points,:)*w_increment;
        (obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi3(current_points,:)*w_increment ]+...
        delt*(right+left-Q(:,current_points)*(obj.prob.SVol(current_points).*obj.prob.dx(current_points))');
    
    derivQ3=zeros(3,obj.trunc);
    for ii= points(end-1):obj.prob.nVol-1 %points(end-1):obj.prob.nVol
        derivQ3=derivQ3+squeeze(dQ(:,:,ii))*[Phi1(ii,:);Phi2(ii,:);Phi3(ii,:)]*(obj.prob.SVol(ii).*obj.prob.dx(ii));
    end
    
    derJL=-obj.prob.S(start+1) * droeF(:,:,start); % 3x6 , first dimension  is w.r.t. (rho, rho*u, e); second dimension is w.r.t. cells j and j+1
    derJR= obj.prob.S(finish+1)* droeF(:,:,finish);
    Dcnstr((obj.ncell-1)*3+1:obj.ncell*3,:)= [(obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi1(current_points,:);
        (obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi2(current_points,:);
        (obj.prob.SVol(current_points).*obj.prob.dx(current_points))*Phi3(current_points,:)]+...
        delt*(derJL*Phi(3*(start-1)+1:3*(start+1),:)+ ...
        derJR*Phi(3*(finish-1)+1:3*(finish+1),:)- derivQ3);
    
    
    
end

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function  [] = NewtonRaphsonLocal(obj)
            %This function solves systems of the form:  f(x) = 0 using
            %Newton-Raphson's method
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - ROM object
            %
            %Outputs:
            %--------
            %No outputs. The state vector at the current time
            %(indicated by the cTimeIter property) is stored in the FOM
            %handle class.
            %--------------------------------------------------------------
            
            %Determine time iteration number plus 1
            itnump1 = obj.cTimeIter + 1;
            
            %Extract current basis
            cLocBasis = obj.LocBasisHist(itnump1-1,1);
            
            %On first time step, extract appropriate basis from Master Basis
            if obj.cTimeIter == 1
                obj.phi = obj.phi0{cLocBasis};
                obj.UrLoc{1,1} = zeros(obj.nY(cLocBasis),1);
            end
            
            %%%%%%UPDATE BASIS ON THE FLY%%%%%%
            if strcmpi(obj.basisUpdate,'border_exact') && ((itnump1 > 2 &&  cLocBasis ~= obj.LocBasisHist(itnump1-2,1)) || (itnump1-1 == 1 && ~strcmpi(obj.snaps.initref,'ic')))
                %If it is AFTER the first timestep and the current basis
                %differs from the previous basis OR if it is the first time
                %step and the initial snapshot reference was NOT the
                %initial condition, perform exact SVD update
                obj.svdUpdateData(cLocBasis).snapMAT = bsxfun(@minus, obj.svdUpdateData(cLocBasis).snapMAT, obj.sv(:,itnump1-1) - obj.svdUpdateData(cLocBasis).prevRef);
                obj.svdUpdateData(cLocBasis).prevRef = obj.sv(:,itnump1-1);
                obj.phi = pod(obj.svdUpdateData(cLocBasis).snapMAT,obj.nY(cLocBasis),[]);
                obj.numswitch = obj.numswitch + 1;
                obj.UrLoc{obj.numswitch,1} = zeros(obj.nY(cLocBasis),1);
            elseif (strcmpi(obj.basisUpdate,'border') || strcmpi(obj.basisUpdate,'border_fast_approx')) && ((itnump1 > 2 && cLocBasis ~= obj.LocBasisHist(itnump1-2,1)) || ( itnump1-1 == 1 && ~strcmpi(obj.snaps.initref,'ic')))
                %If it is AFTER the first timestep and the current basis
                %differs from the previous basis OR if it is the first time
                %step and the initial snapshot reference was NOT the
                %initial condition, perform thin SVD update
                n = obj.nY(cLocBasis)+obj.svdUpdateData(cLocBasis).bufferSize;
                
                vec = obj.svdUpdateData(cLocBasis).wRef - obj.sv(:,itnump1-1);
                U = ThinSVD_AddVec2Cols_Fast([obj.phi0{cLocBasis,1},obj.svdUpdateData(cLocBasis).bufferedLSingVecs],...
                    obj.SingVals{cLocBasis}(1:n,1),obj.svdUpdateData(cLocBasis).ColSumRSVtrunc,...
                    obj.svdUpdateData(cLocBasis).normOneMinusVVt,vec);
                %obj.phi{cLocBasis,1} = U(:,1:obj.nY(cLocBasis));
                obj.phi = U(:,1:obj.nY(cLocBasis));
                obj.FullBasis=U;
                obj.numswitch = obj.numswitch + 1;
                obj.UrLoc{obj.numswitch,1} = zeros(obj.nY(cLocBasis),1);
            elseif strcmpi(obj.basisUpdate,'border_fast') && ((itnump1 > 2 && cLocBasis ~= obj.LocBasisHist(itnump1-2,1))||( itnump1-1 == 1 && ~strcmpi(obj.snaps.initref,'ic')))
                %If it is AFTER the first timestep and the current basis
                %differs from the previous basis OR if it is the first time
                %step and the initial snapshot reference was NOT the
                %initial condition, perform fast SVD update
                obj.efficientOnlineSVDupdate(cLocBasis,obj.sv(:,itnump1-1));
            elseif (strcmpi(obj.basisUpdate,'border_exact') || strcmpi(obj.basisUpdate,'border') || strcmpi(obj.basisUpdate,'border_fast') || strcmpi(obj.basisUpdate,'border_fast_approx')) && strcmpi(obj.snaps.initref,'ic') && itnump1 == 2
                %If we are using fast SVD updating and the initial
                %reference is the initial condition, advance numswitch
                %because we want to skip the first update.
                obj.numswitch = obj.numswitch + 1;
            elseif strcmpi(obj.basisUpdate,'none') && itnump1 > 2 && cLocBasis ~= obj.LocBasisHist(itnump1-2,1)
                %If it is AFTER the first timestep and the current basis
                %differs from the previous basis AND we are not doing
                %updating, simply load the current basis from the Master
                %Basis.
                obj.phi = obj.phi0{cLocBasis,1};
            end
            
            if strcmpi(obj.basisUpdate,'none')
                %Since numswitch is used to pull information out of Ur_Loc
                %(i.e. the q's from David's paper), we need to set it to be
                %the local basis in the non-updating case.
                obj.numswitch=cLocBasis;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %Determine the residual and jacobian based on initial guess
            t = obj.time.T(1) + obj.time.dt*obj.cTimeIter;
            [R,J] = obj.TimeScheme.TimeIntNLFunc(obj.prob,obj.sv(:,itnump1-1),obj.sv(:,itnump1-1),t);
            
            %Determine convergence criteria to use
            [tol,tolIt] = determineConvergeCriterion(obj,norm(R,2));
            
            indexAdj = 1;
            for i_N = 1:obj.newt.maxIter
                %Solve for the search direction and update state vector
                if lower(obj.Rtype(1)) == 'g'
                    newR = obj.phi'*R;
                    du = -((obj.phi'*J*obj.phi)\newR);
                    conv = norm(newR,2);
                else
                    newJ = J*obj.phi;
                    du = -(newJ\R);
                    conv = norm(newJ'*R,2);
                end
                
                %If nan or complex encountered -> kill simulation
                if sum(isnan(du)) || sum(real(du)~=du)
                    obj.killflag=true;
                    return;
                end
                
                p = obj.phi*du;
                %Linesearch if necessary
                if obj.newt.linesrch.status
                    alpha = linesearch(@(alpha) obj.TimeScheme.TimeIntNLFuncDir(obj.prob,obj.sv(:,itnump1-indexAdj),obj.sv(:,itnump1-1),t,alpha,p),obj.newt.linesrch.prop);
                    p = alpha*p;
                end
                
                obj.sv(:,itnump1) = obj.sv(:,itnump1-indexAdj) + p;
                obj.UrLoc{obj.numswitch,1} = obj.UrLoc{obj.numswitch,1} + du;
                
                %Save residual and jacobian
                if i_N > 1
                    obj.writeNonlinearSnapshot(R,J,p,cLocBasis);
                end
                
                %Compute residual and jacobian with updated vector for next iteration
                [R,J] = obj.TimeScheme.TimeIntNLFunc(obj.prob,obj.sv(:,itnump1),obj.sv(:,itnump1-1),t);
                %                 obj.writeNonlinearSnapshot(R,J,p,cLocBasis);
                
                indexAdj = 0;
                %Stop iterating once convergence criteria is met
                if checkConverge(obj,conv,p,tol,tolIt)
                    break;
                end
            end
            %Write the last nonlinear snapshot.  We note that we are using
            %J^(k+1)*0 if using Snapshot 1.5 or 2 (we don't have a
            %search direction p^(k+1) because we converged!).  It shouldn't
            %matter much because ||p|| should be small since we converged,
            %i.e. the next step shouldn't take us away from the solution!
            obj.writeNonlinearSnapshot(R,J,zeros(size(p)),cLocBasis);
            
            %Store the number of newton iterations at this time step
            obj.newt.iter(obj.cTimeIter) = i_N;
            
            %If the maxiumum newton iterations are reached, warn the user
            if ~obj.newt.quiet && (i_N==obj.newt.maxIter)
                disp(['*** Warning; Newton solver reached max number of iterations before convergence : res = ',...
                    num2str(norm(R,2)/tol),'  eps = ',num2str(obj.newt.eps), '***']);
            end
        end %Done
        
        function  [] = NewtonRaphson_PrecompGal(obj)
            %This function solves systems of the form:  f(x) = 0 using
            %Newton-Raphson's method
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - ROM object
            %
            %Outputs:
            %--------
            %No outputs. The state vector at the current time
            %(indicated by the cTimeIter property) is stored in the FOM
            %handle class.
            %--------------------------------------------------------------
            
            %Determine time iteration number plus 1
            itnump1 = obj.cTimeIter + 1;
            
            %Extract current basis
            if obj.nBases > 1
                cLocBasis = obj.LocBasisHist(itnump1-1,1);
            else
                cLocBasis = 1;
                obj.numswitch=1;
            end
            
            %%%%%%Basis Updating NOT SUPPORTED here (yet)%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %Determine the residual and jacobian based on initial guess
            t = obj.time.T(1) + obj.time.dt*obj.cTimeIter;
            [R,J] = obj.TimeScheme.TimeIntNLFunc(obj.prob,obj.sv(1:obj.nY(cLocBasis),itnump1-1),obj.sv(1:obj.nY(cLocBasis),itnump1-1),t,true);
            %[R,J] = obj.prob.ResJacPrecomp(obj.sv(1:obj.nY(cLocBasis),itnump1-1),t,cLocBasis);
            
            %Determine convergence criteria to use
            [tol,tolIt] = determineConvergeCriterion(obj,norm(R,2));
            
            indexAdj = 1;
            for i_N = 1:obj.newt.maxIter
                dy = -J\R;
                
                obj.sv(1:obj.nY(cLocBasis),itnump1) = obj.sv(1:obj.nY(cLocBasis),itnump1-indexAdj) + dy;
                obj.UrLoc{obj.numswitch,1} = obj.UrLoc{cLocBasis,1} + dy;
                
                %Compute residual and jacobian with updated vector for next iteration
                %[R,J] = obj.prob.ResJacPrecomp(obj.sv(1:obj.nY(cLocBasis),itnump1),t,cLocBasis);
                [R,J] = obj.TimeScheme.TimeIntNLFunc(obj.prob,obj.sv(1:obj.nY(cLocBasis),itnump1),obj.sv(1:obj.nY(cLocBasis),itnump1-1),t,true);
                
                indexAdj = 0;
                %Stop iterating once convergence criteria is met
                if checkConverge(obj,norm(R),dy,tol,tolIt)
                    break;
                end
            end
            %Store the number of newton iterations at this time step
            obj.newt.iter(obj.cTimeIter) = i_N;
            
            %If the maxiumum newton iterations are reached, warn the user
            if ~obj.newt.quiet && (i_N==obj.newt.maxIter)
                disp(['*** Warning; Newton solver reached max number of iterations before convergence : res = ',...
                    num2str(norm(R,2)/tol),'  eps = ',num2str(obj.newt.eps), '***']);
            end
        end
        
        function  [] = precomputeLocalSimulation(obj)
            %This function precomputes quantities used in the local ROB
            %simulations.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %trainIC - ndof x 1 vector containing the training initial
            %          condition
            %
            %Outputs:
            %--------
            %This function has no outputs.  It stores the quantites f and S
            %in the ROM object, which are precomputed for the local
            %simulation.
            %--------------------------------------------------------------
            
            if obj.nBases == 1, return; end
            
            if strcmpi(obj.basisUpdate,'none')
                for i1 = 1:obj.nBases
                    for i2 = 1:obj.nBases
                        obj.preCompDistQuant(i1,i2).d = norm(obj.clusCenter(:,i1)-obj.prob.ic,2)^2 - norm(obj.clusCenter(:,i2)-obj.prob.ic,2)^2;
                        for j = 1:obj.nBases
                            obj.preCompDistQuant(i1,i2).g{j} = 2*(obj.clusCenter(:,i2)-obj.clusCenter(:,i1))'*obj.phi0{j};
                        end
                    end
                end
            elseif strcmpi(obj.basisUpdate,'border_fast') || strcmpi(obj.basisUpdate,'border_fast_approx')
                for i1 = 1:obj.nBases
                    for i2 = 1:obj.nBases
                        obj.preCompDistQuant(i1,i2).d = norm(obj.clusCenter(:,i1)-obj.prob.ic,2)^2 - norm(obj.clusCenter(:,i2)-obj.prob.ic,2)^2;
                        obj.preCompDistQuant(i1,i2).e = 2*(obj.clusCenter(:,i2)-obj.clusCenter(:,i1))'*obj.prob.ic;
                        for j = 1:obj.nBases
                            obj.preCompDistQuant(i1,i2).f{j} = 2*(obj.clusCenter(:,i2)-obj.clusCenter(:,i1))'*obj.svdUpdateData(j).wRef;
                            obj.preCompDistQuant(i1,i2).g{j} = 2*(obj.clusCenter(:,i2)-obj.clusCenter(:,i1))'*obj.phi0{j};
                        end
                    end
                end
            end
        end
        
        function  [S_kp1] = GN_HessianMod(obj,S_k,s,Jr_kp1,r_kp1,Jr_k,r_k)
            
            
            y = Jr_kp1'*r_kp1 - Jr_k'*r_k;
            y_sh = Jr_kp1'*r_kp1 - Jr_k'*r_kp1;
            
            tau = min(1,abs(s'*y_sh)/abs(s'*S_k*s));
            S_k = tau*S_k;
            
            rho = y'*s;
            y = y/rho;
            v = (y_sh - S_k*s);
            
            S_kp1 = S_k + (v*y' + y*v') - ((v'*s)*y)*y';
            
            
        end
        
        %Open/Close/Read/Write Files
        function  [] = openCloseResJacFiles(obj,flag)
            %This function opens/closes the binary files that will be used
            %to store the residual and jacobian snapshots (if necessary)
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - ROM object
            %flag    - string indicating whether to open or close the file
            %
            %Outputs:
            %--------
            %There are no outputs.
            %--------------------------------------------------------------
            
            %If nonlinear terms are not saved, exit.
            if obj.saveNL == 0
                return;
            end
            
            if strcmpi(flag,'open')
                if obj.nBases == 1
                    %Open res/jac files for the single basis case
                    if obj.saveNL == 1 || obj.saveNL == 1.5
                        obj.res = fopen([obj.fileRoot,'ResSnaps_basis1.bin'],'wb');
                    else
                        obj.res = fopen([obj.fileRoot,'ResSnaps_basis1.bin'],'wb');
                        obj.jac = fopen([obj.fileRoot,'JacSnaps_basis1.bin'],'wb');
                    end
                else
                    %Open res/jac files for the multiple basis case
                    for i = 1:obj.nBases
                        if obj.saveNL == 1 || obj.saveNL == 1.5
                            obj.res{i} = fopen([obj.fileRoot,'ResSnaps_basis',num2str(i),'.bin'],'wb');
                        elseif obj.saveNL == 2 || obj.saveNL == 3
                            obj.res{i} = fopen([obj.fileRoot,'ResSnaps_basis',num2str(i),'.bin'],'wb');
                            obj.jac{i} = fopen([obj.fileRoot,'JacSnaps_basis',num2str(i),'.bin'],'wb');
                        end
                    end
                end
            elseif strcmpi(flag,'close')
                if obj.nBases == 1
                    %Close res/jac files for the single basis case
                    if obj.saveNL == 1 || obj.saveNL == 1.5
                        fclose(obj.res);
                    else
                        fclose(obj.res);
                        fclose(obj.jac);
                    end
                else
                    %Close res/jac files for the multiple basis case
                    for i = 1:obj.nBases
                        if obj.saveNL == 1 || obj.saveNL == 1.5
                            fclose(obj.res{i});
                        elseif obj.saveNL == 2 || obj.saveNL == 3
                            fclose(obj.res{i});
                            fclose(obj.jac{i});
                        end
                    end
                end
            end
        end %Done
        
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
            %res      - matrix containing the residual snapshots (if saveNL
            %           is 1, 2, or 3). Otherwise, res = [].
            %jac      - matrix containing the jacobian snapshots (if saveNL
            %           is 2 or 3).  Otherwise, jac = [].
            %--------------------------------------------------------------
            
            if obj.nBases == 1
                basisNum = 1;
            end
            ny = obj.nY(basisNum);
            
            %Return if saveNL = 0
            if obj.saveNL == 0
                res = [];
                jac = [];
                disp('Nonlinear terms were not saved because saveNL = 0. Exiting.');
                return;
            end
            
            %Read residual/jacobian snapshots for each NLSnapColl case
            if NLSnapColl == 1 || NLSnapColl == 1.5
                fidR = fopen([obj.fileRoot,'ResSnaps_basis',num2str(basisNum),'.bin'],'r');
                
                res = reshape(fread(fidR,obj.ndof*obj.numres(basisNum),'double'),obj.ndof,obj.numres(basisNum));
                jac = [];
                
                fclose(fidR);
            elseif NLSnapColl == 2
                fidR = fopen([obj.fileRoot,'ResSnaps_basis',num2str(basisNum),'.bin'],'r');
                fidJ = fopen([obj.fileRoot,'JacSnaps_basis',num2str(basisNum),'.bin'],'r');
                
                res = reshape(fread(fidR,obj.ndof*obj.numres(basisNum),'double'),obj.ndof,obj.numres(basisNum));
                jac = reshape(fread(fidJ,obj.ndof*obj.numres(basisNum),'double'),obj.ndof,obj.numres(basisNum));
                
                fclose(fidR);
                fclose(fidJ);
            elseif NLSnapColl == 3
                fidR = fopen([obj.fileRoot,'ResSnaps_basis',num2str(basisNum),'.bin'],'r');
                fidJ = fopen([obj.fileRoot,'JacSnaps_basis',num2str(basisNum),'.bin'],'r');
                
                res = reshape(fread(fidR,obj.ndof*obj.numres(basisNum),'double'),obj.ndof,obj.numres(basisNum));
                jac = reshape(fread(fidJ,obj.ndof*ny*obj.numres(basisNum),'double'),obj.ndof,ny*obj.numres(basisNum));
                
                fclose(fidR);
                fclose(fidJ);
            end
        end %Done
        
        function  [] = writeNonlinearSnapshot(obj,R,J,p,cbase)
            %This function writes the nonlinear snapshots to the appropriate
            %file.
            
            %Handle the case with only 1 basis separately
            if obj.saveNL == 0
                %Note this part of the if statement was inserted to speed up the computation.
                %If we are not saving R, J then this line will be executed but the elseifs won't.
            elseif obj.saveNL == 1
                %Snapshot collection procedure 1 from Carlberg paper
                if obj.nBases == 1
                    fwrite(obj.res,R,'double');
                else
                    fwrite(obj.res{cbase},R,'double');
                end
                obj.numres(cbase)=obj.numres(cbase)+1;
            elseif obj.saveNL == 1.5
                if obj.nBases == 1
                    fwrite(obj.res,R,'double');
                    fwrite(obj.res,J*p,'double');
                else
                    fwrite(obj.res{cbase},R,'double');
                    fwrite(obj.res{cbase},J*p,'double');
                end
                obj.numres(cbase)=obj.numres(cbase)+2;
            elseif obj.saveNL == 2
                %Snapshot collection procedure 2 from Carlberg paper
                if obj.nBases == 1
                    fwrite(obj.res,R,'double');
                    fwrite(obj.jac,J*p,'double');
                else
                    fwrite(obj.res{cbase},R,'double');
                    fwrite(obj.jac{cbase},J*p,'double');
                end
                obj.numres(cbase)=obj.numres(cbase)+1;
            elseif obj.saveNL == 3
                %Snapshot collection procedure 3 from Carlberg paper
                if obj.nBases == 1
                    fwrite(obj.res,R,'double');
                    fwrite(obj.jac,J*obj.phi,'double');
                else
                    fwrite(obj.res{cbase},R,'double');
                    fwrite(obj.jac{cbase},J*obj.phi,'double');
                end
                obj.numres(cbase)=obj.numres(cbase)+1;
            end
            
        end %Done
        
        function  [X] = readClusSnapFiles(obj,basisNum,rootFName)
            %This function reads clustered snapshots from files.
            
            %If no rootFName specific, use default rootname
            if nargin < 3 || isempty(rootFName)
                rootFName = obj.fileRoot;
            end
            
            %Open file, read first 2 entries as snapshot matrix size, and
            %read the remainder of the file as snapshots
            fid = fopen([rootFName,'clusteredSnaps',num2str(basisNum),'.bin'],'r');
            snapSize = fread(fid,2,'double');
            X = reshape(fread(fid,prod(snapSize),'double'),snapSize(1),snapSize(2));
            fclose(fid);
        end %Done
        
        function  [U,Sdiag,V] = readSVDFiles(obj,basisNum,rootFName)
            %This function reads SVD from files.
            V=[];
            
            %If rootFName not specified, use default rootname
            if nargin < 3 || isempty(rootFName)
                rootFName = obj.fileRoot;
            end
            
            %Read size of U and U itself from the binary file
            fid = fopen([rootFName,'svdFactU',num2str(basisNum),'.bin'],'r');
            sizeUsvd = fread(fid,2,'double');
            U = reshape(fread(fid,prod(sizeUsvd),'double'),sizeUsvd(1),sizeUsvd(2));
            fclose(fid);
            
            %Read size of Sigma and Sigma itself from the binary file
            fid = fopen([rootFName,'svdFactS',num2str(basisNum),'.bin'],'r');
            sizeSsvd = fread(fid,2,'double');
            Sdiag = reshape(fread(fid,prod(sizeSsvd),'double'),sizeSsvd(1),sizeSsvd(2));
            fclose(fid);
            
            %Only if nBases > 1, read the size of V and V itself from the binary file
            if obj.nBases > 1
                fid = fopen([rootFName,'svdFactV',num2str(basisNum),'.bin'],'r');
                sizeVsvd = fread(fid,2,'double');
                V = reshape(fread(fid,prod(sizeVsvd),'double'),sizeVsvd(1),sizeVsvd(2));
                fclose(fid);
                
                %CLUSTER CENTERS!!!!
                
                %                 X = obj.readClusSnapFiles(basisNum);
                %                 obj.clusCenter(:,basisNum) = mean(X,2);
            end
        end %Done
        
        %Basis Updating
        function  [] = initializeROBcompon(obj,cbase)
            
            obj.ROBcompon(1).alpha = zeros(obj.nY(cbase),1);
            for i = 1:obj.nBases
                obj.ROBcompon(1).beta{i}  = zeros(obj.nY(cbase),1);
                if i == cbase
                    obj.ROBcompon(1).N{i} = eye(obj.nY(i));
                else
                    obj.ROBcompon(1).N{i} = zeros(obj.nY(i),obj.nY(cbase));
                end
            end
            
        end
        
        function  [] = efficientOnlineSVDupdate(obj,newBasisNum,wf_switch)
            obj.numswitch=obj.numswitch+1;
            
            obj.UrLoc{obj.numswitch,1} = zeros(obj.nY(newBasisNum),1);
            switch obj.numswitch
                case 1
                case 2
                    %obj.wr_switch{1} = obj.phi{prevBasisNum,1}'*(wf_switch - obj.prob.ic);
                    obj.wr_switch{1} = obj.phi'*(wf_switch - obj.prob.ic);
                otherwise
                    %obj.wr_switch{obj.numswitch-1} = obj.phi{prevBasisNum,1}'*(wf_switch-obj.prevSwitchState);
                    obj.wr_switch{obj.numswitch-1} = obj.phi'*(wf_switch-obj.prevSwitchState);
            end
            obj.prevSwitchState = wf_switch;
            
            jsubl = newBasisNum;
            
            alpha_switch = 1;
            for p = 1:obj.numswitch-1
                alpha_switch = alpha_switch + obj.ROBcompon(p).alpha'*obj.wr_switch{p};
            end
            
            beta_switch = zeros(obj.nBases,1);
            n_switch  = cell(obj.nBases,1);
            for i = 1:obj.nBases
                n_switch{i}=zeros(obj.nY(i),1);
                for p = 1:obj.numswitch-1
                    beta_switch(i) = beta_switch(i) + obj.ROBcompon(p).beta{i}'*obj.wr_switch{p};
                    n_switch{i} = n_switch{i} + obj.ROBcompon(p).N{i}*obj.wr_switch{p};
                end
            end
            
            alpha_a = -alpha_switch;
            
            beta_a  = zeros(obj.nBases,1);
            n_a     = cell(obj.nBases,1);
            for i = 1:obj.nBases
                beta_a(i) = -beta_switch(i) + (i == jsubl);
                n_a{i} = -n_switch{i};
            end
            
            m = alpha_a*obj.svdUpdateData(jsubl).d;
            for i = 1:obj.nBases
                m = m + beta_a(i)*obj.svdUpdateData(jsubl).e{i} + obj.svdUpdateData(jsubl).F{i}*n_a{i};
            end
            
            r = alpha_a^2*obj.svdUpdateData(jsubl).a^2;
            for i = 1:obj.nBases
                r = r + 2*alpha_a*beta_a(i)*obj.svdUpdateData(i).c;
                if i == jsubl
                    r = r + 2*alpha_a*obj.svdUpdateData(i).d'*(n_a{i} - m);
                    for p = 1:obj.nBases
                        r = r + obj.svdUpdateData(i).g(p)^2*beta_a(i)*beta_a(p);
                        
                        if p == jsubl
                            r = r + (n_a{jsubl}-m)'*obj.svdUpdateData(jsubl).F{jsubl}*(n_a{jsubl}-m) ...
                                + 2*beta_a(i)*obj.svdUpdateData(p).e{i}'*(n_a{p} - m);
                        else
                            r = r + (n_a{i}-m)'*obj.svdUpdateData(i).F{p}*n_a{p} ...
                                + 2*beta_a(i)*obj.svdUpdateData(p).e{i}'*n_a{p};
                        end
                    end
                else
                    r = r + 2*alpha_a*obj.svdUpdateData(i).d'*n_a{i};
                    for p = 1:obj.nBases
                        r = r + obj.svdUpdateData(i).g(p)^2*beta_a(i)*beta_a(p) ...
                            + 2*beta_a(i)*obj.svdUpdateData(p).e{i}'*n_a{p};
                        if p == jsubl
                            r = r + n_a{i}'*obj.svdUpdateData(i).F{p}*(n_a{p}-m) ...
                                + 2*beta_a(i)*obj.svdUpdateData(p).e{i}'*(n_a{p} - m);
                        else
                            r = r + n_a{i}'*obj.svdUpdateData(i).F{p}*n_a{p} ...
                                + 2*beta_a(i)*obj.svdUpdateData(p).e{i}'*n_a{p};
                        end
                    end
                end
            end
            
            if r<0% && abs(r) < 1e-12
                fprintf('Warning...Negative r.  Skipping Update.');
                r = abs(sqrt(r));
                
                
                obj.ROBcompon(obj.numswitch).alpha = zeros(obj.nY(newBasisNum),1);
                
                obj.ROBcompon(obj.numswitch).beta  = cell(1,obj.nBases);
                obj.ROBcompon(obj.numswitch).N     = cell(1,obj.nBases);
                for i = 1:obj.nBases
                    obj.ROBcompon(obj.numswitch).beta{i} = zeros(obj.nY(newBasisNum),1);
                    if i == newBasisNum
                        obj.ROBcompon(obj.numswitch).N{i} = eye(obj.nY(newBasisNum));
                    else
                        obj.ROBcompon(obj.numswitch).N{i} = zeros(obj.nY(i),obj.nY(newBasisNum));
                    end
                end
                obj.phi = obj.phi0{newBasisNum};
                return;
            else
                r=sqrt(r);
            end
            
            %             ae = obj.svdUpdateData(jsubl).wRef-wf_switch;
            %             me = obj.phi0{jsubl}'*ae;
            %             pe = ae - obj.phi0{jsubl}*me;
            %             rcorrect = norm(pe);
            
            %%%%%%
            
            alpha_p = alpha_a/r;
            beta_p  = zeros(obj.nBases,1);
            n_p     = cell(obj.nBases,1);
            for i = 1:obj.nBases
                beta_p(i) = beta_a(i)/r;
                if i == jsubl
                    n_p{i} = (n_a{i} - m)/r;
                else
                    n_p{i} = n_a{i}/r;
                end
            end
            
            %             a=obj.prob.ic*alpha_a;
            %             for i = 1:obj.nBases
            %                 a = a + obj.svdUpdateData(i).wRef*beta_a(i) + obj.phi0{i}*n_a{i};
            %             end
            %             ae = obj.svdUpdateData(jsubl).wRef-wf_switch;
            %
            %             p=obj.prob.ic*alpha_p;
            %             for i = 1:obj.nBases
            %                 p = p + obj.svdUpdateData(i).wRef*beta_p(i) + obj.phi0{i}*n_p{i};
            %             end
            %             me = obj.phi0{jsubl}'*ae;
            %             pe = ae - obj.phi0{jsubl}*me;% pe=pe/norm(pe,2);
            
            K = [diag(obj.SingVals{jsubl}(1:obj.nY(jsubl))),zeros(obj.nY(jsubl),1);zeros(1,obj.nY(jsubl)),0] + ...
                [m;r]*[obj.svdUpdateData(jsubl).ColSumRSVtrunc',obj.svdUpdateData(jsubl).normOneMinusVVt];
            [C,~,~] = svd(K,0);
            C11 = C(1:end-1,1:end-1);
            c21 = C(end,1:end-1)';
            
            obj.ROBcompon(obj.numswitch).alpha = alpha_p*c21;
            
            obj.ROBcompon(obj.numswitch).beta  = cell(obj.nBases,1);
            obj.ROBcompon(obj.numswitch).N     = cell(obj.nBases,1);
            for i = 1:obj.nBases
                obj.ROBcompon(obj.numswitch).beta{i} = beta_p(i)*c21;
                if i == jsubl
                    obj.ROBcompon(obj.numswitch).N{i} = C11 + n_p{i}*c21';
                else
                    obj.ROBcompon(obj.numswitch).N{i} =  n_p{i}*c21';
                end
            end
            
            obj.phi = obj.prob.ic*obj.ROBcompon(obj.numswitch).alpha';
            for i = 1:obj.nBases
                obj.phi = obj.phi + ...
                    obj.svdUpdateData(i).wRef*obj.ROBcompon(obj.numswitch).beta{i}' + ...
                    obj.phi0{i}*obj.ROBcompon(obj.numswitch).N{i};
            end
        end
        
        %External Functions
        function  [precompROM] = createCopy4Precomp(obj,precompflag,philam,phimu,phirho)
            precompROM = ROM();
            
            props = {'fileRoot','Rtype','id','cfgid','snaps','time','newt','nY',...
                'nYflag','saveNL','nBases','TimeScheme','augment','basisUpdate',...
                'numswitch'};
            for i = 1:length(props)
                precompROM.(props{i})=obj.(props{i});
            end
            precompROM.precompFlag = precompflag;
            
            precompROM.prob = obj.prob.createReducedCopy(obj,precompflag,philam,phimu,phirho);
            
            precompROM.basisUpdate = 'none';
            
            precompROM.svdUpdateData=[];
            precompROM.ROBcompon = [];
            precompROM.ClusterIndex = [];
            precompROM.SingVals = [];
            precompROM.res = [];
            precompROM.jac = [];
            precompROM.UrLoc = cell(precompROM.nBases,1);
            for i = 1:precompROM.nBases
                precompROM.UrLoc{i} = zeros(precompROM.nY(i),1);
            end
            
            precompROM.phi0 = [];
            precompROM.phi  = [];
            precompROM.prevSwitchState = [];
            precompROM.clusCenter = [];
            
            precompROM.sv = zeros(max(precompROM.nY),precompROM.time.nstep+1);
        end
        
        function  [] = computePOD(obj,flag,svdOptions,rootFname)
            %Build nLocal local POD bases.
            
            %Set defaults
            if nargin < 3
                svdOptions=[];
            end
            if nargin < 4 || isempty(rootFname)
                rootFname = obj.fileRoot;
            end
            
            if obj.nBases == 1
                switch lower(flag)
                    case 'fromclustsnaps'
                        %Read clustered snapshots, compute POD (using SVD
                        %options specified; default is exact SVD)
                        X = obj.readClusSnapFiles(1);
                        [LSingVecs,SVs,RSingVecs]  = pod(X,[],svdOptions);
                        obj.SingVals = diag(SVs); clear SVs;
                        
                        %Write SVD factors to files (not V because we don't
                        %need it in the global case)
                        fid = fopen([rootFname,'svdFactU1.bin'],'wb');
                        fwrite(fid,size(LSingVecs)','double');
                        fwrite(fid,LSingVecs,'double');
                        fclose(fid);
                        
                        fid = fopen([rootFname,'svdFactS1.bin'],'wb');
                        fwrite(fid,size(obj.SingVals)','double');
                        fwrite(fid,obj.SingVals,'double');
                        fclose(fid);
                    case 'fromsvd'
                        %If we want to compute POD directly from an SVD we
                        %already have available, read the SVD files
                        [LSingVecs,obj.SingVals,RSingVecs]  = obj.readSVDFiles(1,rootFname);
                end
                
                if obj.nYflag == 1 %(nYflag = 1 means we use specified nY)
                    %If nY specified, just use standard POD and
                    %return
                    newNY = min(obj.nY,size(LSingVecs,2));
                    if newNY ~= obj.nY
                        warning(['nY changed from ',num2str(obj.nY),' to ',num2str(newNY),' due to available number of snapshots!']);
                        obj.nY = newNY;
                    end
                    obj.phi = LSingVecs(:,1:obj.nY);
                else
                    
                    %Once SVD is computed or read, find relative energy and set
                    %nY (if we have made it to this point, nY was not already
                    %specified).
                    relenergy = cumsum(obj.SingVals)/sum(obj.SingVals);
                    obj.nY = find(relenergy > obj.nYrelEnergy,1);
                    
                    %Make sure nY is in the allowable range, if one is
                    %specified
                    if ~isempty(obj.nYMin)
                        if obj.nY < obj.nYMin
                            obj.nY = obj.nYMin;
                        end
                    end
                    if ~isempty(obj.nYMax)
                        if obj.nY > obj.nYMax
                            obj.nY = obj.nYMax;
                        end
                    end
                    
                    %Extract the ROB
                    obj.phi = LSingVecs(:,1:obj.nY);
                end
                
                %Store necessary data for using border
                %Rank of SVD update (including buffer)
                if isempty(obj.svdUpdateData), obj.svdUpdateData.bufferSize=0; end;
                n = obj.nY+obj.svdUpdateData.bufferSize;
                
                %Store the vector to initially reference the
                %snapshots
                if strcmpi(obj.snaps.initref,'ic')
                    obj.svdUpdateData.wRef          = obj.prob.ic;
                elseif strcmpi(obj.snaps.initref,'zero')
                    obj.svdUpdateData.wRef          = zeros(size(obj.prob.ic));
                else
                    %error('When using SVD updating, initref can be ''ic'', or ''zero''');
                end
                %Compute values required for online svd updating
                if ~isempty(RSingVecs)
                    obj.svdUpdateData.ColSumRSVtrunc    = sum(RSingVecs(:,1:n),1)';
                    obj.svdUpdateData.normOneMinusVVt   = norm(ones(size(RSingVecs,1),1) - RSingVecs(:,1:n)*obj.svdUpdateData.ColSumRSVtrunc,2);
                end
                obj.svdUpdateData.bufferedLSingVecs = LSingVecs(:,obj.nY+1:n);
            else
                for iLocBasis = 1:obj.nBases
                    switch lower(flag)
                        case 'fromclustsnaps'
                            %Read clustered snapshots and compute POD
                            X = obj.readClusSnapFiles(iLocBasis);
                            [LSingVecs,SVs,RSingVecs] = pod(X,[],svdOptions);
                            obj.SingVals{iLocBasis} = diag(SVs); clear SVs;
                            
                            %Write SVD factors to files
                            fid = fopen([rootFname,'svdFactU',num2str(iLocBasis),'.bin'],'wb');
                            fwrite(fid,size(LSingVecs)','double');
                            fwrite(fid,LSingVecs,'double');
                            fclose(fid);
                            
                            fid = fopen([rootFname,'svdFactS',num2str(iLocBasis),'.bin'],'wb');
                            fwrite(fid,size(obj.SingVals{iLocBasis})','double');
                            fwrite(fid,obj.SingVals{iLocBasis},'double');
                            fclose(fid);
                            
                            fid = fopen([rootFname,'svdFactV',num2str(iLocBasis),'.bin'],'wb');
                            fwrite(fid,size(RSingVecs)','double');
                            fwrite(fid,RSingVecs,'double');
                            fclose(fid);
                        case 'fromsvd'
                            %Read SVD factors from file
                            [LSingVecs,obj.SingVals{iLocBasis},RSingVecs] = obj.readSVDFiles(iLocBasis,rootFname);
                    end
                    
                    %Determine nY(iLocBasis) and the ROB
                    if obj.nYflag == 2 %(nYflag = 2 means we use energy to determine nY)
                        relenergy = cumsum(obj.SingVals{iLocBasis})/sum(obj.SingVals{iLocBasis});
                        obj.nY(iLocBasis,1) = find(relenergy > obj.nYrelEnergy,1);
                        obj.nY(iLocBasis,1) = min(obj.nY(iLocBasis),size(LSingVecs,2));
                        
                        %Make sure nY is in the allowable range, if one is
                        %specified
                        if ~isempty(obj.nYMin)
                            if obj.nY(iLocBasis) < obj.nYMin(iLocBasis)
                                obj.nY(iLocBasis,1) = obj.nYMin(iLocBasis);
                            end
                        end
                        if ~isempty(obj.nYMax)
                            if obj.nY(iLocBasis) > obj.nYMax(iLocBasis)
                                obj.nY(iLocBasis,1) = obj.nYMax(iLocBasis);
                            end
                        end
                    end
                    obj.nY(iLocBasis) = min(obj.nY(iLocBasis),size(LSingVecs,2));
                    %Master Basis (i.e. the original local bases)
                    obj.phi0{iLocBasis,1} = LSingVecs(:,1:obj.nY(iLocBasis));
                    
                    %Store snapshots if using border_exact
                    if strcmpi(obj.basisUpdate,'border_exact')
                        obj.svdUpdateData(iLocBasis).snapMAT = X;
                        if strcmpi(obj.snaps.initref,'ic')
                            obj.svdUpdateData(iLocBasis).prevRef   = obj.prob.ic;
                        elseif strcmpi(obj.snaps.initref,'clustcenters')
                            obj.svdUpdateData(iLocBasis).prevRef   = obj.clusCenter(:,iLocBasis);
                        elseif strcmpi(obj.snaps.initref,'zero')
                            obj.svdUpdateData(iLocBasis).prevRef   = zeros(size(obj.prob.ic));
                        else
                            error('When using SVD updating, initref can be ''ic'', ''clustcenters'', or ''zero''');
                        end
                    end
                    
                    %Store necessary data for using border
                    if strcmpi(obj.basisUpdate,'border')
                        %Rank of SVD update (including buffer)
                        n = obj.nY(iLocBasis)+obj.svdUpdateData(iLocBasis).bufferSize;
                        
                        %Store the vector to initially reference the
                        %snapshots
                        if strcmpi(obj.snaps.initref,'ic')
                            obj.svdUpdateData(iLocBasis).wRef          = obj.prob.ic;
                        elseif strcmpi(obj.snaps.initref,'clustcenters')
                            obj.svdUpdateData(iLocBasis).wRef          = obj.clusCenter(:,iLocBasis);
                        elseif strcmpi(obj.snaps.initref,'zero')
                            obj.svdUpdateData(iLocBasis).wRef          = zeros(size(obj.prob.ic));
                        else
                            error('When using SVD updating, initref can be ''ic'', ''clustcenters'', or ''zero''');
                        end
                        %Compute values required for online svd updating
                        obj.svdUpdateData(iLocBasis).ColSumRSVtrunc    = sum(RSingVecs(:,1:n),1)';
                        obj.svdUpdateData(iLocBasis).normOneMinusVVt   = norm(ones(size(RSingVecs,1),1) - RSingVecs(:,1:n)*obj.svdUpdateData(iLocBasis).ColSumRSVtrunc,2);
                        obj.svdUpdateData(iLocBasis).bufferedLSingVecs = LSingVecs(:,obj.nY(iLocBasis)+1:n);
                    end
                    
                    if strcmpi(obj.basisUpdate,'border_fast')
                        %Rank of SVD update (including buffer)
                        n = obj.nY(iLocBasis)+obj.svdUpdateData(iLocBasis).bufferSize;
                        
                        %Store the vector to initially reference the
                        %snapshots
                        if strcmpi(obj.snaps.initref,'ic')
                            obj.svdUpdateData(iLocBasis).wRef          = obj.prob.ic;
                        elseif strcmpi(obj.snaps.initref,'clustcenters')
                            obj.svdUpdateData(iLocBasis).wRef          = obj.clusCenter(:,iLocBasis);
                        elseif strcmpi(obj.snaps.initref,'zero')
                            obj.svdUpdateData(iLocBasis).wRef          = zeros(size(obj.prob.ic));
                        else
                            error('When using SVD updating, initref can be ''ic'', ''clustcenters'', or ''zero''');
                        end
                        
                        %Compute values required for online svd updating
                        obj.svdUpdateData(iLocBasis).ColSumRSVtrunc     = sum(RSingVecs(:,1:n),1)';
                        obj.svdUpdateData(iLocBasis).normOneMinusVVt    = norm(ones(size(RSingVecs,1),1) - RSingVecs(:,1:n)*obj.svdUpdateData(iLocBasis).ColSumRSVtrunc,2);
                        obj.svdUpdateData(iLocBasis).bufferedLSingVecs  = LSingVecs(:,obj.nY(iLocBasis)+1:n);
                        
                        %Pre-compute quantities for fast (exact) SVD
                        %updating
                        obj.svdUpdateData(iLocBasis).a = norm(obj.prob.ic);
                        obj.svdUpdateData(iLocBasis).c = obj.prob.ic'*obj.svdUpdateData(iLocBasis).wRef;
                        obj.svdUpdateData(iLocBasis).d = obj.phi0{iLocBasis,1}'*obj.prob.ic;
                    end
                    
                    %Fast approximate border referencing
                    if strcmpi(obj.basisUpdate,'border_fast_approx')
                        %Rank of SVD update (including buffer)
                        n = obj.nY(iLocBasis)+obj.svdUpdateData(iLocBasis).bufferSize;
                        
                        if strcmpi(obj.snaps.initref,'ic')
                            obj.svdUpdateData(iLocBasis).wRef          = obj.prob.ic;
                        elseif strcmpi(obj.snaps.initref,'clustcenters')
                            obj.svdUpdateData(iLocBasis).wRef          = obj.clusCenter(:,iLocBasis);
                        elseif strcmpi(obj.snaps.initref,'zero')
                            obj.svdUpdateData(iLocBasis).wRef          = zeros(size(obj.prob.ic));
                        else
                            error('When using SVD updating, initref can be ''ic'', ''clustcenters'', or ''zero''');
                        end
                        obj.svdUpdateData(iLocBasis).ColSumRSVtrunc     = sum(RSingVecs(:,1:n),1)';
                        obj.svdUpdateData(iLocBasis).normOneMinusVVt    = norm(ones(size(RSingVecs,1),1) - RSingVecs(:,1:n)*obj.svdUpdateData(iLocBasis).ColSumRSVtrunc,2);
                        obj.svdUpdateData(iLocBasis).bufferedLSingVecs  = LSingVecs(:,obj.nY(iLocBasis)+1:n);
                        
                        
                        %Compute truncation error without buffer
                        if obj.nY(iLocBasis) == length(obj.SingVals{iLocBasis})
                            obj.svdUpdateData(iLocBasis).TruncErrorWOBuffer = 0;
                        else
                            obj.svdUpdateData(iLocBasis).TruncErrorWOBuffer   = obj.SingVals{iLocBasis}(obj.nY(iLocBasis)+1);
                        end
                        %Compute truncation error with buffer
                        if n == length(obj.SingVals{iLocBasis})
                            obj.svdUpdateData(iLocBasis).TruncErrorWBuffer = 0;
                        else
                            obj.svdUpdateData(iLocBasis).TruncErrorWBuffer    = obj.SingVals{iLocBasis}(n+1);
                        end
                    end
                    
                    %If SVD updating not used, we know apriori how many
                    %different bases will be used (nBases), so we can
                    %initialize UrLoc.
                    if strcmpi(obj.basisUpdate,'none')
                        obj.UrLoc{iLocBasis,1} = zeros(obj.nY(iLocBasis),1);
                    else
                        %Needs to be computed after the first basis to use is
                        %determined
                    end
                    clear X LSingVecs RSingVecs;
                end
                
                %Pre-compute quantities for fast (exact) SVD
                %updating
                if strcmpi(obj.basisUpdate,'border_fast')
                    for iLocBasis = 1:obj.nBases
                        for jLocBasis = 1:obj.nBases
                            obj.svdUpdateData(iLocBasis).g(jLocBasis) = sqrt(obj.svdUpdateData(iLocBasis).wRef'*obj.svdUpdateData(jLocBasis).wRef);
                            obj.svdUpdateData(iLocBasis).e{jLocBasis} = obj.phi0{iLocBasis,1}'*obj.svdUpdateData(jLocBasis).wRef;
                            obj.svdUpdateData(iLocBasis).F{jLocBasis} = obj.phi0{iLocBasis,1}'*obj.phi0{jLocBasis,1};
                        end
                    end
                end
                obj.numswitch = 0;
                
                disp([' Generated local POD bases of size ',num2str(obj.nY(:)'),', respectively']);
                
                %Pre-Compute Distance Quantities
                obj.precomputeLocalSimulation();
            end
        end %Done
        
        function  [] = computePODadaptsize(obj,fom,trainIC,nYBnds,tol,tolEngy,varargin)
            
            obj.nY = nYBnds(2);
            obj.clusterSnapsOnly(fom,trainIC);
            obj.computePOD('fromClustSnaps',varargin{:});
            
            newNY = 0;
            for p = 1:length(fom)
                projerrNY=obj.computeProjErrorVsNY(fom(p),false,2);
                
                smallbasistol = find(projerrNY < tol,1,'first');
                if isempty(smallbasistol), smallbasistol=0; end;
                newNY = max(newNY,smallbasistol);
                %Make sure newNY in bounds
                newNY = min(newNY,nYBnds(2));
                newNY = max(newNY,nYBnds(1));
            end
            engy = 1-cumsum(obj.SingVals)/sum(obj.SingVals);
            newNY2 = find(engy<tolEngy,1,'first');
            if isempty(newNY2)
                newNY2=nYBnds(2);
            end
            
            obj.nY = max(newNY,newNY2);
            
            obj.computePOD('fromsvd');
            
        end
        
        function  [] = appendPODvecs2base(obj,X,flag,val)
            
            if obj.nBases ~= 1 && strcmpi(obj.snaps.initref,'ic')
                error('This function only works with 1 basis and referencing the IC');
            end
            
            Xa = bsxfun(@minus,X,obj.prob.ic);
            [U,S]=svd(Xa,0);
            S=diag(S);
            N = size(U,2);
            
            if strcmpi(flag,'size')
                nYp = min(N,val);
                if nYp < val
                    warning([num2str(val),' LI vectors not available, using ', num2str(N),' instead.']);
                end
            end
            
            if strcmpi(flag,'projerr') || strcmpi(flag,'projtruncerr')
                Ut_times_Xa = U'*Xa;
                UUt_times_Xa = zeros(size(Xa));
                cwn_Xa = ColumnwiseNorm(Xa,2);
                projerrNY=zeros(N,1);
                for i = 1:N
                    UUt_times_Xa = UUt_times_Xa + U(:,i)*Ut_times_Xa(i,:);
                    projerrNY(i) = max(ColumnwiseNorm(Xa - UUt_times_Xa,2)./cwn_Xa);
                end
                
                tolP = val(1);
                nYp = find(projerrNY < tolP,1,'first');
                if isempty(nYp), nYp = N; end;
            end
            
            if strcmpi(flag,'truncerr') || strcmpi(flag,'projtruncerr')
                tolt = val(end);
                
                engy = 1 - cumsum(S)./sum(S);
                newNY = find(engy < tolt,1,'first');
                if isempty(newNY), newNY = N; end;
                
                if exist('nYp','var')
                    nYp = max(nYp,newNY);
                else
                    nYp = newNY;
                end
            end
            
            oldNY   = obj.nY;
            obj.nY  = obj.nY + nYp;
            obj.phi = [obj.phi, U(:,1:nYp)];
            for i = oldNY+1:obj.nY
                scal = norm(obj.phi(:,i));
                vec = (obj.phi(:,1:i-1)'*(obj.phi(:,i)/scal))';
                vec2 = obj.phi(:,i)/scal - sum(bsxfun(@times,obj.phi(:,1:i-1),vec),2);
                obj.phi(:,i) = scal*vec2/norm(vec2,2);
            end
            
            r=rank(obj.phi);
            if r < size(obj.phi,2)
                obj.phi = obj.phi(:,1:r);
            end
        end
        
        function [] = modifyInitialCondition(obj,probobj)
            
            obj.snaps.initref = 'somevec';
            obj.prob.modifyInitialCondition(probobj);
            if (obj.nBases>1 && strcmpi(obj.basisUpdate,'border_fast') )
                for iLocBasis = 1:obj.nBases
                    obj.svdUpdateData(iLocBasis).a = norm(probobj.ic,2);
                    obj.svdUpdateData(iLocBasis).c = probobj.ic'*obj.svdUpdateData(iLocBasis).wRef;
                    obj.svdUpdateData(iLocBasis).d = obj.phi0{iLocBasis,1}'*probobj.ic;
                end
            end
            obj.precomputeLocalSimulation();
        end
        
        function [] = modifyTimeProperties(obj,fomobj)
            
            obj.time = fomobj.time;
            % obj.prob.config.modifyTimeProperties(fomobj.prob.config);
            obj.prob.config.setTimeObj('T',fomobj.prob.config.time.T);
            obj.prob.config.setTimeObj('dt',fomobj.prob.config.time.dt);
            obj.prob.config.setTimeObj('nstep',fomobj.prob.config.time.nstep);
            obj.prob.config.setTimeObj('quiet',fomobj.prob.config.time.quiet);
            obj.TimeScheme = fomobj.TimeScheme;
        end
        
        function [] = clearOnlineUpdateProperties(obj)
            
            obj.ROBcompon = [];
            obj.numswitch = 0;
            obj.prevSwitchState = [];
            obj.wr_switch = [];
            obj.phi = [];
            obj.UrLoc = [];
            obj.LocBasisHist = [];
            
        end
        
        function  [] = clusterPODprecompLoc(obj,trainFOM,trainIC,flag,svdOptions,rootFname)
            %This function combines the clustering and POD steps.
            
            obj.clusterSnapsOnly(trainFOM,trainIC);
            if nargin == 4
                obj.computePOD(flag);
            elseif nargin == 5
                obj.computePOD(flag,svdOptions);
            elseif nargin == 6
                obj.computePOD(flag,svdOptions,rootFname);
            end
            
        end %Done
        
        function  [] = clusterSnapsOnly(obj,trainFOM,trainIC)
            %This function clusters the snapshots from the FOM objects
            %(snapshots extracted according to the ROM file).  The
            %snapshots are then saved in a binary file named
            %[rom.folderRoot,'procSnaps',num2str(basisNum),'.bin'].
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj       - ROM object
            %trainFOM  - array of FOM objects that will contain
            %trainIC   - training initial condition
            %
            %Outputs:
            %--------
            %There are no outputs.
            %--------------------------------------------------------------
            
            %Need fom array to use this function
            if nargin < 2
                error('ROM.computePOD requires fom array input');
            end
            
            if nargin < 3 || isempty(trainIC)
                fprintf('Warning: Did not specify trainIC in call to computePOD() in workflow file! \n');
                fprintf('If operations requiring trainIC are required (i.e. using initial condition to reference snapshots;\n');
                fprintf('the initial condition of the first FOM in the array trainFOM will be used.\n');
                
                trainIC = trainFOM(1).prob.ic;
            end
            
            if obj.nBases == 1
                %Extract processed snapshots from FOM objects.
                %We take into account the fact that two
                %consecutive raw snapshots may not correspond to
                %consecutive snapshots in a simulation.
                ProcessedSnaps = [];
                for i = 1:length(trainFOM)
                    %                     %Extract every raw snapshot from the current FOM.
                    %                     AllRawStates = obj.extractRawSnapshots(trainFOM(i),0);
                    %Extract all desired snapshots (i.e. only in interval,
                    %from distribution, etc)
                    [DesiredRawStates,indDesired] = obj.extractRawSnapshots(trainFOM(i),1);
                    %Extract the "previous state" snapshots with respect to the desired states
                    PrevSnaps = obj.extractRawSnapshots(trainFOM(i),0,indDesired-1);
                    %PrevSnaps = AllRawStates(:,indDesired-1); clear AllRawStates;
                    
                    if strcmpi(obj.snaps.initref,'clustcenters')
                        ProcessedSnaps = [ProcessedSnaps,snapshotCollection(DesiredRawStates,'zero',trainIC,PrevSnaps)]; %#ok<AGROW>
                    else
                        ProcessedSnaps = [ProcessedSnaps,snapshotCollection(DesiredRawStates,obj.snaps.initref,trainIC,PrevSnaps)]; %#ok<AGROW>
                    end
                end
                %If we are referencing the cluster centers, we must
                %complete the above loop because we need to know all
                %snapshots before we can know the center
                if strcmpi(obj.snaps.initref,'clustcenters')
                    ProcessedSnaps = snapshotCollection(ProcessedSnaps,obj.snaps.initref,[],[],mean(ProcessSnaps,2));
                end
                nSnaps = size(ProcessedSnaps,2);
                obj.snaps.nSnapsPerBase = nSnaps;
                
                %Write clustered snapshots to file
                fid = fopen([obj.fileRoot,'clusteredSnaps1.bin'],'wb');
                fwrite(fid,size(ProcessedSnaps)','double');
                fwrite(fid,ProcessedSnaps,'double');
                fclose(fid);
            else
                %Extract raw snapshots and processed snapshots from FOM
                %objects.  We take into account the fact that two
                %consecutive raw snapshots may not correspond to
                %consecutive snapshots in a simulation.
                RawSnaps = []; PrevSnaps = []; FOMNUM = []; SVNUM = [];
                for i = 1:length(trainFOM)
                    %Extract every raw snapshot from the current FOM.
                    AllRawStates = obj.extractRawSnapshots(trainFOM(i),0);
                    %Extract all desired snapshots (i.e. only in interval,
                    %from distribution, etc)
                    [DesiredRawStates,indDesired] = obj.extractRawSnapshots(trainFOM(i),1);
                    RawSnaps  = [RawSnaps, DesiredRawStates]; %#ok<AGROW>
                    %Extract the "previous state" snapshots with respect to
                    %the desired states
                    PrevSnaps = [PrevSnaps, AllRawStates(:,indDesired-1)]; %#ok<AGROW>
                    %Keep track of which FOM simulation each snapshot comes
                    %from and which state within each FOM the snapshot
                    %corresponds to.
                    FOMNUM = [FOMNUM, i*ones(1,size(DesiredRawStates,2))]; %#ok<AGROW>
                    SVNUM = [SVNUM, indDesired]; %#ok<AGROW>
                end
                nSnaps = size(RawSnaps,2);
                
                clear DesiredRawStates indDesired AllRawStates;
                
                %Cluster the snapshots in nLocalBases clusters
                switch obj.clustering
                    case 'k_means'
                        [closClusIndex,obj.clusCenter,additionalPoints] = kMeansClustering(obj.nBases,RawSnaps,obj.addElemTol);
                    case 'a_priori'
                        [closClusIndex,obj.clusCenter,additionalPoints] = aPrioriClustering(obj.nBases,RawSnaps,obj.addElemTol);
                end
                
                %Open binary files for writing
                for i = 1:obj.nBases
                    fid(i) = fopen([obj.fileRoot,'clusteredSnaps',num2str(i),'.bin'],'wb');
                end
                
                %Store the time steps corresponding to each cluster if the
                %basis is built from only 1 training simulation.
                for i = 1:obj.nBases
                    obj.ClusterIndex.true{i}  = find(closClusIndex==i);
                    if obj.addElemTol > 0
                        obj.ClusterIndex.added{i} = additionalPoints(i).Indices;
                    end
                end
                obj.ClusterIndex.fomnum = FOMNUM;
                obj.ClusterIndex.svnum = SVNUM;
                clear FOMNUM SVNUM;
                
                
                nPODLoc = zeros(obj.nBases,1);
                for iSnap=1:nSnaps
                    %Process each snapshot and write it to the appropriate
                    %cluster file.  Also, keep track of the number of
                    %snapshots in each cluster.
                    fwrite(fid(closClusIndex(iSnap,1)),snapshotCollection(RawSnaps(:,iSnap),obj.snaps.initref,trainIC,PrevSnaps(:,iSnap),obj.clusCenter(:,closClusIndex(iSnap,1))),'double');
                    nPODLoc(closClusIndex(iSnap,1)) = nPODLoc(closClusIndex(iSnap,1)) + 1;
                end
                
                %Process and write the snapshots that are "shared" to the
                %appropriate files.
                if obj.addElemTol > 0
                    for iLocBasis=1:obj.nBases
                        Indices = additionalPoints(iLocBasis).Indices;
                        for iAddElem = 1:length(Indices)
                            fwrite(fid(iLocBasis),snapshotCollection(RawSnaps(:,Indices(iAddElem)),obj.snaps.initref,trainIC,PrevSnaps(:,Indices(iAddElem)),obj.clusCenter(:,iLocBasis)),'double');
                            nPODLoc(iLocBasis) = nPODLoc(iLocBasis) + 1;
                        end
                    end
                end
                clear RawSnaps PrevSnaps;
                
                %Include the number of snapshots in the file.
                for i = 1:obj.nBases
                    %%Need to insert the size AT THE BEGINNING of the file.
                    %%Do this by reading the file to tmp and then overwriting it.
                    fclose(fid(i));
                    fid(i) = fopen([obj.fileRoot,'clusteredSnaps',num2str(i),'.bin'],'r');
                    tmp = fread(fid(i),obj.ndof*nPODLoc(i),'double');
                    fclose(fid(i));
                    fid(i) = fopen([obj.fileRoot,'clusteredSnaps',num2str(i),'.bin'],'wb');
                    fwrite(fid(i),[obj.ndof;nPODLoc(i)],'double');
                    fwrite(fid(i),tmp,'double');
                    fclose(fid(i));
                end
            end
        end %Done
        
        function  [] = visualizeClustering(obj,subplotshape)
            %This function visualizes the clustering.  If the snapshots
            %were generated from multiple FOMs, the snapshots for each FOM
            %are put in different subplots.
            
            %Determine number of FOMs used to generate snapshots
            numFom  = length(unique(obj.ClusterIndex.fomnum));
            %If not subplot layout specified, just use a single column
            if nargin == 1
                subplotshape = [numFom,1];
            end
            
            for j = 1:numFom
                subplot(subplotshape(1),subplotshape(2),j);
                title(['FOM ',num2str(j)],'interpreter','latex'); hold on;
                %Determine the snapshot indices that correspond to the
                %current FOM
                ind = find(obj.ClusterIndex.fomnum == j);
                for i = 1:obj.nBases
                    %Determine the "true" indicies that correspond to the
                    %current FOM
                    indt = intersect(obj.ClusterIndex.true{i},ind);
                    %Determine the "added" indicies that correspond to the
                    %current FOM
                    inda = intersect(obj.ClusterIndex.added{i},ind);
                    
                    %Plot all indicies in cluster (true and added) in large
                    %black markers and plot only the true clustered indices
                    %in small yellow markers
                    plot(obj.ClusterIndex.svnum(1,union(indt,inda))-1,i*ones(1,length(union(indt,inda))),'ko','markerfacecolor','k','markersize',6); hold on;
                    plot(obj.ClusterIndex.svnum(1,indt)-1,i*ones(1,length(indt)),'yo','markerfacecolor','y','markersize',3); hold on;
                end
            end
        end %Done
        
function  [] = executeModel(obj,restart,augment)
    %This function performs the time-stepping for the ROM
    %simulation defined in this object.
    %--------------------------------------------------------------
    %Inputs:
    %-------
    %obj     - ROM object
    %
    %Outputs:
    %--------
    %There are no outputs.
    %--------------------------------------------------------------
    obj.trunc = min(obj.trunc, obj.nY);

    obj.numExecute = obj.numExecute+1;

    if nargin == 3 && ~isempty(augment)
        obj.augment = augment;
    else
        obj.augment=false;
    end

    obj.killflag = false;
    %nstep = obj.prob.config.time.nstep;
    nstep = obj.time.nstep;

    %Determine whether to run simulation with precomputations being
    %carried out before hand
    if obj.precompFlag
        obj.sv=zeros(obj.nY,nstep+1);
        obj.sv(:,1) = obj.prob.ic;

        obj.curr_param = obj.prob.p;

        tROMstart = tic;
        for i_t = 1:nstep %Loop over each time step to determine all state vectors
            if ~obj.time.quiet && rem(i_t,round(nstep*0.1)) == 0
                %Generate output so user can see progress
                fprintf('%4.0f of %4.0f Timesteps Complete (%2.0f%%)\n',i_t,nstep,100*i_t/nstep)
            end

            if obj.TimeScheme.explicit
                obj.sv(:,i_t+1) = obj.TimeScheme.ExplicitStep(obj.prob,obj.sv(:,i_t),obj.time.T(1)+obj.TimeScheme.dt*i_t,true);
                if sum(isnan(obj.sv(:,i_t+1)))>0
                    fprintf('NaN encountered on time step %i!\n',i_t);
                    obj.killflag=true;
                    return;
                end
                continue;
            end

            %Set current time iteration (so model can keep track)
            obj.cTimeIter = i_t;

            %Use Newton's Method to solve nonlinear system of equations
            %and store: state vector, residual and jacobian snapshots,
            %and number of newton iterations required
            if obj.nBases>1
                obj.LocBasisHist(obj.cTimeIter,1) = computeDistanceToCenters(obj,obj.UrLoc,obj.cTimeIter);
            end

            obj.NewtonRaphson_PrecompGal();
            %                     obj.NewtonRaphson();

            if obj.killflag, warning('Encountered NaN at Time Step %i. Simulation Killed',obj.cTimeIter); return; end;
        end
        obj.ontime = toc(tROMstart); %Record simulation time
        return;
    end

    %Determine restart status and prepare to start simulation
    if nargin == 2 && ~isempty(restart) && restart == 1
        OldTime = obj.ontime;
        firstStep = obj.cTimeIter-1;
    else
        OldTime = 0;
        firstStep = 1;

        obj.numres = zeros(obj.nBases,1);
        obj.numswitch = 0;
        %Set up state vector and store the initial condition.
        obj.sv = zeros(obj.ndof,obj.time.nstep+1);
        obj.sv(:,1) = obj.prob.ic;
    end

    %Set the current parameter (from prob object) to be used in
    %NAND to know which parameter values have been used to generate
    %the current state vectors
    obj.curr_param = obj.prob.p;
    obj.phi= obj.phi(:,1:obj.trunc);
    obj.openCloseResJacFiles('open');
    switch obj.nBases
        case 1
            tROMstart = tic;
            for i_t = firstStep:nstep %Loop over each time step to determine all state vectors
                if ((obj.printLevel>0)&&(rem(i_t,round(nstep*0.1)) == 0)) || (obj.printLevel>1)
                    fprintf('-------------------------- Time Step %4i of %4i (%2i%%) -------------------------- \n',i_t,nstep,ceil(100*i_t/nstep));
                end
                %                         if ~obj.time.quiet && rem(i_t,round(nstep*0.1)) == 0
                %                             %Generate output so user can see progress
                %                             fprintf('%4.0f of %4.0f Timesteps Complete (%2.0f%%)\n',i_t,nstep,100*i_t/nstep)
                %                         end

                %Set current time iteration (so model can keep track)
                obj.cTimeIter = i_t;

                if ~isempty(obj.time.cfl)
                    obj.time.curr_cfl = obj.time.cfl(i_t,obj);
                    obj.time.dt = obj.time.curr_cfl*obj.prob.computeCFL(obj.sv(:,i_t));
                    %                           obj.time.dt = obj.time.cfl(i_t,obj.sv(:,i_t))*obj.prob.computeCFL(obj.sv(:,i_t));
                    obj.TimeScheme.setDT(obj.time.dt);
                end


                obj.chooseSolver()

                if i_t==nstep
                    NewFrBasis=[];
                    NewFlBasis=[];
                    soln=obj.sv;
                    for jj=1:size(soln,2)

                        U = reshape(soln(:,jj),3,obj.prob.nVol);
                        [rho,u,P,c] = obj.prob.conservativeToPrimitive(U(:,2:end-1));
                        rho = [U(1,1),rho,U(1,end)]; %Density
                        u   = [U(2,1),u,U(2,end)]; %Velocity
                        P   = [U(3,1),P,U(3,end)]; %Pressure
                        c   = [sqrt(obj.prob.gamma*P(1)/rho(1)),c,sqrt(obj.prob.gamma*P(end)/rho(end))]; %Speed of sound
                        e   = [P(1)/(obj.prob.gamma-1)+rho(1)*u(1)^2/2,U(3,2:end-1),P(end)/(obj.prob.gamma-1)+rho(end)*u(end)^2/2];
                        dc_cons = [0.5*obj.prob.gamma./(c.*rho).*(0.5*(obj.prob.gamma-1)*u.*u - P./rho);...
                            -0.5*obj.prob.gamma*(obj.prob.gamma-1)*u./(rho.*c);...
                            0.5*obj.prob.gamma*(obj.prob.gamma-1)./(rho.*c)]';
                        [roeF,droeF] = obj.prob.roeFlux(rho,u,P,c,e,dc_cons);
                        [Q,dQ]=obj.prob.forceTerm(u,P);

                        dUdV=[1,0,0;u(1),rho(1),0;0.5*u(1)*u(1),rho(1)*u(1),1/(obj.prob.gamma-1)];
                        droeF(:,1:3,1)=droeF(:,1:3,1)*dUdV;
                        dUdV=[1,0,0;u(end),rho(end),0;0.5*u(end)*u(end),rho(end)*u(end),1/(obj.prob.gamma-1)];
                        droeF(:,4:6,end)=droeF(:,4:6,end)*dUdV;

                        roeF1(:,jj) = roeF(1,:)';
                        roeF2(:,jj) = roeF(2,:)';
                        roeF3(:,jj) = roeF(3,:)';

                        Fright(:,jj)= reshape(bsxfun(@times,roeF(:,2:end),obj.prob.S(3:end-1)), size(roeF,1)*size(roeF(:, 2:end),2),1);
                        Fleft(:,jj) = reshape(-bsxfun(@times,roeF(:,1:end-1),obj.prob.S(2:end-2)),size(roeF,1)*size(roeF(:,2:end),2),1);
                        forceQ(:,jj)= reshape(Q, size(Q,1)*size(Q,2), 1);
                        J2L=zeros(3*(obj.prob.nVol-2), 3*obj.prob.nVol);
                        J2R=J2L;
                        for kl= 2:obj.prob.nVol-1
                            i=kl-1;
                            J2L(3*(i-1)+1:3*i,3*(kl-2)+1:3*kl) = [-obj.prob.S(kl)*droeF(:,1:3,kl-1), -obj.prob.S(kl)*droeF(:,4:6,kl-1)];
                            J2R(3*(i-1)+1:3*i,3*(kl-1)+1:3*(kl+1)) = [obj.prob.S(kl+1)*droeF(:,1:3,kl), obj.prob.S(kl+1)*droeF(:,4:6,kl)];
                        end


                        if jj==1

                            sq=Q(:,2:end-1)*(obj.prob.SVol(2:end-1).*obj.prob.dx(2:end-1))';
                            sdq=zeros(3,obj.trunc);
                            Phi = obj.phi(:,1:obj.trunc);
                            for i=2:obj.prob.nVol-1
                                sdq=sdq+squeeze(dQ(:,:,i))*Phi(3*i-2:3*i,:)*(obj.prob.SVol(i).*obj.prob.dx(i));
                            end
                            sq_true = sq;
                            sdq_true = sdq;
                            save sq_true sq_true
                            save sdq_true sdq_true

                            J2L_true=J2L;
                            J2R_true=J2R;
                            DforceQ=dQ;
                            romdroeF= droeF;
                            Q2_true=Q(2,:);
                            
                            save romroeF roeF
                            save Fright Fright
                            save Fleft Fleft
                            save forceQ forceQ

                            save DforceQ DforceQ
                            save J2R_true J2R_true
                            save Q2_true Q2_true
                            save J2L_true J2L_true
                            save romdroeF romdroeF
                        end

                        NewFrBasis = [NewFrBasis, Fright(:,jj), J2R*obj.phi];
                        NewFlBasis = [NewFlBasis, Fleft(:,jj), J2L*obj.phi];
%                         NewQBasis  = [NewQBasis, forceQ(:,jj), dQ*obj.phi];
                        romroeF(:,jj)=roeF(:);
                        if imag(Fright(:,jj))~=0
                            disp(['righ flux is complex for snapshot number', num2str(jj)])
                        elseif imag(Fleft(:,jj))~=0
                            disp(['left flux is complex for snapshot number', num2str(jj)])
                        elseif imag(forceQ)~=0
                            disp(['force term is complex for snapshot number', num2str(jj)])
                        end
                    end

                    %save romroeF romroeF
                    %keyboard

                    [ur2,sr2,~] = svd(NewFrBasis,0);
                    [ul2,sl2,~] = svd(NewFlBasis,0);
                    [ur,sr,~]   = svd(Fright,0);
                    [ul, sl, ~] = svd(Fleft,0);
                    [uQ,sQ,~]   = svd(forceQ,0);
                    [uroe1,s1,~] = svd(roeF1,0);
                    [uroe2,s2,~] = svd(roeF2,0);
                    [uroe3,s3,~] = svd(roeF3,0);

                    % mFr=min(find(cumsum(diag(sr2))/sum(diag(sr2))>0.9999));
                    % mF=min(find(cumsum(diag(sr))/sum(diag(sr))>0.9999));
                    % mQ=min(find(cumsum(diag(sQ))/sum(diag(sQ))>0.9999));
                    % mr=min(find(cumsum(diag(s1))/sum(diag(s1))>0.9999));

                    %keyboard

                    obj.phiFrnew=ur2;%(:,1:mFr);
                    obj.phiFlnew=ul2;%(:,1:mFr);
                    obj.phiFright=ur;%(:,1:mF); %99.9%
                    obj.phiFleft=ul;%(:,1:mF);  %99.9%
                    obj.phiRoeF1=uroe1;%(:,1:mr);
                    obj.phiRoeF2=uroe2;%(:,1:mr);
                    obj.phiRoeF3=uroe3;%(:,1:mr);
                    obj.phiQ=uQ;%(:,1:mQ);
                end


                if obj.killflag, warning('Encountered NaN at Time Step %i. Simulation Killed',obj.cTimeIter); return; end;

            end
            obj.ontime = toc(tROMstart)+OldTime; %Record simulation time
        otherwise
            tROMstart = tic;
            for i_t = firstStep:nstep %Loop over each time step to determine all state vectors
                %                         if ~obj.time.quiet && rem(i_t,round(nstep*0.1)) == 0
                %                             %Generate output so user can see progress
                %                             fprintf('%4.0f of %4.0f Timesteps Complete (%2.0f%%)\n',i_t,nstep,100*i_t/nstep)
                %                         end

                %Set current time iteration (so model can keep track)
                obj.cTimeIter = i_t;

                if ~isempty(obj.time.cfl)
                    obj.time.dt = obj.time.cfl(i_t,obj)*obj.prob.computeCFL(obj.sv(:,i_t));
                    %                             obj.time.dt = obj.time.cfl(i_t,obj.sv(:,i_t))*obj.prob.computeCFL(obj.sv(:,i_t));
                    obj.TimeScheme.setDT(obj.time.dt);
                end

                %Use Newton's Method to solve nonlinear system of equations
                %and store: state vector, residual and jacobian snapshots,
                %and number of newton iterations required
                %                         switch obj.clustering
                %                             case 'kmeans'
                obj.LocBasisHist(obj.cTimeIter,1) = computeDistanceToCenters(obj,obj.UrLoc,obj.cTimeIter);
                %                             case 'a_priori'
                %                                 obj.LocBasisHist(obj.cTimeIter,1) = floor(obj.nBases*obj.cTimeIter/nstep)+1;
                %                         end
                cbase=obj.LocBasisHist(obj.cTimeIter,1);

                if ((obj.printLevel>0)&&(rem(i_t,round(nstep*0.1)) == 0)) || ((obj.printLevel>1)&&(rem(i_t,round(nstep*0.05)) == 0)) || (obj.printLevel>2)
                    fprintf('-------------------------- Time Step %4i of %4i (%2i): Basis %i -------------------------- \n',i_t,nstep,100*i_t/nstep,cbase);
                end

                %Determine initial ROB components if using fast
                %updating
                if obj.cTimeIter == 1 && strcmpi(obj.basisUpdate,'border_fast')
                    obj.initializeROBcompon(cbase);
                end
                obj.NewtonRaphsonLocal();
                if obj.killflag, warning('Encountered NaN at Time Step %i. Simulation Killed',obj.cTimeIter); return; end;

                %Check for steady convergence
                if obj.time.steadyconverge
                    if norm(obj.sv(:,i_t+1)-obj.sv(:,i_t)) < obj.time.steadyconverge*norm(obj.sv(:,i_t+1))
                        obj.ontime = toc(tROMstart)+OldTime;
                        obj.time.nstep=i_t;
                        obj.time.T(2) = obj.time.T(1) + i_t*obj.time.dt;
                        obj.sv(:,i_t+2:end)=[];
                        obj.openCloseResJacFiles('close');
                        obj.newt.avgIter = mean(obj.newt.iter);
                        return;
                    end
                end
            end
            obj.ontime = toc(tROMstart)+OldTime; %Record simulation time
    end
    obj.openCloseResJacFiles('close');
    obj.newt.avgIter = mean(obj.newt.iter);
end %Done
        
        function  [] = computeResNormFromSV(obj,flag,nType)
            %This function computes the nType-norm of the residual (from
            %rom.sv).
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj   - ROM object
            %flag  - string indicating whether to return a relative or
            %        absolute norm ('rel' vs 'abs').
            %nType - positive scalar indicating norm to use
            %--------------------------------------------------------------
            
            %If nType not specific, use 2-norm as default
            if nargin == 2
                nType = 2;
            elseif nargin < 2 || nargin > 3
                error('computeResNormFromSV.m requires the ROM object and at least 1 additional input');
            end
            
            obj.newt.resConv = zeros(obj.time.nstep,1);
            
            %Compute the residual nType norm.  SPECIFIC TO BACKWARD EULER!
            for i = 1:obj.time.nstep
                t = obj.time.T(1) + i*obj.time.dt;
                Rtemp = obj.prob.ResJac(obj.sv(:,i+1),t);
                resVec = Rtemp - (1/obj.time.dt)*(obj.sv(:,i+1) - obj.sv(:,i));
                if strcmpi(flag(1),'r')
                    %Notice that we are normalizing by the residual at the
                    %beginning of the Newton iteration assuming the initial
                    %guess is the previous time step.
                    obj.newt.resConv(i) = norm(resVec,nType)/norm(obj.prob.ResJac(obj.sv(:,i),t),nType);
                elseif strcmpi(flag(1),'a')
                    obj.newt.resConv(i) = norm(resVec,nType);
                else
                    error('flag must be ''abs'' or ''rel'' to indicate absolute or relative residual norm');
                end
            end
        end %Done
        
        function  [] = computeResNormBound(obj,fomobj,flag,nType)
            %This function computes the residual bound for each time step.
            %SPECIFIC TO BACKWARD EULER and ref_init and rom time struct =
            %fom time struct.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj    - ROM object
            %fomobj - FOM object
            %flag   - string indicating whether to return a relative or
            %         absolute norm ('rel' vs 'abs').
            %nType  - positive scalar indicating norm to use.  If not
            %         specified, 2-norm is the default
            %
            %Outputs:
            %--------
            %There are no outputs.
            %--------------------------------------------------------------
            
            if ~(strcmpi(class(obj.TimeScheme),'BackwardEuler') && strcmpi(obj.snaps.initref,'ic') && strcmpi(obj.basisUpdate,'none'))
                disp('A priori residual error bound only works for BACKWARD EULER time scheme and ic snapshot reference.  Exiting...');
                return;
            end
            
            if ~(fomobj.time.nstep == obj.time.nstep && sum(fomobj.time.T == obj.time.T) == 2)
                %If FOM and ROM do not have same time structure, warn
                %user and don't try to compute error analysis.
                disp('Warning.....Your FOM and ROM time structure are not the same.  Residual error estimation is not supported for this case.');
                return;
            end
            
            if obj.nBases ~= 1
                disp('Only works for single basis.');
                return;
            end
            
            %If nType not specific, use 2-norm as default
            if nargin == 3
                nType = 2;
            elseif nargin < 3 || nargin > 4
                error('computeResNormBound.m requires the ROM object and at least 2 additional inputs');
            end
            
            obj.resBound = zeros(obj.time.nstep,1);
            %Compute a bound on the residual.
            for i = 1:obj.time.nstep
                t = obj.time.T(1) + i*obj.time.dt;
                Xfull = [fomobj.sv(:,i) - obj.prob.ic, fomobj.sv(:,i+1) - obj.prob.ic];
                
                w_tm1 = obj.prob.ic + obj.phi(:,:)*(obj.phi(:,:)'*Xfull(:,1));
                w = obj.prob.ic + obj.phi(:,:)*(obj.phi(:,:)'*Xfull(:,2));
                
                if strcmpi(flag(1),'r')
                    obj.resBound(i,1) = norm(obj.prob.ResJac(w,t) - (1/obj.time.dt)*(w - w_tm1),nType)/norm(obj.prob.ResJac(w_tm1,t),nType);
                elseif strcmpi(flag(1),'a')
                    obj.resBound(i,1) = norm(obj.prob.ResJac(w,t) - (1/obj.time.dt)*(w - w_tm1),nType);
                else
                    error('flag must be ''abs'' or ''rel'' to indicate absolute or relative residual norm');
                end
            end
        end %Done
        
        function  [rB] = computeResNormBoundAtTimeVsBasisSize(obj,tsteps,fomobj,flag,nType)
            %This function computes the residual bound for each time step
            %and a range of basis sizes. SPECIFIC TO BACKWARD EULER and
            %ref_init and rom time struct = fom time struct.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - ROM object
            %tsteps  - time steps at which norm should be computed
            %fomobj  - FOM object
            %flag    - string indicating whether to return a relative or
            %          absolute norm ('rel' vs 'abs')
            %
            %Outputs:
            %--------
            %There are no outputs.
            %--------------------------------------------------------------
            
            rB = [];
            
            if ~strcmpi(class(obj.TimeScheme),'BackwardEuler') || ~(strcmpi(obj.snaps.initref,'ic')&&strcmpi(obj.basisUpdate,'none'))
                disp('A priori residual error bound only works for BACKWARD EULER time scheme and ic snapshot reference.  Exiting...');
                return;
            end
            
            if ~(fomobj.time.nstep == obj.time.nstep && sum(fomobj.time.T == obj.time.T) == 2)
                %If FOM and ROM do not have same time structure, warn
                %user and don't try to compute error analysis.
                disp('Warning.....Your FOM and ROM time structure are not the same.  Residual error estimation is not supported for this case.');
                return;
            end
            
            if obj.nBases ~= 1
                return;
            end
            
            %If nType not specific, use 2-norm as default
            if nargin == 4
                nType = 2;
            elseif nargin < 4 || nargin > 5
                error('computeResNormBound.m requires the ROM object and at least 2 additional inputs');
            end
            
            rB = zeros(length(tsteps),obj.nY);
            statevecs    = bsxfun(@minus,fomobj.sv(:,tsteps),obj.prob.ic);
            pstatevecs   = bsxfun(@minus,fomobj.sv(:,tsteps-1),obj.prob.ic);
            %Compute a bound on the residual.
            for i = 1:length(tsteps)
                w_tm1 = obj.prob.ic;
                w = obj.prob.ic;
                for j = 1:obj.nY
                    t = obj.time.T(1) + tsteps(i)*obj.time.dt;
                    
                    w_tm1 = w_tm1 + obj.phi(:,j)*(obj.phi(:,j)'*pstatevecs(:,i));
                    w     = w + obj.phi(:,j)*(obj.phi(:,j)'*statevecs(:,i));
                    
                    if strcmpi(flag(1),'r')
                        rB(i,j) = norm(obj.prob.ResJac(w,t) - (1/obj.time.dt)*(w - w_tm1),nType)/norm(obj.prob.ResJac(w_tm1,t),nType);
                    elseif strcmpi(flag(1),'a')
                        rB(i,j) = norm(obj.prob.ResJac(w,t) - (1/obj.time.dt)*(w - w_tm1),nType);
                    else
                        error('flag must be ''abs'' or ''rel'' to indicate absolute or relative residual norm');
                    end
                end
            end
        end %Done
        
        function  [] = computeProjError(obj,fomobj,plotFlag,normType,normalizeFlag)
            %Computes the projection error for each basis and store results
            %in obj.projerr.  If plotFlag = true, plots will be generated.
            
            %Default norm is the 2 norm
            if nargin < 4, normType = 2; end;
            if nargin < 5, normalizeFlag = false; end;
            
            obj.projerr = zeros(fomobj.time.nstep+1,obj.nBases);
            if obj.nBases == 1
                %Compute projection error for single basis case
                vecs = bsxfun(@minus,fomobj.sv,obj.prob.ic);
                obj.projerr = ColumnwiseNorm(vecs - obj.phi*(obj.phi'*vecs),normType);
                if normalizeFlag
                    obj.projerr = obj.projerr./ColumnwiseNorm(fomobj.sv,2);
                end
                %Plot results
                if plotFlag
                    %figure;
                    plot(linspace(fomobj.time.T(1),fomobj.time.T(2),fomobj.time.nstep+1),obj.projerr,...
                        'k','linewidth',2);
                end
                xlabel('time','interpreter','latex'); ylabel('Projection Error','interpreter','latex');
                return;
            end
            
            %Copmute projection error for multiple basis case
            for i = 1:obj.nBases
                vecs = bsxfun(@minus,fomobj.sv,obj.prob.ic);
                obj.projerr(:,i) = ColumnwiseNorm(vecs - obj.phi0{i}*(obj.phi0{i}'*vecs),normType);
                if normalizeFlag
                    obj.projerr(:,i) = obj.projerr(:,i)./ColumnwiseNorm(fomobj.sv,2)';
                end
            end
            
            %Plot projection error for each local basis
            if plotFlag
                figure; ax = axes;
                cmap = colormap(ax);
                
                %Determine color indexing
                ind = round(linspace(1,size(cmap,1),obj.nBases));
                str = [];
                for i = 1:obj.nBases
                    %Plot projection error vs. time
                    plot(linspace(fomobj.time.T(1),fomobj.time.T(2),fomobj.time.nstep+1),obj.projerr(:,i),...
                        'linewidth',2,'color',cmap(ind(i),:));
                    hold on;
                    
                    str{i} = ['Basis ',num2str(i)];
                end
                legend(str{:});
                xlabel('time','interpreter','latex'); ylabel('Projection Error','interpreter','latex');
            end
            
        end %Done
        
        function  [projerrNY] = computeProjErrorVsNY(obj,fomobj,plotFlag,normType)
            %Computes the projection error for each basis and store results
            %in obj.projerr.  If plotFlag = true, plots will be generated.
            
            %Default norm is the 2 norm
            if nargin < 4
                normType = 2;
            end
            
            projerrNY = zeros(obj.nY,obj.nBases);
            if obj.nBases == 1
                %Compute projection error for single basis case
                vecs = bsxfun(@minus,fomobj.sv,obj.prob.ic);
                phiT_times_vecs = obj.phi'*vecs;
                phiphiT_times_vecs = zeros(size(vecs));
                cwn_fomsv= ColumnwiseNorm(fomobj.sv,2);
                for i = 1:obj.nY
                    phiphiT_times_vecs = phiphiT_times_vecs + obj.phi(:,i)*phiT_times_vecs(i,:);
                    projerrNY(i) = max(ColumnwiseNorm(vecs - phiphiT_times_vecs,normType)./cwn_fomsv);
                end
                %Plot results
                if plotFlag
                    figure;
                    plot(1:obj.nY,projerrNY,'k','linewidth',2);
                    xlabel('time','interpreter','latex'); ylabel('Projection Error','interpreter','latex');
                end
                return;
            end
            
            %             %Copmute projection error for multiple basis case
            %             for i = 1:obj.nBases
            %                 vecs = bsxfun(@minus,fomobj.sv,obj.prob.ic);
            %                 obj.projerr(:,i) = ColumnwiseNorm(vecs - obj.phi0{i}*(obj.phi0{i}'*vecs),normType);
            %                 if normalizeFlag
            %                     obj.projerr(:,i) = obj.projerr(:,i)./ColumnwiseNorm(fomobj.sv,2)';
            %                 end
            %             end
            %
            %             %Plot projection error for each local basis
            %             if plotFlag
            %                 figure; ax = axes;
            %                 cmap = colormap(ax);
            %
            %                 %Determine color indexing
            %                 ind = round(linspace(1,size(cmap,1),obj.nBases));
            %                 str = [];
            %                 for i = 1:obj.nBases
            %                     %Plot projection error vs. time
            %                     plot(linspace(fomobj.time.T(1),fomobj.time.T(2),fomobj.time.nstep+1),obj.projerr(:,i),...
            %                         'linewidth',2,'color',cmap(ind(i),:));
            %                     hold on;
            %
            %                     str{i} = ['Basis ',num2str(i)];
            %                 end
            %                 legend(str{:});
            %                 xlabel('time','interpreter','latex'); ylabel('Projection Error','interpreter','latex');
            %             end
        end
        
        function  [] = setPhi(obj,cbase)
            %Sets phi to the master basis specified by cbase
            
            obj.phi = obj.phi0{cbase};
        end %Done
        
        function  [] = setNumSwitch(obj,num)
            %Sets numswitch to the appropriate number
            
            obj.numswitch=num;
        end %Done
        
        function  [newobj] = hardCopy(obj,idnum)
            % Construct a new object based on a deep copy of an old object
            % of this class by copying properties over.
            
            newobj=ROM();
            
            props = properties(obj);
            for i = 1:length(props)
                % Use Dynamic Expressions to copy the required property.
                % For more info on usage of Dynamic Expressions, refer to
                % the section "Creating Field Names Dynamically" in:
                % web([docroot '/techdoc/matlab_prog/br04bw6-38.html#br1v5a9-1'])
                newobj.(props{i}) = obj.(props{i});
            end
            
            %Set id and fileRoot of new object to appropriate values
            newobj.id = idnum;
            str=['_ID',num2str(idnum),'_'];
            startind = regexp(newobj.fileRoot,'_ID[0-9]*_','start');
            endind   = startind + length(str) - 1;
            newobj.fileRoot(startind:endind) = ['_ID',num2str(idnum),'_'];
        end %Done
        
        function  [] = modifyICandTranslateBasisRk1(obj,newIC)
            
            switch obj.nBases
                case 1
                    if norm(newIC-obj.svdUpdateData.wRef) < 1e-4
                        fprintf('Skipping basis update...');
                        return;
                    end
                    
                    obj.prob.setProperty('ic',newIC);
                    obj.phi=ThinSVD_AddVec2Cols_Fast(obj.phi,obj.SingVals(1:obj.nY),...
                        obj.svdUpdateData.ColSumRSVtrunc,...
                        obj.svdUpdateData.normOneMinusVVt,...
                        obj.svdUpdateData.wRef-newIC);
                    obj.phi=obj.phi(:,1:end-1);
                    obj.phi=obj.phi(:,1:obj.trunc);
                    obj.svdUpdateData.wRef=newIC;
                    if length(obj.SingVals)==obj.nY, err=0; else err = 100*obj.SingVals(obj.nY+1)/obj.SingVals(1); end;
                    fprintf('Relative 2-norm error in update = %f%%\n',err);
                otherwise
                    for nn = 1:obj.nBases
                        if norm(newIC-obj.svdUpdateData(nn).wRef) < 1e-4
                            fprintf('Skipping basis update...');
                            return;
                        end
                        
                        obj.prob.setProperty('ic',newIC);
                        obj.phi0{nn}=ThinSVD_AddVec2Cols_Fast(obj.phi0{nn},obj.SingVals{nn}(1:obj.nY),...
                            obj.svdUpdateData(nn).ColSumRSVtrunc,...
                            obj.svdUpdateData(nn).normOneMinusVVt,...
                            obj.svdUpdateData(nn).wRef-newIC);
                        obj.phi0{nn}=obj.phi0{nn}(:,1:end-1);
                        
                        obj.svdUpdateData(nn).wRef=newIC;
                        if length(obj.SingVals{nn})==obj.nY, err=0; else err = 100*obj.SingVals{nn}(obj.nY+1)/obj.SingVals{nn}(1); end;
                        fprintf('Relative 2-norm error in update = %f%%\n',err);
                    end
            end
        end
        
        function  [] = addFOMsample(obj,sample)
            obj.fom_samples=[obj.fom_samples,sample];
        end
        
        function  [] = append2phi(obj,X)
            
            obj.phi = [obj.phi, X - obj.phi*(obj.phi'*X)];
            obj.phi(:,end) = obj.phi(:,end)/norm(obj.phi(:,end));
            %obj.phi=qr([obj.phi,X],0);
            obj.setProperty('nY',[],size(obj.phi,2));
            
        end
        
        %Optimization
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
            
            %Compute the residual sensitivities
            % [R,dRdw,dRdp]=obj.TimeScheme.TimeIntNLFuncSens(obj.prob,obj.sv(:,end),obj.sv(:,end),inf);
            if obj.killflag
                f = nan;
                DfDp = nan(size(z));
                return;
            end
            
            if strcmpi(class(obj.prob),'structuralFEM')
                [R,dRdw,dRdp] = obj.prob.ResJacPrecompSens(obj.sv(:,end),1,1);
            else
                [R,dRdw,dRdp] = obj.prob.ResSens(obj.sv(:,end));
            end
            
            %Evaluate objective and sensitivities
            [f,dfdw,dfdp] = obj.objectiveNAND(obj.sv(:,end),z);
            
            if strcmpi(sensMethod,'findiff')
                DfDp=[];
                return;
            end
            
            switch lower(obj.Rtype)
                case 'petrov-galerkin'
                    [d2Rdu_r,d2Rdudz_r] = SecondDerivsPG(obj,R,dRdw,dRdp,obj.phi,1e-5,inf,obj.sv(:,end),2);
                    
                    Jr = dRdw*obj.phi;
                    dRrdy = Jr'*Jr + d2Rdu_r;
                    dRrdp = Jr'*dRdp + d2Rdudz_r;
                    
                    %                     if ~isempty(d2RdwdpR)
                    %                         [~,Rj] = qr(dRdw,0);
                    %
                    %                         dRrdy = obj.phi'*(d2Rdw2R + Rj'*Rj)*obj.phi;
                    %                         dRrdp = obj.phi'*(d2RdwdpR + dRdw'*dRdp);
                    %                     else
                    %                         dRrdy = obj.phi'*(dRdw*obj.phi);
                    %                         dRrdp = obj.phi'*dRdp;
                    %                     end
                    
                case 'galerkin'
                    
                    if strcmpi(class(obj.prob),'structuralFEM')
                        dRrdy = dRdw; dRrdp = dRdp;
                    else
                        dRrdy = obj.phi'*(dRdw*obj.phi);
                        dRrdp = obj.phi'*dRdp;
                    end
            end
            
            if strcmpi(class(obj.prob),'structuralFEM')
                dfdy=dfdw;
            else
                dfdy=obj.phi'*dfdw;
            end
            
            switch sensMethod
                case 'direct'
                    %Using the direct method, solve for the state
                    %sensitivity (w.r.t. parameter) and compute the reduced
                    %gradient
                    DwDp = -dRrdy\dRrdp;
                    DfDp = dfdp + DwDp'*dfdy;
                    %DwDp = -dRdw\dRdp;
                    %DfDp = dfdp + DwDp'*dfdw;
                    lambda = [];
                case 'adjoint'
                    %Using the adjoint method, solve for the dual variable
                    %and compute the reduced gradient
                    lambda = -(dRrdy'\dfdy);
                    DfDp = dfdp + dRrdp'*lambda;
                    %lambda = -dRdw'\dfdw;
                    %DfDp = dfdp + dRdp'*lambda;
            end
            
            %             switch lower(obj.Rtype)
            %                 case 'petrov-galerkin'
            %                     %Compute the residual sensitivities
            %                     [dRdw,dRdp,d2RdwdpR,d2Rdw2R] = obj.prob.ResSens(obj.sv(:,2));
            %                     %Evaluate objective and sensitivities
            %                     [f,dfdw,dfdp] = obj.objective(obj.sv(:,2),z);
            %
            %                     [~,Rj] = qr(dRdw,0);
            %
            %                     dRrdy = obj.phi'*(d2Rdw2R + Rj'*Rj)*obj.phi;
            %                     dRrdp = obj.phi'*(d2RdwdpR + dRdw'*dRdp);
            %                     switch sensMethod
            %                         case 'direct'
            %                             %Using the direct method, solve for the state
            %                             %sensitivity (w.r.t. parameter) and compute the reduced
            %                             %gradient
            %                             DwDp = -dRrdy\dRrdp;
            %                             DfDp = dfdp + DwDp'*(obj.phi'*dfdw);
            %                             %DwDp = -dRdw\dRdp;
            %                             %DfDp = dfdp + DwDp'*dfdw;
            %                             lambda = [];
            %                         case 'adjoint'
            %                             %Using the adjoint method, solve for the dual variable
            %                             %and compute the reduced gradient
            %                             lambda = -(dRrdy'\(obj.phi'*dfdw));
            %                             DfDp = dfdp + dRrdp'*lambda;
            %                             %lambda = -dRdw'\dfdw;
            %                             %DfDp = dfdp + dRdp'*lambda;
            %                     end
            %                 case 'galerkin'
            %                     %Compute the residual sensitivities
            %                     [dRdw,dRdp] = obj.prob.ResSens(obj.sv(:,2));
            %                     %Evaluate objective and sensitivities
            %                     [f,dfdw,dfdp] = obj.objective(obj.sv(:,2),z);
            %
            %                     switch sensMethod
            %                         case 'direct'
            %                             %Using the direct method, solve for the state
            %                             %sensitivity (w.r.t. parameter) and compute the reduced
            %                             %gradient
            %                             DwDp = -(obj.phi'*dRdw*obj.phi)\(obj.phi'*dRdp);
            %                             DfDp = dfdp + (obj.phi*DwDp)'*dfdw;
            %                             %DwDp = -dRdw\dRdp;
            %                             %DfDp = dfdp + DwDp'*dfdw;
            %                             lambda = [];
            %                         case 'adjoint'
            %                             %Using the adjoint method, solve for the dual variable
            %                             %and compute the reduced gradient
            %                             lambda = -(obj.phi'*dRdw*obj.phi)'\(obj.phi'*dfdw);
            %                             DfDp = dfdp + (dRdp'*obj.phi)*lambda;
            %                             %lambda = -dRdw'\dfdw;
            %                             %DfDp = dfdp + dRdp'*lambda;
            %                     end
            %             end
        end
        
        function   [f,df] = ObjFMINCON_SAND(obj,u)
            
            w=u(1:obj.nY,1);
            p=u(obj.nY+1:end,1);
            
            obj.prob.updateParameters(p);
            
            [f,dfdw,dfdp] = obj.objectiveSAND(w,p);
            df = [obj.phi'*dfdw;dfdp];
            
        end
        
        function   [f,df] = ObjFMINCON_SAND_ROMVAR(obj,u)
            
            w=u(1:obj.nY,1);
            p=u(obj.nY+1:end,1);
            
            obj.prob.updateParameters(p);
            
            [f,dfdw,dfdp] = obj.objectiveSAND_ROMvar(w,p);
            
            df = [dfdw;dfdp];
        end
        
        function   [f,df] = ObjFMINCON_NAND(obj,p,sensMethod,flag)
            
            if norm(p - obj.curr_param) > 0
                switch class(obj.prob)
                    case {'SteadyNozzle','quasi1dEuler'}
                        obj.prob.updateParameters(p);
                        
                        R(1)=norm(obj.prob.ResJac(obj.svdUpdateData.wRef));
                        R(2)=inf;%norm(obj.prob.ResJac(obj.sv(:,end)));
                        for i = 1:size(obj.fom_samples,2)
                            R(2+i)=norm(obj.prob.ResJac(obj.fom_samples(:,i)));
                        end
                        [~,min_ind] = min(R);
                        
                        if (min_ind==2), obj.modifyICandTranslateBasisRk1(obj.sv(:,end)); end;
                        if (min_ind>2) , obj.modifyICandTranslateBasisRk1(obj.fom_samples(:,min_ind-2)); end;
                        %if R >=R0, obj.modifyICandTranslateBasisRk1(obj.svdUpdateData.wRef); end;
                        
                        if (min_ind==2), obj.modifyICandTranslateBasisRk1(obj.sv(:,end)); end;
                        if (min_ind>2) , obj.modifyICandTranslateBasisRk1(obj.fom_samples(:,min_ind-2)); end;
                        %if R >=R0, obj.modifyICandTranslateBasisRk1(obj.svdUpdateData.wRef); end;
                    case 'structuralFEM'
                        n1 = obj.prob.PrecompROM.n_lam; n2=obj.prob.PrecompROM.n_mu; n3=obj.prob.PrecompROM.n_rho;
                        obj.prob.PrecompFromRedMatPrecomp(p(1:n1),p(n1+1:n1+n2),p(n1+n2+1:end));
                end
                obj.executeModel();
            end
            if nargin < 4 || isempty(flag)
                [df,f] = obj.reducedGrad(p,sensMethod);
            elseif strcmpi(flag,'objective')
                f = obj.objectiveNAND(obj.sv(:,end),p);
                %[~,f] = obj.reducedGrad(p,sensMethod);
                if isnan(f)
                    f=1e6;
                end;
                df=[];
            elseif strcmpi(flag,'gradient')
                [f,~] = obj.reducedGrad(p,sensMethod);
                if sum(isnan(f))>0
                    disp('NAN in GRADIENT!  Changing to zeros!');
                    f=zeros(size(f));
                end
                df=[];
            end
            %             [df,f] = obj.reducedGrad(p,sensMethod);
        end
        
        function   [cineq,ceq,dcineq,dceq] = NLConstFMINCON_SAND(obj,u,varargin)
            
            w=u(1:obj.nY,1);
            p=u(obj.nY+1:end,1);
            obj.prob.updateParameters(p);
            
            if strcmpi(obj.Rtype,'Galerkin')
                [R,dRdw,dRdp] = obj.prob.ResSens(obj.prob.ic + obj.phi*w,[]);
                
                ceq = obj.phi'*R;
                dceq = [obj.phi'*dRdw*obj.phi, obj.phi'*dRdp]';
                
                cineq=[];
                dcineq=[];
            else
                [R,dRdw,dRdp] = obj.prob.ResSens(obj.prob.ic + obj.phi*w,[]);
                %                 wF = obj.prob.ic + obj.phi*w;
                %                 [R,dRdw] = obj.prob.ResSens(obj.prob,wF,wF,inf);
                
                ceq = obj.phi'*dRdw'*R;
                
                [d2Rdu_r,d2Rdudz_r] = SecondDerivsPG(obj,R,dRdw,dRdp,obj.phi,1e-5,inf,2);
                Jr = dRdw*obj.phi;
                dRrdy = Jr'*Jr + d2Rdu_r;
                dRrdp = Jr'*dRdp + d2Rdudz_r;
                dceq = [dRrdy,dRrdp]';
                
                %dRdwPhi = dRdw*obj.phi;
                %dceq = [(dRdwPhi'*dRdwPhi), dRdwPhi'*dRdp]';
                %dceq = [];
                
                cineq=[];
                dcineq=[];
            end
            
            %Include additional constraints
            for i = 1:length(varargin)
                [tmp,tmpeq,dtmp,dtmpeq] = varargin{i}(w,p);
                
                cineq=[cineq;tmp];
                ceq=[ceq;tmpeq];
                
                if size(dtmp,1) == length(p)
                    dtmp = [zeros(obj.nY,size(dtmp,2));dtmp];
                end
                
                if size(dtmpeq,1) == length(p)
                    dtmpeq = [zeros(obj.nY,size(dtmpeq,2));dtmpeq];
                end
                
                dcineq=[dcineq,dtmp];
                dceq=[dceq,dtmpeq];
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
            
            obj.objectiveSAND = eval(['@(u,z) obj.prob.',txt,'(obj.prob.ic + obj.phi*u,z)']);
            obj.objectiveSAND_ROMvar = eval(['@(u,z) obj.prob.',txt,'(u,z)']);
            obj.objectiveNAND = eval(['@(u,z) obj.prob.',txt,'(u,z)']);
        end
        
        function  [] = setOptConstraints(obj,txt)
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
            
            obj.constraints = eval(['@(u,z) obj.prob.',txt,'(u,z)']);
        end
        
        function  [A,c] = reducedConstJac(obj,z,sensMethod)
            %This function computes the reduced Jacobian of the constraints
            %using either the direct or adjoint method to compute the
            %sensitivity of the state with respect to the parameter.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj    - ROM object
            %z      - parameter at the current parameter iterate
            %Outputs:
            %--------
            %A      - Jacobian of constraint vector
            %c      - value of constraint vector
            %--------------------------------------------------------------
            
            %INEFFICIENT - REDO THIS BY TAKING ADVANTAGE OF THE FACT THAT
            %dudw was previously computed
            
            %Compute the residual sensitivities
            [dRdw,dRdp] = obj.prob.ResSens(obj.sv(:,2));
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
            %ROM object
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - ROM object
            %
            %Outputs:
            %--------
            %svF     - full state vector matrix
            %--------------------------------------------------------------
            
            svF = obj.sv;
        end %Done
        
        %Unused?
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
            
            all2clear4PP = {'res','jac','clusCenter','UrLoc'};
            
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
        
        function  [] = copyOfflineResults(obj,romobj)
            %This function sets the phi property of the current object to
            %ROBphi.  This will be useful in the event that have several
            %ROMs with the same ROB (i.e. you can compute the ROB once and
            %then just set the ROB of the remaining ROMs to this computed
            %ROB).
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj    - ROM object
            %romobj - ROM object whose offline computations will be copied
            %         into obj (instead of recomputing).
            %
            %Outputs:
            %--------
            %There are no outputs, obj.phi is set.
            %--------------------------------------------------------------
            
            %Error checking on the dimension of phi
            if size(romobj.sv,1) ~= size(obj.sv,1)
                error('Must pass in a ROM object whose FOM has the same number of DOFs as the FOM of the current object');
            end
            if romobj.nBases ~= obj.nBases
                error('Must pass in a ROM object that has the same number of bases as the current object');
            end
            
            obj.nY = romobj.nY;
            obj.phi = romobj.phi;
            if obj.nBases > 1
                obj.clusCenter = romobj.clusCenter;
                obj.UrLoc = romobj.UrLoc;
                obj.f = romobj.f;
                obj.S = romobj.S;
                %Do not want to copy res and jac from some romobj because
                %they might not
                obj.res = cell(obj.nBases,1);
                obj.jac = cell(obj.nBases,1);
            end
        end
        
        function  [] = setProperty(obj,prop1,prop2,val)
            %Should include some checking to:
            %1) make sure prop is a string corresponding to a property name
            %2) val is valid for the specific property
            %3) the user isn't trying to change a property that shouldn't
            %be changed
            if isempty(prop2)
                obj.(prop1)=val;
            else
                obj.(prop1).(prop2)=val;
            end
        end
        
        function  [] = setFieldsFromOldROM(obj,oldROM,varargin)
            
            % Construct a new object based on a deep copy of an old object
            % of this class by copying properties over.
            if isa(oldROM,'ROM')
                props = properties(oldROM);
                [~,ind] = intersect(props,varargin);
                for i = ind(:)'
                    % Use Dynamic Expressions to copy the required property.
                    % For more info on usage of Dynamic Expressions, refer to
                    % the section "Creating Field Names Dynamically" in:
                    % web([docroot '/techdoc/matlab_prog/br04bw6-38.html#br1v5a9-1'])
                    obj.(props{i}) = oldROM.(props{i});
                end
                return;
            end
            
        end
        
        function  [snapMat] = determineSnapshots(obj,fomobj,trainIC)
            %This function determines the state vector snapshot matrix from
            %FOM object input while adhereing to the user requests
            %indicated in the snaps structure
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - ROM object
            %fomobj  - FOM object to extract state vectors from
            %trainIC - ndof x 1 vector defining the training initial
            %          condition.  This is to be used if multiple
            %          simulations comprise the snapshots for a Local POD
            %          simulation.  If this is empty, the initial condition
            %          of the fomobj is used.
            %
            %Outputs:
            %--------
            %snapMat - snapshot matrix
            %--------------------------------------------------------------
            
            %Handle the steady state case.  Only use the 2nd column of the
            %state vector matrix (the first column is the initial guess).
            if fomobj.prob.staticFlag
                %snapMat = fomobj.sv(:,2);
                if fomobj.saveAllIt
                    if ~isempty(trainIC)
                        snapMat = snapshotCollection(fomobj.svIter,obj.snaps.initref,trainIC);
                    else
                        snapMat = snapshotCollection(fomobj.svIter,obj.snaps.initref,fomobj.sv(:,1));
                    end
                else
                    if ~isempty(trainIC)
                        snapMat = snapshotCollection(fomobj.sv,obj.snaps.initref,trainIC);
                    else
                        snapMat = snapshotCollection(fomobj.sv,obj.snaps.initref,fomobj.sv(:,1));
                    end
                end
                return;
            end
            
            %Determine the all possible snapshots from the state vector
            %matrix without considering the snapshot distribution,
            %interval, etc.
            if ~isempty(trainIC)
                Fsv = snapshotCollection(fomobj.sv,obj.snaps.initref,trainIC);
            else
                Fsv = snapshotCollection(fomobj.sv,obj.snaps.initref,fomobj.sv(:,1));
            end
            %Determine the time vector from the fomobj time structure
            %and remove the first time (corresponds to ic)
            t = linspace(fomobj.time.T(1),fomobj.time.T(2),fomobj.time.nstep+1);
            t(:,1) = [];
            
            %Determine the indicies of t that are valid for use as
            %snapshots (i.e. in the snapshot interval) and the maximum
            %number of snapshots in the interval specified.
            if ~isempty(obj.snaps.int)
                indValid = find((t >= obj.snaps.int(1,1)) & (t <= obj.snaps.int(1,2)));
                nsnapMAX = numel(indValid);
            else
                %if obj.snaps.int is empty, the entire time domain is used
                indValid = 1:fomobj.time.nstep;
                nsnapMAX = fomobj.time.nstep;
            end
            
            if isempty(obj.snaps.num) || (ischar(obj.snaps.num) && strcmpi(obj.snaps.num,'all'))
                %If obj.snaps.num is a char array whose value is 'all',
                %use all snapshots in the indicated interval
                numsnaps = nsnapMAX;
            elseif ischar(obj.snaps.num) && ~strcmpi(obj.snaps.num,'all')
                %If obj.snaps.num is a char array whose value is not 'all',
                %return an error
                error('Only valid keyword in the nsnap property is all');
            else
                %Otherwise, set the number of snapshots to the number
                %indicated in the input file (as long as obj.snaps.num
                %< nstep)
                if obj.snaps.num > fomobj.time.nstep
                    fprintf('Warning:  Too many snapshots specified.  Using nsnapMAX = %d snapshots',nsnapMAX);
                    numsnaps = nsnapMAX;
                else
                    numsnaps = obj.snaps.num;
                end
            end
            
            if numsnaps == nsnapMAX
                %If all snapshots in the interval are to be used, simply
                %return those indices.
                ind = indValid;
            else
                %Otherwise, compute them from the specified
                %distribution.
                if ~isempty(obj.snaps.dist) && ~isempty(obj.snaps.distseed) && ~isempty(obj.snaps.distparams)
                    ind = computeSnapshotFromDist(numsnaps,obj.snaps.distparams,...
                        obj.snaps.dist,obj.snaps.distseed,indValid(1),indValid(end));
                else
                    ind = computeSnapshotFromDist(numsnaps,indValid(end),'uniform',1,[],[]);
                end
            end
            %Return the appropriate snapshots
            snapMat = Fsv(:,ind);
        end %Replaced with extractRawSnapshots?
    end
end

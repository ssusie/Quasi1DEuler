classdef locGNAT < handle
    
    properties (SetAccess = private, GetAccess = public)
        %Quantities read from input file
        id;
        prob;
        probGNAT;
        nY;
        nR;
        nJ;
        nI;
        newt;
        time;
        nSample; %number of sample nodes to use
        nSampleProd;
        nGreed; %number of basis vectors to use in greed node selection
        nBases;
        addInd;
        
        %Compute quantities
        TimeScheme;
        ontime;
        LocBasisHist;
        
        curr_param;
        
        S; %Local POD precomputed
        f; %Local POD precomputed
        Q;
        small_mask;
        smallMaskIndexMatchForQ;
        localMaskIndexMatchForQ;
        reconstructedLocalROB;
        phi0tmp;
        wReftmp;
        icTmp;
        
        basisUpdate;
        initref;
        svdUpdateData;
        wr_switch;
        SingVals;
        ROBcompon;
        numswitch;
        switchHist;
        prevSwitchState;
        pprevBasisNum;
        
        preCompDistQuant
        
        sv; %State vector (reduced) nY x (nstep+1) matrix
        A; %Online matrix A from Carlberg et. al. 2011
        B; %Online matrix B from Carlberg et. al. 2011
        sampleNodes; %Vector containing the sample node numbers
        sampleInd; %Vector containing the sample indice numbers
        %(will be the same as sampleNodes if there is only 1 unknown per node)
        
        irstart; %vector of indices pointing to locations in jrow to indicate the start of a new indice
        jrow; %vector of indices indicating the indices of the full state vector where the state vector needs to be evaluated
        jdiag; %boolean vector the same size as jrow indicating whether or not that component corresponds to a diagonal entry
        phiYhat; %the partial reduced order basis (i.e. phiY(jrow,:))
        phiYhat0;
        cTimeIter; %current time iteration number
        
        numExecute=0;
    end
    
    properties (Hidden=true, SetAccess = private, GetAccess = public)
        %Temporary properties
        partialULoc; %partial state vector at the current time step
        partialUprevLoc; %partial state vector at the previous time step
        JhatTemp; %matrix of size nI x number of necessary state vectors used to store the jacobian
        %Phi_hat from Carlberg et. al 2010 (online algorithm)
        reconstJhatInd; %vector that is used take the partial Jacobian from vector form into matrix form (JhatTemp),
        %i.e. JhatTemp(reconstJhatInd) = Jhat will take the vector Jhat into the appropriate matrix JhatTemp
        uniqueJROWind;
        jdiagHat;
        cLocBasis; %Local basis currently being used
        Ur_Loc;
        phiYhatUnion; %This is the unique union of the rows of all phiYhat
        indOfSvPerBasis;
        GVARtxt;
        
        killflag;
    end
    
    methods
        %Constructor
        function [obj] = locGNAT(ROMfile,romobj,id,oldobj)
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
            if nargin == 4 && isa(oldobj,'locGNAT')
                copyProperties(obj,oldobj);
                return;
            end
            
            %Copy the time/newt structure and store in the class instance
            %The default time structure is whatever is stored in ROM
            obj.time = romobj.time;
            obj.newt = struct('maxIter',[],'eps',[],'iter',[],'quiet',[]);
            
            %Copy other properties from the ROM
            obj.nY = romobj.nY;
            obj.nBases = romobj.nBases;
            %             obj.S = romobj.S;
            %             obj.f = romobj.f;
            %             obj.phi = romobj.phi;
            
            obj.svdUpdateData = romobj.svdUpdateData;
            %obj.preCompDistQuant = romobj.preCompDistQuant;
            obj.numswitch=0;
            obj.basisUpdate = romobj.basisUpdate;
            obj.initref = romobj.snaps.initref;
            
            %Extract parameters from ROM file
            obj.GVARtxt = romobj.GVARtxt;
            VARtext = readInFile('VAR',ROMfile,1);
            %Extract the GNAT text from the rom file
            GNATtext = readInFile('GNAT',ROMfile,1);
            
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
            if size(obj.nR,1) == 1
                obj.nR = obj.nR*ones(obj.nBases,1);
            elseif ~isempty(obj.nR) && size(obj.nR,1) ~= obj.nBases
                error('Size of nR and nBases not compatible');
            end
            obj.nJ = extractInputROM(GNATtext,N,obj.id,'nJ',[]); obj.nJ = obj.nJ(:);
            if size(obj.nJ,1) == 1
                obj.nJ = obj.nJ*ones(obj.nBases,1);
            elseif ~isempty(obj.nJ) && size(obj.nJ,1) ~= obj.nBases
                error('Size of nJ and nBases not compatible');
            end
            obj.nI = extractInputROM(GNATtext,N,obj.id,'nI',[]); obj.nI = obj.nI(:);
            if size(obj.nI,1) == 1
                obj.nI = obj.nI*ones(obj.nBases,1);
            elseif ~isempty(obj.nI) && size(obj.nI,1) ~= obj.nBases
                error('Size of nI and nBases not compatible');
            end
            obj.nSample = extractInputROM(GNATtext,N,obj.id,'nSample',[]); obj.nSample = obj.nSample(:);
            if size(obj.nSample,1) == 1
                obj.nSample = obj.nSample*ones(obj.nBases,1);
            elseif ~isempty(obj.nSample) &&  size(obj.nSample,1) ~= obj.nBases
                error('Size of nSample and nBases not compatible');
            end
            obj.nGreed  = extractInputROM(GNATtext,N,obj.id,'nGreed',[]); obj.nGreed = obj.nGreed(:);
            if size(obj.nGreed,1) == 1
                obj.nGreed = obj.nGreed*ones(obj.nBases,1);
            elseif ~isempty(obj.nGreed) && size(obj.nGreed,1) ~= obj.nBases
                error('Size of nGreed and nBases not compatible');
            end
            obj.nSampleProd = extractInputROM(GNATtext,N,obj.id,'nSampleProd',[]);
            
            obj.sampleNodes = cell(obj.nBases,1);
            obj.sampleInd = cell(obj.nBases,1);
            obj.Ur_Loc = cell(obj.nBases,1);
            obj.A = cell(obj.nBases,1);
            obj.B = cell(obj.nBases,1);
            obj.reconstJhatInd = cell(obj.nBases,1);
            obj.JhatTemp = cell(obj.nBases,1);
            obj.jdiag = cell(obj.nBases,1);
            obj.jrow = cell(obj.nBases,1);
            obj.irstart = cell(obj.nBases,1);
            
            %%%%%%%%% Determine the 'free' indices to add %%%%%%%%%
            obj.addInd = extractInputROM(GNATtext,N,obj.id,'addInd','none');
            if ischar(obj.addInd)
                obj.addInd = {obj.addInd};
            end
            if size(obj.addInd,1) == 1
                obj.addInd = repmat(obj.addInd,obj.nBases,1);
            elseif size(obj.addInd,1) ~= obj.nBases
                error('addInd is not compatible with nBases');
            end
            %%%%%%%%% Determine time structure %%%%%%%%%
            obj.time.T     = extractInputROM(GNATtext,N,obj.id,'T',obj.time.T);
            obj.time.dt    = extractInputROM(GNATtext,N,obj.id,'dt',obj.time.dt);
            obj.time.nstep = extractInputROM(GNATtext,N,obj.id,'nstep',obj.time.nstep);
            obj.time.quiet = extractInputROM(GNATtext,N,obj.id,'timeQuiet',obj.time.quiet);
            %disp('Warning: put the next 3 lines back on in localGNAT !!');
            obj.time.steadyconverge = extractInputROM(GNATtext,N,obj.id,'steadyconverge',obj.time.steadyconverge);
            obj.time.cfl = extractInputRobust(GNATtext,'cfl',obj.time.cfl);
            
        %    if isempty(obj.time.dt) && isempty(obj.time.cfl)
            if isempty(obj.time.dt)
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
        function  [] = createGNAT(obj,probobj,romobj,phiR,phiJ,basisNum)
            %This function performs all of the offline GNAT computations.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - instance of GNAT class
            %probobj - Problem object
            %phiR    - nBases x 1 cell array of local bases of subspace
            %          residual is constrained to lie in
            %phiJ    - nBases x 1 cell array of local bases of subspace
            %          jacobian is constrained to lie in
            %phiY    - nBases x 1 cell array of local bases of subspace
            %          model is constrained to lie in (reduced order basis)
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
                computeSampleIndices(obj,probobj,phiR,phiJ,basisNum);
                
                %Compute additional nodes required in the reduced mesh in order
                %to compute the residual and jacobian at the sample nodes.
                determineAdditionalNodes(obj,probobj,basisNum);
            elseif ~isempty(obj.nSample) && ~isempty(obj.nGreed)
                %Compute the sample nodes based on a greedy node selection
                %algorithm
                computeSampleNodes(obj,probobj,phiR,phiJ,basisNum);
                
                %Compute additional nodes required in the reduced mesh in order
                %to compute the residual and jacobian at the sample nodes.
                determineAdditionalNodes(obj,probobj,basisNum);
                obj.nI(basisNum) = length(obj.sampleInd{basisNum});
            end
            
            %Mask SVD update data, if border_fast
            if strcmpi(romobj.basisUpdate,'border_fast') || strcmpi(romobj.basisUpdate,'border_fast_approx')
                obj.phi0tmp=romobj.phi0;
                obj.wReftmp{basisNum}=romobj.svdUpdateData(basisNum).wRef;
                obj.icTmp = romobj.prob.ic;
                
                obj.SingVals{basisNum} = romobj.SingVals{basisNum};
                obj.svdUpdateData(basisNum).wRef=zeros(length(unique(obj.jrow{basisNum,1})),obj.nBases);
                for kk = 1:obj.nBases
                    %For wRef_i with mask j, use obj.svdUpdateData(j).wRef(:,i)
                    obj.svdUpdateData(basisNum).wRef(:,kk) = romobj.svdUpdateData(kk).wRef(unique(obj.jrow{basisNum,1}),1);
                end
            end
            
            for kk = 1:obj.nBases
                %Mask jj, basis kk
                if strcmpi(romobj.basisUpdate,'border_fast')% || strcmpi(romobj.basisUpdate,'border_exact')
                    obj.phiYhat0{basisNum,kk} = romobj.phi0{kk,1}(unique(obj.jrow{basisNum,1}),:);
                else
                    obj.phiYhat0{basisNum,kk} = romobj.phi0{kk,1}(unique(obj.jrow{basisNum,1}),:);
                end
            end
            
            %Determine nI based on the number of sample indices
            obj.nI(basisNum) = length(obj.sampleInd{basisNum,1});
            
            %Check to ensure that theoretical restrictions on nR, nJ, nI,
            %nY are met.
            if (obj.nR(basisNum) < obj.nY(basisNum)) || (obj.nJ(basisNum) < obj.nY(basisNum))
                error(['To ensure uniqueness in the reconstructed gappy quantities,',...
                    'nR >= nY and nJ >= nY are required ']);
            end
            
            if (obj.nI(basisNum)~=0) && ((obj.nI(basisNum) < obj.nR(basisNum)) || (obj.nI(basisNum) < obj.nJ(basisNum)))
                error(['To ensure uniqueness in the least squares gappy reconstruction,',...
                    'nI >= nR and nI >= nJ are required ']);
            end
            
            %Compute the online matrices A and B
            if ~isempty(phiR) && ~isempty(phiJ)
                phiRtrunc = phiR(:,1:obj.nR(basisNum)); clear phiR;
                phiJtrunc = phiJ(:,1:obj.nJ(basisNum)); clear phiJ;
                computeOnlineMatrices(obj,phiRtrunc,phiJtrunc,basisNum);
            end
            
            %Initialize the state vector matrix (reduced coordinates)
            obj.sv = zeros(max(obj.nY),obj.time.nstep+1); %Reduced ic = 0
            
            %Setup Ur_Loc for local GNAT computations
            obj.Ur_Loc{basisNum,1} = zeros(obj.nY(basisNum),1);
            
            %Determine phiYhatInd.  Let X be a full state vector, then
            %Xpartial = X(phiYhatInd,:), where Xpartial has no repeats.
            obj.reconstJhatInd{basisNum,1} = [];
            [~,~,J] = unique(obj.jrow{basisNum,1});
            for i = 1:length(obj.irstart{basisNum,1})-1
                obj.reconstJhatInd{basisNum,1} = [obj.reconstJhatInd{basisNum,1}, i + obj.nI(basisNum)*(J(obj.irstart{basisNum,1}(i):obj.irstart{basisNum,1}(i+1)-1)'-1)];%):J(obj.irstart(i+1)-1)]-1)];
            end
            obj.uniqueJROWind{basisNum,1} = J;
            
            %Initialize (sparse) JhatTemp
            obj.JhatTemp{basisNum,1} = spalloc(obj.nI(basisNum),length(unique(obj.jrow{basisNum,1})),length(obj.reconstJhatInd{basisNum,1}));
            
            %Set up the vector containing the indices of the partial
            %vectors that correspond to diagonal entries
            obj.jdiagHat{basisNum,1} = obj.uniqueJROWind{basisNum,1}(logical(obj.jdiag{basisNum,1}));
        
            %Store handle to problem in instance
            obj.probGNAT = probobj.createCopy4GNATloc(obj);
            
            % delete prob properties
            obj.prob = [];
            
            %Add GNAT parameters to time scheme
            obj.TimeScheme.addGNAT(obj);
        end
        
%         function  [] = reduceProbObj(obj)
%             
%             %Store handle to problem in instance
%             obj.probGNAT = probobj.createCopy4GNATloc(obj);
%             %Add GNAT parameters to time scheme
%             obj.TimeScheme.addGNAT(obj);
%             
%         end
        
        function  [] = executeModel(obj)
            %This function performs the time-stepping for the GNAT
            %simulation defined in this object.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - GNAT object
            %
            %Outputs:
            %--------
            %This function has no outputs.
            %--------------------------------------------------------------
            obj.numExecute = obj.numExecute+1;
            obj.killflag = false;
            
            nstep = obj.time.nstep; %Extract the number of time steps from the object
            
            obj.LocBasisHist = zeros(nstep,1);
            obj.switchHist   = zeros(nstep,1);
            
            obj.numswitch = 0;
            
            
            %Ensure sv = 0 initially!
            obj.sv = zeros(max(obj.nY),obj.time.nstep+1);
            
            %Set the current parameter (from prob object) to be used in
            %NAND to know which parameter values have been used to generate
            %the current state vectors
            obj.curr_param = obj.probGNAT(1).p;
            
            tLocGNATstart = tic;
            for i_t = 1:nstep %Loop over each time step to determine all state vectors
                % Compute reference state distances to bases centers
                obj.cLocBasis = computeDistanceToCenters(obj,obj.Ur_Loc,i_t);
                obj.LocBasisHist(i_t,1) = obj.cLocBasis;
                
                if i_t == 1
                    %Make sure partialU and partialUprev are properly set whenever this
                    %function is called
                    obj.setUpartial2ic;
                end
                
                if ~obj.time.quiet && rem(i_t,round(nstep*0.1)) == 0
                    %Generate output so user can see progress
                    fprintf('%4.0f of %4.0f Timesteps Complete (%2.0f%%)\n',i_t,nstep,100*i_t/nstep)
                end
                
                %Set current time iteration (so model can keep track)
                obj.cTimeIter = i_t;
                
                %Determine initial ROB components if using fast
                %updating
                if obj.cTimeIter == 1 && strcmpi(obj.basisUpdate,'border_fast')
                    obj.ROBcompon(1).alpha = zeros(obj.nY(obj.cLocBasis),1);
                    for i = 1:obj.nBases
                        obj.ROBcompon(1).beta{i}  = zeros(obj.nY(obj.cLocBasis),1);
                        if i == obj.cLocBasis
                            obj.ROBcompon(1).N{i} = eye(obj.nY(i));
                        else
                            obj.ROBcompon(1).N{i} = zeros(obj.nY(i),obj.nY(obj.cLocBasis));
                        end
                    end
                end
                
                %Use Newton's Method to solve nonlinear system of equations
                %and store: state vector, residual and jacobian snapshots,
                %and number of newton iterations required

                obj.NewtonRaphsonLocal;
                if obj.killflag, warning('Encountered NaN at Time Step %i. Simulation Killed',obj.cTimeIter); return; end;
                
                %Check for steady convergence
            %    if (i_t==1)
             %      disp('Warning: put the next  line back on in localGNAT !!');
              %  end
               % if 1==0
                if obj.time.steadyconverge
                    if norm(obj.sv(:,i_t+1)-obj.sv(:,i_t)) < obj.time.steadyconverge*norm(obj.sv(:,i_t+1))
                        obj.ontime = toc(tLocGNATstart); %Record simulation time
                        obj.time.nstep=i_t;
                        obj.time.T(2) = obj.time.T(1) + i_t*obj.time.dt;
                        obj.sv(:,i_t+2:end)=[];
                        obj.newt.avgIter = mean(obj.newt.iter);
                        return;
                    end
                end
            end
                            obj.time
                obj.TimeScheme
            obj.ontime = toc(tLocGNATstart); %Record simulation time
            
            obj.newt.avgIter = mean(obj.newt.iter);
        end
        
        function  [] = associateFullProblem(obj,probobj)
            %This function associates the full problem (i.e. not the sample
            %mesh) with the locGNAT object (i.e. it is stored in the prob
            %property)
            %---------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - locGNAT object
            %probobj - Problem object
            %
            %Outputs:
            %--------
            %There are no outputs.
            %---------------------------------------------------------------
            
            obj.prob = probobj;
        end
        
        function  [svF] = reconstructFullState(obj,romobj,finalState)
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
            
            if nargin < 3 || isempty(finalState)
                finalState = obj.time.nstep;
            end
            
            svF = zeros(size(romobj.phi0{1},1),finalState+1);
            %svF = zeros(size(romobj.phi,1),romobj.time.nstep+1);
            %First column is the initial condition.
            svF(:,1) = obj.prob.ic;
            for i = 2:finalState+1
                %The newton method that was implemented for the GNAT
                %computes the current reduced state vector relative to the
                %previous state vector.  Therefore, we reconstructed the
                %full state by looping through the time steps and adding
                %the reconstructed state difference to the previous state
                %vector.
                dw = obj.sv(1:size(romobj.phi0{obj.LocBasisHist(i-1),1},2),i);
                if strcmpi(obj.basisUpdate,'border_fast')
                    svF(:,i) = svF(:,i-1) + obj.computeFullPhiTimesVecFromComp(romobj,obj.switchHist(i-1,1),dw);
                elseif strcmpi(obj.basisUpdate,'border_fast_approx')
                    if (i==2)
                        switchNum = 1;
                        reconstructLocalROB(obj,romobj,[],switchNum,obj.LocBasisHist(1));
                    elseif (i>2 && obj.LocBasisHist(i-1)~=obj.LocBasisHist(i-2))
                        switchNum = switchNum + 1;
                        reconstructLocalROB(obj,romobj,svF(:,i-1),switchNum,obj.LocBasisHist(i-1));
                    end
                    svF(:,i) = svF(:,i-1) + obj.reconstructedLocalROB{switchNum}*dw;
                elseif strcmpi(obj.basisUpdate,'none')
                    svF(:,i) = svF(:,i-1) + romobj.phi0{obj.LocBasisHist(i-1),1}*dw;
                end
            end
        end
        
        function [] = modifyTimeProperties(obj,fomobj)
            
            obj.time = fomobj.time;
            obj.prob.config.setTimeObj('T',fomobj.prob.config.time.T);
            obj.prob.config.setTimeObj('dt',fomobj.prob.config.time.dt);
            obj.prob.config.setTimeObj('nstep',fomobj.prob.config.time.nstep);
            obj.prob.config.setTimeObj('quiet',fomobj.prob.config.time.quiet);
            obj.TimeScheme = fomobj.TimeScheme;
            obj.TimeScheme.addGNAT(obj);

        end
        
        function [] = modifyInitialCondition(obj,romobj,probobj)
            
            obj.initref = 'somevec';
            [obj.probGNAT] =  probobj.updateIcCopy4GNATloc(obj,obj.probGNAT);
            for iLocBasis = 1:obj.nBases
                obj.svdUpdateData(iLocBasis).a = norm(probobj.ic,2);
                obj.svdUpdateData(iLocBasis).c = probobj.ic'*romobj.svdUpdateData(iLocBasis).wRef;
                obj.svdUpdateData(iLocBasis).d = romobj.phi0{iLocBasis,1}'*probobj.ic;
            end
             obj.precomputeDistanceQuantities(romobj);
        end
        
       function [] = clearOnlineUpdateProperties(obj)
           
           obj.ROBcompon = [];
           obj.numswitch = 0;
           obj.prevSwitchState = [];
           obj.wr_switch = [];
           obj.LocBasisHist = [];
            
        end   
        
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
            obj.partialULoc = obj.probGNAT(obj.cLocBasis).ic;
            obj.partialUprevLoc = obj.partialULoc;
        end
        
        function  [] = computeOnlineMatrices(obj,phiR,phiJ,basisNum)
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
            
            
            if ~isempty(phiR) && ~isempty(phiJ)
                obj.A{basisNum,1} = pinv(phiJ(obj.sampleInd{basisNum,1},:));
                obj.B{basisNum,1} = phiJ'*phiR*pinv(phiR(obj.sampleInd{basisNum,1},:));
            end
        end
        
        function  [] = computeSampleNodes(obj,probobj,phiR,phiJ,cLocBasis)
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
            
            obj.sampleNodes{cLocBasis,1} = [];
            obj.sampleInd{cLocBasis,1} = [];
            
            %tempNodes = cell(obj.nBases,1);
            %tempInd = cell(obj.nBases,1);
            tempNodes = [];
            tempInd = [];
            
            %Make sure that we have not specified too few sample nodes.
            if obj.nSample(cLocBasis,1)*probobj.ndim < obj.nGreed(cLocBasis,1)
                nSampleOLD = obj.nSample(cLocBasis,1);
                obj.nSample(cLocBasis) = ceil(obj.nGreed(cLocBasis,1)/probobj.ndim);
                fprintf('Warning: not enough sample nodes! Increasing number of sample nodes from %d to %d\n',nSampleOLD,obj.nSample(cLocBasis));
                clear nSampleOLD;
            end
            
            P = min(obj.nSample(cLocBasis,1),obj.nGreed(cLocBasis,1));
            Q = floor(obj.nGreed(cLocBasis,1)/P)*ones(P,1);
            ind = [1:P] < mod(obj.nGreed(cLocBasis,1),P);
            Q(ind) = Q(ind) + 1;
            s = ceil(obj.nGreed(cLocBasis,1)/obj.nSample(cLocBasis,1));
            SS = floor(obj.nSample(cLocBasis,1)*s/obj.nGreed(cLocBasis,1))*ones(P,1);
            ind = [1:P] < mod(obj.nSample(cLocBasis,1),obj.nGreed(cLocBasis,1));
            SS(ind) = SS(ind)+1;
            
            %for kk = 1:obj.nBases
            QQ = 0;
            
            if isempty(phiR) || isempty(phiJ)
                return;
            end
            
            %if isempty(phiR{kk,1}) || isempty(phiJ{kk,1})
            %   continue;
            %end
            %Lines 6 and 7 from Algorithm 4 of Carlberg's thesis
            %R = phiR{kk,1}(:,1:max(Q));
            %J = phiJ{kk,1}(:,1:max(Q));
            R = phiR(:,1:max(Q));
            J = phiJ(:,1:max(Q));
            
            for p = 1:P
                QQ = QQ + Q(p);
                for s = 1:SS(p)
                    %Determine the set of indices not already in our sample mesh
                    %validNodes = setdiff(1:length(probobj.mesh.node),tempNodes{kk,1}');
                    validNodes = setdiff(1:length(probobj.mesh.node),tempNodes');
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
                    tempNodes = [tempNodes;n];
                    tempInd = [tempInd;probobj.node2ind(n)];
                    %tempNodes{kk,1} = [tempNodes{kk,1};n];
                    %tempInd{kk,1} = [tempInd{kk,1};probobj.node2ind(n)];
                end
                %Lines 13-16 of Algorithm 4 of Carlberg thesis
                for q = 1:Q(p)
                    a = phiR(tempInd,1:QQ)\phiR(tempInd,QQ+q);
                    b = phiJ(tempInd,1:QQ)\phiJ(tempInd,QQ+q);
                    
                    R(:,q) = phiR(:,QQ+q) - phiR(:,1:QQ)*a;
                    J(:,q) = phiJ(:,QQ+q) - phiJ(:,1:QQ)*b;
                    
                    %a = phiR{kk,1}(tempInd{kk,1},1:QQ)\phiR{kk,1}(tempInd{kk,1},QQ+q);
                    %b = phiJ{kk,1}( tempInd{kk,1},1:QQ)\phiJ{kk,1}(tempInd{kk,1},QQ+q);
                    
                    %R(:,q) = phiR{kk,1}(:,QQ+q) - phiR{kk,1}(:,1:QQ)*a;
                    %J(:,q) = phiJ{kk,1}(:,QQ+q) - phiJ{kk,1}(:,1:QQ)*b;
                end
            end
            
            %Sort the sample nodes and indices
            obj.sampleInd{cLocBasis,1} = unique(tempInd);
            obj.sampleNodes{cLocBasis,1} = unique(tempNodes);
            %obj.sampleInd{kk,1} = unique(tempInd{kk,1});
            %obj.sampleNodes{kk,1} = unique(tempNodes{kk,1});
            %end
        end
        
        function  [] = determineAdditionalNodes(obj,probobj,basisNum)
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
            
            obj.irstart{basisNum,1} = [];
            obj.jrow{basisNum,1} = [];
            obj.jdiag{basisNum,1} = [];
            
            if ~isempty(obj.sampleInd{basisNum,1})
                
                obj.irstart{basisNum,1} = 1;
                for i = 1:length(obj.sampleInd{basisNum,1})
                    %Find the nonzero entries in the appropriate row of Jstruct
                    temp = find(probobj.Jstruct(obj.sampleInd{basisNum,1}(i,1),:) ~= 0)';
                    %Store these indices in jrow
                    obj.jrow{basisNum,1} = [obj.jrow{basisNum,1};temp];
                    %Determine which of these values appended onto jrow
                    %corresponds to the diagonal.
                    obj.jdiag{basisNum,1} = [obj.jdiag{basisNum,1};temp == obj.sampleInd{basisNum,1}(i)];
                    %Store the index of the start of the next node
                    obj.irstart{basisNum,1} = [obj.irstart{basisNum,1}, length(obj.jrow{basisNum,1})+1];
                end
                %Make jdiag a column vector
                obj.jdiag{basisNum,1} = obj.jdiag{basisNum,1}(:);
            end
            
        end
        
        function  [] = computeSampleIndices(obj,probobj,phiR,phiJ,cLocBasis)
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
            
            obj.sampleNodes{cLocBasis,1} = [];
            obj.sampleInd{cLocBasis,1}   = [];
            
            %for kk = 1:obj.nBases
            if isempty(phiR) || isempty(phiJ)
                return;
            end
            %if isempty(phiR{kk,1}) || isempty(phiJ{kk,1})
            %    continue;
            %end
            % Algorithm 5 of Carlberg's 2010 Paper
            R = phiR(:,1);
            J = phiJ(:,1);
            
            nbarI = 0;
            m = 1;
            
            while nbarI < obj.nI(cLocBasis)
                %Determine the set of indices not already in our sample
                %indices
                validNodes = setdiff(1:length(R),obj.sampleInd{cLocBasis,1}');
                %Add up the contribution of the residual and
                %jacobian for each valid index
                temp = R(validNodes,1).^2 + J(validNodes,1).^2;
                %Determine the node with the max. value
                [~,ID] = max(temp);
                n = validNodes(ID);
                %Add the computed indices to the set
                obj.sampleInd{cLocBasis,1} = [obj.sampleInd{cLocBasis,1};n];
                K = checkResEvaluatedForFree(obj,probobj,n,nbarI,cLocBasis);
                obj.sampleInd{cLocBasis,1} = unique([obj.sampleInd{cLocBasis,1};K(:)]);
                
                nbarI = nbarI + 1 + length(K);
                m = m+1;
                pR = min(m-1,obj.nR(cLocBasis)); pJ = min(m-1,obj.nJ(cLocBasis));
                %Lines 10-11 of Algorithm 5 of Carlberg paper 2010

                a = phiR(obj.sampleInd{cLocBasis,1},1:pR)\phiR(obj.sampleInd{cLocBasis,1},m);
                b = phiJ(obj.sampleInd{cLocBasis,1},1:pJ)\phiJ(obj.sampleInd{cLocBasis,1},m);
                
                R = phiR(:,m) - phiR(:,1:pR)*a;
                J = phiJ(:,m) - phiJ(:,1:pJ)*b;
            end
            %Sort the sample nodes and indices
            obj.sampleInd{cLocBasis,1} = unique(obj.sampleInd{cLocBasis,1});
            %end
        end
        
        function  [K] = checkResEvaluatedForFree(obj,probobj,greedInd,iIndex,iBases)
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
            %iBases     - current local bases
            %
            %Outputs:
            %--------
            %K          - vector containing the set of indices that can be
            %             added for free
            %--------------------------------------------------------------
            %Code from David Amsallem 2010, with modifications made by
            %Matthew Zahr 2011 for integration with MORTestbed.
            
            switch lower(obj.addInd{iBases})
                case 'none'
                    K= [];
                case {'samevar', 'all'}
                    K = [];
                    switch obj.addInd{iBases}
                        case 'samevar'
                            Jref = find(probobj.Jstruct(greedInd,:));
                        case 'all'
                            determineAdditionalNodes(obj,probobj,iBases);
                            Jref = obj.jrow{iBases,1};
                    end
                    IndexComp = setdiff(1:probobj.config.ndof,obj.sampleInd{iBases,1});
                    
                    len = iIndex+1;
                    for i=1:length(IndexComp);
                        k = IndexComp(i);
                        J = find(probobj.Jstruct(k,:));
                        flag = true;
                        for j = 1:length(J)
                            if ~ismember(J(j),Jref)
                                flag = false;
                                break;
                            end
                        end
                        %J = find(probobj.Jstruct(k,:));
                        if flag
                            K = [K k];
                            len = len+1;
                        end
                        if (len>=obj.nI(iBases)) % indices limit is attained
                            if (len>obj.nI(iBases))
                                K(:,end) = [];
                            end
                            break;
                        end
                    end
                    
                    %                     for i=1:length(IndexComp);
                    %                         k = IndexComp(i);
                    %                         J = find(probobj.Jstruct(k,:));
                    %                         if (isempty(setdiff(J,Jref)))
                    %                             K = [K k];
                    %                         end
                    %                         if (iIndex+1+length(K)>=obj.nI(iBases)) % indices limit is attained
                    %                             if (iIndex+1+length(K)>obj.nI(iBases))
                    %                                 K(:,end) = [];
                    %                             end
                    %                             break;
                    %                         end
                    %                     end
            end
        end
        
        function  [] = NewtonRaphsonLocal(obj)
            %This function solves systems of the form:  f(x) = 0 using
            %Newton-Raphson's method for the GNAT class
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - GNAT object
            %
            %Outputs:
            %--------
            %There are no outputs.  The state vector, sv, is updated in the
            %GNAT handle class
            %--------------------------------------------------------------
            
            %Determine time iteration number plus 1
            itnump1 = obj.cTimeIter + 1;
            cbase = obj.cLocBasis;
            t = obj.time.T(1) + obj.time.dt*obj.cTimeIter;
            
            if obj.cTimeIter == 1
                obj.phiYhat = obj.phiYhat0(cbase,:)';
            end
            
            %On the first iteration, we determine if we need to do an
            %update (will be necessary if initial snapshot reference is not
            %the ic).  If we don't need to do an initial update, we advance
            %numswitch.
            if obj.cTimeIter == 1 && strcmpi(obj.basisUpdate,'border_fast') && ~strcmpi(obj.initref,'ic')
                obj.efficientOnlineSVDupdate(obj.cLocBasis,0);
            elseif obj.cTimeIter == 1 && strcmpi(obj.basisUpdate,'border_fast_approx') && ~strcmpi(obj.initref,'ic')
                obj.efficientOnlineSVDupdateApprox(obj.cLocBasis,0);
            elseif itnump1-1 == 1
                obj.numswitch = obj.numswitch + 1;
                if obj.cTimeIter == 1 && strcmpi(obj.basisUpdate,'border_fast_approx') && strcmpi(obj.initref,'ic')
                    [phiYhat_small] = reconstructSmallMaskforQ(obj,{obj.phiYhat0{:,cbase}}');
                    
                    for i1 = 1:obj.nBases
                        for i2 = 1:obj.nBases
                            if (i2>i1)
                                
                                obj.preCompDistQuant(i1,i2).h{obj.numswitch} = obj.preCompDistQuant(i1,i2).c'*(obj.Q*phiYhat_small);
                            end
                        end
                    end
                    for p=1:obj.nBases
                        obj.ROBcompon(obj.numswitch).phiYhat{p} = obj.phiYhat0{p,cbase};
                    end
                end
            end
            
            if obj.cTimeIter > 1
                if (obj.cLocBasis==obj.LocBasisHist(obj.cTimeIter-1))
                    %If the current basis is the same as the previous, set
                    %the initial guess to be the solution at the previous
                    %timestep (everything is masked).
                    obj.partialUprevLoc = obj.partialULoc;
                else
                    if strcmpi(obj.basisUpdate,'border_fast')
                        %If the basis changed and we are using fast
                        %updating, we need to reconstruct the solution of
                        %the previous timestep using the MASK of the NEW
                        %basis (which will be our initial guess for
                        %Newton's method).  Since we are using fast
                        %updating, the bases need to be reconstructed as
                        %well.
                        obj.partialUprevLoc = obj.probGNAT(cbase).ic;
                        for kk=1:obj.numswitch
                            if ~isempty(obj.ROBcompon(kk).alpha)
                                obj.partialUprevLoc = obj.partialUprevLoc + obj.computePhihatTimesVecFromComp(kk,cbase,obj.Ur_Loc{kk,1});
                            end
                        end
                        obj.partialULoc = obj.partialUprevLoc;
                        
                        %Update the (MASKED) bases
                        obj.efficientOnlineSVDupdate(obj.cLocBasis,obj.LocBasisHist(obj.cTimeIter-1));
                    elseif strcmpi(obj.basisUpdate,'border_fast_approx')
                        %If the basis changed and we are using fast
                        %updating, we need to reconstruct the solution of
                        %the previous timestep using the MASK of the NEW
                        %basis (which will be our initial guess for
                        %Newton's method).  Since we are using fast
                        %updating, the bases need to be reconstructed as
                        %well.
                        obj.partialUprevLoc = obj.probGNAT(cbase).ic;
                        for kk=1:obj.numswitch
                            obj.partialUprevLoc = obj.partialUprevLoc + obj.ROBcompon(kk).phiYhat{cbase}*obj.Ur_Loc{kk,1};
                        end
                        obj.partialULoc = obj.partialUprevLoc;
                        
                        %Update the (MASKED) bases
                        obj.efficientOnlineSVDupdateApprox(obj.cLocBasis,obj.LocBasisHist(obj.cTimeIter-1));
                    elseif strcmpi(obj.basisUpdate,'none')
                        %If the basis changed, we need to reconstruct the solution of
                        %the previous timestep using the MASK of the NEW
                        %basis (which will be our initial guess for
                        %Newton's method
                        obj.partialUprevLoc = obj.probGNAT(cbase).ic;
                        for kk=1:obj.nBases
                            %if ~isempty(obj.phiYhat{cbase,kk})
                            if ~isempty(obj.phiYhat{kk})
                                obj.partialUprevLoc = obj.partialUprevLoc + obj.phiYhat0{cbase,kk}*obj.Ur_Loc{kk,1};
                            end
                        end
                        obj.partialULoc = obj.partialUprevLoc;
                        
                        %Without basis updating, the new basis is simply
                        %the original basis at all masks.
                        obj.phiYhat = obj.phiYhat0(cbase,:)';
                    end
                end
            end
            if strcmpi(obj.basisUpdate,'none')
                obj.numswitch=cbase;
            end
            obj.switchHist(obj.cTimeIter) = obj.numswitch;
            
            %Set initial guess as last time step (note that sv is the
            %partial state). Initial guess for the sv = 0 for GNAT, but we
            %do need to update the partial state initial guess.
            
            %Determine the residual and jacobian based on initial guess
            [Rhat,JVhat] = obj.TimeScheme.TimeIntNLFuncGNATloc(obj.probGNAT,obj.partialULoc,obj.partialUprevLoc,t,cbase);
            %[Rhat,JVhat] = BackwardEulerNLFuncGNATloc(obj);
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
                %(unit step length)
                p = -RAbar\E;
                %If we use the below line, we have to precompute B*Rhat for
                %convergence check, then we don't need the
                %updateLocOperator calls.
                %p =-(obj.A{obj.cLocBasis,1}*JVhat)\(obj.B{obj.cLocBasis}*Rhat);
                %Linesearch if necessary
                
                if sum(isnan(p)) > 0 || sum(real(p)~=p) || rcond(RAbar) < 1e-16
                    obj.killflag=true;
                    return;
                end
                
                if obj.newt.linesrch.status
                    alpha = linesearch(@(alpha) obj.TimeScheme.TimeIntNLFuncGNATlocDir(obj.probGNAT,obj.partialULoc,obj.partialUprevLoc,t,cbase,alpha,p),obj.newt.linesrch.prop);
                    p = alpha*p;
                end
                
                
                obj.sv(1:length(p),itnump1) = obj.sv(1:length(p),itnump1) + p;
                obj.partialULoc = obj.partialUprevLoc + obj.phiYhat{cbase}*obj.sv(1:length(p),itnump1);
                %obj.partialULoc = obj.partialUprevLoc + obj.phiYhat{cbase,cbase}*obj.sv(1:length(p),itnump1);
                
                %Compute residual and jacobian with updated vector for next iteration
                [Rhat,JVhat] = obj.TimeScheme.TimeIntNLFuncGNATloc(obj.probGNAT,obj.partialULoc,obj.partialUprevLoc,t,cbase);
                %[Rhat,JVhat] = BackwardEulerNLFuncGNATloc(obj);
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
            obj.Ur_Loc{obj.numswitch,1} = obj.Ur_Loc{obj.numswitch} + obj.sv(1:length(p),itnump1);
            %obj.Ur_Loc{cbase,1} = obj.Ur_Loc{cbase,1} + obj.sv(1:length(p),itnump1);
            
            %Store the number of newton iterations at this time step
            obj.newt.iter(obj.cTimeIter+1) = i_N;
            
            %If the maxiumum newton iterations are reached, warn the user
            if ~obj.newt.quiet && (i_N==obj.newt.maxIter)
                disp(['*** Warning; Newton solver reached max number of iterations before convergence : res = ',...
                    num2str(norm(E,2)/tol),'  eps = ',num2str(obj.newt.eps), '***']);
            end
        end
        
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
            
            [Q,RAbar] = qr(obj.A{obj.cLocBasis,1}*JVhat,0);
            E = Q'*(obj.B{obj.cLocBasis}*Rhat);
            
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
            cbase = obj.cLocBasis;
            obj.JhatTemp{cbase}(obj.reconstJhatInd{cbase,1}) = Jhat;
            %JVhat = obj.JhatTemp{cbase}*obj.phiYhat{cbase,cbase};
            JVhat = obj.JhatTemp{cbase}*obj.phiYhat{cbase};
            
        end
        
        %Basis Updating
        function  [] = efficientOnlineSVDupdate(obj,newBasisNum,prevBasisNum)
            obj.numswitch=obj.numswitch+1;
            
            obj.Ur_Loc{obj.numswitch} = zeros(obj.nY(newBasisNum),1);
            if obj.numswitch > 1
                firstIndOfPrevBasis = find(prevBasisNum == obj.LocBasisHist,1,'first')+1;
                lastIndOfPrevBasis  = find(newBasisNum == obj.LocBasisHist,1,'last');
%                 obj.wr_switch{obj.numswitch-1}=obj.sv(1:obj.nY(prevBasisNum),firstIndOfPrevBasis);
%                 for i=firstIndOfPrevBasis+1:lastIndOfPrevBasis
%                     obj.wr_switch{obj.numswitch-1} = obj.wr_switch{obj.numswitch-1} + obj.sv(1:obj.nY(prevBasisNum),i);
%                 end
            end
            %if obj.numswitch == 2
            %obj.wr_switch{1} = obj.svdUpdateData(prevBasisNum).d - w_switch;
            %else
            %    obj.wr_switch{obj.numswitch-1} = w_switch - obj.svdUpdateData(prevBasisNum).F{obj.pprevBasisNum}*obj.prevSwitchState;
            %end
            %obj.prevSwitchState = w_switch;
            %obj.pprevBasisNum=prevBasisNum;
            
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

            if r<0% && abs(r) < 1e-7
                fprintf('Warning...Negative r. Skipping Update.');
                r = sqrt(abs(r))
                
                obj.ROBcompon(obj.numswitch).alpha = zeros(obj.nY(newBasisNum),1);
                
                obj.ROBcompon(obj.numswitch).beta  = cell(obj.nBases,1);
                obj.ROBcompon(obj.numswitch).N     = cell(obj.nBases,1);
                for i = 1:obj.nBases
                    obj.ROBcompon(obj.numswitch).beta{i} = zeros(obj.nY(newBasisNum),1);
                    obj.ROBcompon(obj.numswitch).N{i} = zeros(obj.nY(i),obj.nY(newBasisNum));
                    
                    if i == newBasisNum
                        obj.ROBcompon(obj.numswitch).N{i} = eye(obj.nY(newBasisNum));
                    else
                        obj.ROBcompon(obj.numswitch).N{i} = zeros(obj.nY(i),obj.nY(newBasisNum));
                    end
                    
                    obj.phiYhat{i} = obj.phiYhat0{i,newBasisNum};
                end
                return;
            else
                r=sqrt(r)
            end
            %             ae = obj.wReftmp{jsubl} - obj.icTmp;
            %             for i = 1:obj.cTimeIter
            %                 cbase=obj.LocBasisHist(i);
            %                 ae = ae - obj.phi0tmp{cbase}*obj.sv(1:obj.nY(cbase),i);
            %             end
            %             me = obj.phi0tmp{jsubl}'*ae;
            %             pe = ae - obj.phi0tmp{jsubl}*me;
            %             r = norm(pe);
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
            
            %             a=obj.icTmp*alpha_a;
            %             for i = 1:obj.nBases
            %                 a = a + obj.wReftmp{jsubl}*beta_a(i) + obj.phi0tmp{i}*n_a{i};
            %             end
            %             ae = obj.wReftmp{jsubl} - obj.icTmp;
            %             for i = 1:obj.cTimeIter
            %                 cbase=obj.LocBasisHist(i);
            %                 ae = ae - obj.phi0tmp{cbase}*obj.sv(1:obj.nY(cbase),i);
            %             end
            % %             ae = obj.wReftmp{jsubl}-obj.phi0tmp{prevBasisNum}*w_switch;
            %
            %             p=obj.icTmp*alpha_p;
            %             for i = 1:obj.nBases
            %                 p = p + obj.wReftmp{i}*beta_p(i) + obj.phi0tmp{i}*n_p{i};
            %             end
            %             me = obj.phi0tmp{jsubl}'*ae;
            %             pe = ae - obj.phi0tmp{jsubl}*me;
            %             re = norm(pe,2);
            %             pe=pe/norm(pe,2);
            
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
            
            for p = 1:obj.nBases
                obj.phiYhat{p} = obj.probGNAT(p).ic*obj.ROBcompon(obj.numswitch).alpha';
                for i = 1:obj.nBases
                    obj.phiYhat{p} = obj.phiYhat{p} + ...
                        obj.svdUpdateData(p).wRef(:,i)*obj.ROBcompon(obj.numswitch).beta{i}' + ...
                        obj.phiYhat0{p,i}*obj.ROBcompon(obj.numswitch).N{i};
                end
            end
        end
        
        function  [] = efficientOnlineSVDupdateApprox(obj,newBasisNum,prevBasisNum)%
            obj.numswitch=obj.numswitch+1;
            
            obj.Ur_Loc{obj.numswitch} = zeros(obj.nY(newBasisNum),1);
            if obj.numswitch > 1
                firstIndOfPrevBasis = find(prevBasisNum == obj.LocBasisHist,1,'first')+1;
                lastIndOfPrevBasis  = find(newBasisNum == obj.LocBasisHist,1,'last');
                obj.wr_switch{obj.numswitch-1}=obj.sv(1:obj.nY(prevBasisNum),firstIndOfPrevBasis);
                for i=firstIndOfPrevBasis+1:lastIndOfPrevBasis
                    obj.wr_switch{obj.numswitch-1} = obj.wr_switch{obj.numswitch-1} + obj.sv(1:obj.nY(prevBasisNum),i);
                end
            end
            
            jsubl = newBasisNum;
            
            ahat  = cell(obj.nBases,1);
            for p=1:obj.nBases
                ahat{p} = obj.svdUpdateData(p).wRef(:,jsubl) - obj.probGNAT(p).ic;
                for i=1:obj.numswitch-1
                    ahat{p} = ahat{p} - obj.ROBcompon(i).phiYhat{p}*obj.wr_switch{i};
                end
            end
            [ahat_small] = reconstructSmallMaskforQ(obj,ahat);
            [phiYhat0_small] = reconstructSmallMaskforQ(obj,{obj.phiYhat0{:,jsubl}}');
            
            m = phiYhat0_small'*(obj.Q*ahat_small);
            
            phat  = cell(obj.nBases,1);
            for p=1:obj.nBases
                phat{p} = ahat{p} - obj.phiYhat0{p,jsubl}*m;
            end
            [phat_small] = reconstructSmallMaskforQ(obj,phat);
            r = sqrt(phat_small'*(obj.Q*phat_small));
            
            for p=1:obj.nBases
                phat{p} = phat{p}/r;
            end
            
            K = [diag(obj.SingVals{jsubl}(1:obj.nY(jsubl))),zeros(obj.nY(jsubl),1);zeros(1,obj.nY(jsubl)),0] + ...
                [m;r]*[obj.svdUpdateData(jsubl).ColSumRSVtrunc',obj.svdUpdateData(jsubl).normOneMinusVVt];
            [C,~,~] = svd(K,0);
            C11 = C(1:end-1,1:end-1);
            c21 = C(end,1:end-1)';
            
            for p=1:obj.nBases
                obj.phiYhat{p} = obj.phiYhat0{p,jsubl}*C11 + phat{p}*c21';
                obj.ROBcompon(obj.numswitch).phiYhat{p} = obj.phiYhat{p};
            end
            
            % update quantities for fast distance computation
            [phiYhat_small] = reconstructSmallMaskforQ(obj,obj.phiYhat);
            for i1 = 1:obj.nBases
                for i2 = 1:obj.nBases
                    if (i2>i1)
                        obj.preCompDistQuant(i1,i2).h{obj.numswitch} = obj.preCompDistQuant(i1,i2).c'*(obj.Q*phiYhat_small);
                    end
                end
            end
            
        end
        
        function  [phiv] = computePhihatTimesVecFromComp(obj,switchNum,maskNum,vec)
            
            phiv = (obj.ROBcompon(switchNum).alpha'*vec)*obj.probGNAT(maskNum).ic;
            for j = 1:obj.nBases
                phiv = phiv + (obj.ROBcompon(switchNum).beta{j}'*vec)*obj.svdUpdateData(maskNum).wRef(:,j) + ...
                    obj.phiYhat0{maskNum,j}*(obj.ROBcompon(switchNum).N{j}*vec);
            end
            
        end
        
        function  [phiv] = computeFullPhiTimesVecFromComp(obj,romobj,switchNum,vec)
            
            phiv = (obj.ROBcompon(switchNum).alpha'*vec)*romobj.prob.ic;
            for j = 1:obj.nBases
                phiv = phiv + (obj.ROBcompon(switchNum).beta{j}'*vec)*romobj.svdUpdateData(j).wRef + ...
                    romobj.phi0{j}*(obj.ROBcompon(switchNum).N{j}*vec);
            end
            
        end
        
        function [] = precomputeDistanceQuantities(obj,romobj)
            % Pre-compute distance quantities
            if strcmpi(obj.basisUpdate,'border_fast_approx')
                precomputeLocalSimulationOnlineSVDupdateApprox(obj,romobj);
            else
                obj.preCompDistQuant = romobj.preCompDistQuant;
            end
        end
        
        function [] = precomputeLocalSimulationOnlineSVDupdateApprox(obj,romobj) %
            %This function precomputes quantities used in the local ROB
            %simulations.
            
            selectMaskForApproximateProduct(obj);
            
            for i1 = 1:obj.nBases
                for i2 = 1:obj.nBases
                    obj.preCompDistQuant(i1,i2).d =  norm(romobj.clusCenter(:,i1)-romobj.prob.ic,2)^2 - norm(romobj.clusCenter(:,i2)-romobj.prob.ic,2)^2;
                    obj.preCompDistQuant(i1,i2).c = 2*(romobj.clusCenter(obj.small_mask,i2)-romobj.clusCenter(obj.small_mask,i1));
                end
            end
            
        end
        
        function [] = createSPDOperator(obj,Snaps)%
            
            TrueProd = Snaps'*Snaps;
            p = length(obj.small_mask);
            PX = Snaps(obj.small_mask,:);
            %PX = Snaps(mask,:);
            clear Snaps;
            mu = 10^(-12);
            cvx_setup
            cvx_begin
            %cvx_solver sdpt3
            cvx_solver sedumi
            variable QQ(p,p) symmetric
            
            minimize( norm(PX'*QQ*PX - TrueProd,'fro')/norm(TrueProd,'fro'));
            
            subject to
            QQ == semidefinite(p) +mu*eye(p);
            cvx_end
            obj.Q = QQ;
            clear QQ;
        end
        
        function [] = setSPDOperator(obj,Q)
            obj.Q = Q;
        end
        
        function [] = selectMaskForApproximateProduct(obj)%
            
            mask = unique(obj.jrow{1,1});
            
            for kk=2:obj.nBases
                mask = unique([mask;obj.jrow{kk,1}]);
            end
            p = min(length(mask),obj.nSampleProd);
            small_mask_indices = ceil(linspace(1,length(mask),p));
            obj.small_mask = mask(small_mask_indices);
            
            for kk=1:obj.nBases
                obj.smallMaskIndexMatchForQ{kk} = [];
                obj.localMaskIndexMatchForQ{kk} = [];
            end
            
            for i=1:length(obj.small_mask)
                found = 0;
                for kk=1:obj.nBases
                    jrow_loc = unique(obj.jrow{kk,1});
                    for iEntry=1:length(jrow_loc)
                        if (obj.small_mask(i)== jrow_loc(iEntry))
                            found = 1;
                            obj.smallMaskIndexMatchForQ{kk} = [obj.smallMaskIndexMatchForQ{kk};i];
                            obj.localMaskIndexMatchForQ{kk} = [obj.localMaskIndexMatchForQ{kk};iEntry];
                            break
                        end
                    end
                    if (found)
                        break;
                    end
                end
                
            end
        end
        
        function [smallMasked] = reconstructSmallMaskforQ(obj,localMasked)%
            
            smallMasked = zeros(length(obj.small_mask),size(localMasked{1},2));
            for i=1:obj.nBases
                smallMasked(obj.smallMaskIndexMatchForQ{i},:) = localMasked{i}(obj.localMaskIndexMatchForQ{i},:);
            end
            
        end
        
        function  [] = reconstructLocalROB(obj,romobj,prevState,switchNum,localROB)
            
            jsubl = localROB;
            if (isempty(prevState))
                a = romobj.svdUpdateData(jsubl).wRef- romobj.prob.ic;
                
            else
                a = romobj.svdUpdateData(jsubl).wRef- prevState;
            end
            if (~isempty(prevState) && norm(a)>0)
                
                m = (romobj.phi0{jsubl}(obj.small_mask,:))'*(obj.Q*a(obj.small_mask));
                
                p = a - romobj.phi0{jsubl}*m;
                r = sqrt(p(obj.small_mask)'*(obj.Q*p(obj.small_mask)));
                
                p = p/r;
                
                K = [diag(obj.SingVals{jsubl}(1:obj.nY(jsubl))),zeros(obj.nY(jsubl),1);zeros(1,obj.nY(jsubl)),0] + ...
                    [m;r]*[obj.svdUpdateData(jsubl).ColSumRSVtrunc',obj.svdUpdateData(jsubl).normOneMinusVVt];
                
                
                [C,~,~] = svd(K,0);
                C11 = C(1:end-1,1:end-1);
                c21 = C(end,1:end-1)';
                
                obj.reconstructedLocalROB{switchNum} = romobj.phi0{jsubl}*C11 + p*c21';
            else
                obj.reconstructedLocalROB{switchNum} = romobj.phi0{jsubl};
            end
        end

        %Unused?
        function  [] = adjustProp_nRnJnI(obj,nRnew,nJnew,nInew,basenum)
            obj.nR(basenum) = nRnew;
            obj.nJ(basenum) = nJnew;
            obj.nI(basenum) = nInew;
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
                'partialU','partialUprev','JhatTemp','phiYhatInd',...
                'reconstJhatInd','uniqueJROWind','jdiagHat','Ur_Loc',...
                'phiYhatUnion','phi'};
            
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
            
            newGnat = locGNAT([],[],[],obj);
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
    end
end

classdef WROM < handle
    %The WROM class is used for the Wavelet-based reduced-order modeling
    %framework developed by the UTEP research team
    %   Detailed explanation goes here    
    properties (SetAccess = private, GetAccess = public)
        %Input properties
        fileRoot = [pwd,SlashBS];
        prob;                   %Problem object
        Rtype = 'pg';           %ROM type: 'G' or 'PG'
        precompFlag = false;
        
        id;         %This is the id of the Rtype ROM module
        cfgid;      %This is the id of the configuration that
        snaps;      %structure with fields: coll, int, num, dist, distparams, distseed
        time; 
        newt;       %structure with fields: maxIter, eps, nIter (average # of newton steps)
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
        
        % Rei: 4/11/14  Wavelet variables
        % For list of Wavelet families visit
        % http://www.mathworks.com/help/wavelet/ref/waveletfamilies.html
        % http://www.mathworks.com/help/wavelet/ref/wfilters.html
        

        waveName;       % name of the wavelet transform
        waveLvl;        % decomposition level of the wavelet transform
        wavenY;         % desired dimension for low-pass filter (nY for wavelet)
        highPass;       % high-pass filter of the wavelet transform
        lowPass;        % high-pass filter of the wavelet transform 
        Wmatrix;        % wavelet matrix
        MaxLevel;       % maximum evel of decomposition

                
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
        LocBasisHist;
        clusCenter;
        preCompDistQuant;
        UrLoc;
        cTimeIter; %current time iteration number
        resBound = [];
        errBound = [];
        
        numExecute=0;
    end
    
    properties (Hidden = true, SetAccess = private, GetAccess = public)
        %Global variable text
        GVARtxt = [];           % First public property
        killflag = false;
        
        numres=0;
        printLevel=1;
        
    end  
    
    methods
        %Constructor
        function [obj] = WROM(ROMfile,romtype,romid,probobj)
            %This is the constructor of the WROM class.  It reads the
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
                if isa(oldobj,'WROM')
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
            if strcmpi(romtype(1),'g')
                obj.Rtype = 'Galerkin';
                ROMtext    = readInFile('GROM',ROMfile,1);
            elseif strcmpi(romtype(1),'p')
                obj.Rtype = 'Petrov-Galerkin';
                ROMtext    = readInFile('PGROM',ROMfile,1);
            end
            

            
            %Extract and store the configuration id to which this ROM
            %instance is associated with.
            obj.cfgid = probobj.config.id;
            
            
            %Set up number of dofs
            obj.ndof = probobj.config.ndof;
            
            %Determine the ROM properties based on the text in the ROM
            %section of the input file
            determineROM(obj,ROMtext,[obj.GVARtxt;VARtext],romid);            

            % 4/29/14 Build wavelet basis
            buildWaveletBasis(obj);

            
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

%           %%%%%%%%% Determine time integration %%%%%%%%%
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
            
            
%           %%%%%%%%% Determine Wavelet-Based ROM parameters %%%%%%%%%
            obj.waveName = extractInputROM(ROMtext,N,obj.id,'waveName','haar');
            obj.waveLvl  = extractInputROM(ROMtext,N,obj.id,'waveLvl',1);
            obj.wavenY   = extractInputROM(ROMtext,N,obj.id,'wavenY',8);
            
           

%             %%%%%%%%% Initialize SVD update structure and determine buffer %%%%%%%%%
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
        
                    
        function [] = buildWaveletBasis(obj)
            % 4/27/14 Rei: determine the Wavelet basis for nBases = 1

            %W   = get_waveMatrix(obj.ndof,obj.waveLvl,obj.waveName);
            waveletFactory(obj);
            W = obj.Wmatrix;
            % extract low pass matrix
%             L = W(1:obj.ndof/2^obj.waveLvl,:);
            L = obj.lowPass;
            PSI = L';
            % extract high pass matrix
%             obj.highPass = W(obj.ndof/2^obj.waveLvl + 1:obj.ndof,:);
            obj.phi      = PSI;%            
            
            str = sprintf('Using wavelet %s with level %d\n',obj.waveName, obj.waveLvl);
            fprintf(str);
            str2= sprintf('Full-space dimension: %d,\tReduced-space dimension: %d\n',obj.ndof, size(PSI,2));
            fprintf(str2);
        end

        
        function waveletFactory(obj)
            % Leobardo update: 5/23/14 (handle odd and even dimensions)
            dimension  = obj.ndof;
            filtername = obj.waveName;
            level      = obj.wavenY;
            lowpassDim = obj.wavenY;
            
            lowpassDim
            
            [Lo_D,Hi_D,Lo_R,Hi_R] = wfilters(filtername);

            k=0;
            while dimension/(2^k) > lowpassDim
                k=k+1;
            end
            
            k
            
            %filter=zeros(1,dimension);
            index=[1:dimension];
            index1=index;
            filter=zeros(1,dimension);
            filter(1:length(Lo_D))=fliplr(Lo_D);

            % I am obtaining the low-pass block matrix
            for i=1:ceil(dimension/2)
                W(i,index)=filter;
                index1(1:(end-2))=index(3:end);
                index1((end-1):end)=index(1:2);
                index=index1;
            end

            index=[1:dimension];
            index1=index;
            filter=zeros(1,dimension);
            filter(1:length(Hi_D))=fliplr(Hi_D);

            % I am obtaining the Hight-pass block matrix
            for i=1:floor(dimension/2)
                W(ceil(dimension/2)+i,index)=filter;
                index1(1:(end-2))=index(3:end);
                index1((end-1):end)=index(1:2);
                index=index1;
            end
            
            
            obj.Wmatrix  = W^k;
            obj.lowPass  = obj.Wmatrix(1:level,:);
            obj.highPass = obj.Wmatrix(level+1:dimension,:);
            
            
%             % Leobardo update: 5/23/14 (handle odd and even dimensions)
%             dimension  = obj.ndof;
%             filtername = obj.waveName;
%             level      = obj.waveLvl;
%             [Lo_D,Hi_D,Lo_R,Hi_R] = wfilters(filtername);
% 
%             %filter=zeros(1,dimension);
%             index=[1:dimension];
%             index1=index;
%             filter=zeros(1,dimension);
%             filter(1:length(Lo_D))=fliplr(Lo_D);
% 
%             % I am obtaining the low-pass block matrix
%             for i=1:ceil(dimension/2)
%                 W(i,index)=filter;
%                 L(i,index) = W(i,index);
%                 index1(1:(end-2))=index(3:end);
%                 index1((end-1):end)=index(1:2);
%                 index=index1;
%             end
%     
%             obj.lowPass = L;
%             
%             index=[1:dimension];
%             index1=index;
%             filter=zeros(1,dimension);
%             filter(1:length(Hi_D))=fliplr(Hi_D);
% 
%             % I am obtaining the High-pass block matrix
%             for i=1:floor(dimension/2)
%                 W(ceil(dimension/2)+i,index)=filter;
%                 H(ceil(dimension/2)+i,index)= W(ceil(dimension/2)+i,index);
%                 index1(1:(end-2))=index(3:end);
%                 index1((end-1):end)=index(1:2);
%                 index=index1;
%             end
%             
%             obj.highPass = H;
% 
%             W=W^level;
%             obj.Wmatrix = W;            
            
            
            
        end
        
        
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
            res0 = norm(R,2);
            [tol,tolIt] = determineConvergeCriterion(obj,res0);
            
            indexAdj=1;
            if obj.printLevel > 1 
                fprintf('---- Newton Step    # --------- ||R|| ------------ ||W''*R|| ----------- ||du|| ------------- tol -------------- tol_du ------ snap_coll ----------\n');
            end
            for i_N = 1:obj.newt.maxIter
                %Solve for the search direction using Newton or Gauss-Newton
                if lower(obj.Rtype(1)) == 'g'
                    newR = obj.phi'*R;
                    du = -((obj.phi'*J*obj.phi)\newR);
                    conv = norm(newR,2);
                    %conv(i_N) = norm(newR,2);
                else
                    
%                     % 4/28/14 Rei: use the Wavelet basis for nBases = 1
%                     (for NLTrans)
%                    
%                     W   = get_waveMatrix(obj.ndof,obj.waveLvl,obj.waveName);
%                     % extract low pass matrix
%                     L = W(1:obj.ndof/2^obj.waveLvl,:);
%                     PSI = L';
%                     % extract high pass matrix
%                     obj.highPass = W(obj.ndof/2^obj.waveLvl + 1:obj.ndof,:);
%                     str = sprintf('Using wavelet %s with level %d\n',obj.waveName, obj.waveLvl);
%                     fprintf(str);
%                     obj.phi = PSI;%


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
%                     alpha = linesrchBackNewton(@(alpha) obj.TimeScheme.TimeIntNLFuncDir(obj.prob,obj.sv(:,itnump1-indexAdj),obj.sv(:,itnump1-1),t,alpha,du,obj));
                    alpha = linesearch(@(alpha) obj.TimeScheme.TimeIntNLFuncDir(obj.prob,obj.sv(:,itnump1-indexAdj),obj.sv(:,itnump1-1),t,alpha,obj.phi*du),obj.newt.linesrch.prop);
                    du = alpha*du;
                end
                %Reconstructed search direction
                p = obj.phi*du;
                
                obj.sv(:,itnump1)=obj.sv(:,itnump1-indexAdj)+p;
                
                if obj.printLevel > 2
                    fprintf('---- Newton Step %4i ----- %10.7e ----- %10.7e ----- %10.7e ----- %10.7e ----- %10.7e ------- %1i --------------\n',i_N,norm(R,2),conv,norm(du,2),tol,tolIt,obj.saveNL);
                end
            
                
                indexAdj=0;
                if checkConverge(obj,conv,du,tol,tolIt)
                    break;
                end
                
                %Compute residual and jacobian with updated vector for next iteration
                [R,J] = obj.TimeScheme.TimeIntNLFunc(obj.prob,obj.sv(:,itnump1),obj.sv(:,itnump1-1),t);
            end
                
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
                
                
                %Compute residual and jacobian with updated vector for next iteration
                [R,J] = obj.TimeScheme.TimeIntNLFunc(obj.prob,obj.sv(:,itnump1),obj.sv(:,itnump1-1),t);
                
                indexAdj = 0;
                %Stop iterating once convergence criteria is met
                if checkConverge(obj,conv,p,tol,tolIt)
                    break;
                end
            end
%             %Write the last nonlinear snapshot.  We note that we are using
%             %J^(k+1)*0 if using Snapshot 1.5 or 2 (we don't have a
%             %search direction p^(k+1) because we converged!).  It shouldn't
%             %matter much because ||p|| should be small since we converged,
%             %i.e. the next step shouldn't take us away from the solution!
%             obj.writeNonlinearSnapshot(R,J,zeros(size(p)),cLocBasis);
%             
            %Store the number of newton iterations at this time step
            obj.newt.iter(obj.cTimeIter) = i_N;
            
            %If the maxiumum newton iterations are reached, warn the user
            if ~obj.newt.quiet && (i_N==obj.newt.maxIter)
                disp(['*** Warning; Newton solver reached max number of iterations before convergence : res = ',...
                    num2str(norm(R,2)/tol),'  eps = ',num2str(obj.newt.eps), '***']);
            end
        end %Done
%         
%         function  [] = NewtonRaphson_PrecompGal(obj)
%             %This function solves systems of the form:  f(x) = 0 using
%             %Newton-Raphson's method
%             %--------------------------------------------------------------
%             %Inputs:
%             %-------
%             %obj     - ROM object
%             %
%             %Outputs:
%             %--------
%             %No outputs. The state vector at the current time
%             %(indicated by the cTimeIter property) is stored in the FOM
%             %handle class.
%             %--------------------------------------------------------------
%             
%             %Determine time iteration number plus 1
%             itnump1 = obj.cTimeIter + 1;
%             
%             %Extract current basis
%             if obj.nBases > 1
%                 cLocBasis = obj.LocBasisHist(itnump1-1,1);
%             else
%                 cLocBasis = 1;
%                 obj.numswitch=1;
%             end
%                         
%             %%%%%%Basis Updating NOT SUPPORTED here (yet)%%%%%%
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             
%             %Determine the residual and jacobian based on initial guess
%             t = obj.time.T(1) + obj.time.dt*obj.cTimeIter;
%             [R,J] = obj.TimeScheme.TimeIntNLFunc(obj.prob,obj.sv(1:obj.nY(cLocBasis),itnump1-1),obj.sv(1:obj.nY(cLocBasis),itnump1-1),t,true);
%             %[R,J] = obj.prob.ResJacPrecomp(obj.sv(1:obj.nY(cLocBasis),itnump1-1),t,cLocBasis);
%             
%             %Determine convergence criteria to use
%             [tol,tolIt] = determineConvergeCriterion(obj,norm(R,2));
% 
%             indexAdj = 1;
%             for i_N = 1:obj.newt.maxIter
%                 dy = -J\R;
%                 
%                 obj.sv(1:obj.nY(cLocBasis),itnump1) = obj.sv(1:obj.nY(cLocBasis),itnump1-indexAdj) + dy;
%                 obj.UrLoc{obj.numswitch,1} = obj.UrLoc{cLocBasis,1} + dy;
% 
%                 %Compute residual and jacobian with updated vector for next iteration
%                 %[R,J] = obj.prob.ResJacPrecomp(obj.sv(1:obj.nY(cLocBasis),itnump1),t,cLocBasis);
%                 [R,J] = obj.TimeScheme.TimeIntNLFunc(obj.prob,obj.sv(1:obj.nY(cLocBasis),itnump1),obj.sv(1:obj.nY(cLocBasis),itnump1-1),t,true);
%             
%                 indexAdj = 0;
%                 %Stop iterating once convergence criteria is met
%                 if checkConverge(obj,norm(R),dy,tol,tolIt)
%                     break;
%                 end
%             end
%             %Store the number of newton iterations at this time step
%             obj.newt.iter(obj.cTimeIter) = i_N;
%             
%             %If the maxiumum newton iterations are reached, warn the user
%             if ~obj.newt.quiet && (i_N==obj.newt.maxIter)
%                 disp(['*** Warning; Newton solver reached max number of iterations before convergence : res = ',...
%                     num2str(norm(R,2)/tol),'  eps = ',num2str(obj.newt.eps), '***']);
%             end
%         end
                       
%         function  [] = precomputeLocalSimulation(obj)
%             %This function precomputes quantities used in the local ROB
%             %simulations.
%             %--------------------------------------------------------------
%             %Inputs:
%             %-------
%             %trainIC - ndof x 1 vector containing the training initial
%             %          condition
%             %
%             %Outputs:
%             %--------
%             %This function has no outputs.  It stores the quantites f and S
%             %in the ROM object, which are precomputed for the local
%             %simulation.
%             %--------------------------------------------------------------
%             
%             if obj.nBases == 1, return; end
%             
%             if strcmpi(obj.basisUpdate,'none')
%                 for i1 = 1:obj.nBases
%                     for i2 = 1:obj.nBases
%                         obj.preCompDistQuant(i1,i2).d = norm(obj.clusCenter(:,i1)-obj.prob.ic,2)^2 - norm(obj.clusCenter(:,i2)-obj.prob.ic,2)^2;
%                         for j = 1:obj.nBases
%                             obj.preCompDistQuant(i1,i2).g{j} = 2*(obj.clusCenter(:,i2)-obj.clusCenter(:,i1))'*obj.phi0{j};
%                         end
%                     end
%                 end
%             elseif strcmpi(obj.basisUpdate,'border_fast') || strcmpi(obj.basisUpdate,'border_fast_approx')
%                 for i1 = 1:obj.nBases
%                     for i2 = 1:obj.nBases
%                         obj.preCompDistQuant(i1,i2).d = norm(obj.clusCenter(:,i1)-obj.prob.ic,2)^2 - norm(obj.clusCenter(:,i2)-obj.prob.ic,2)^2;
%                         obj.preCompDistQuant(i1,i2).e = 2*(obj.clusCenter(:,i2)-obj.clusCenter(:,i1))'*obj.prob.ic;
%                         for j = 1:obj.nBases
%                             obj.preCompDistQuant(i1,i2).f{j} = 2*(obj.clusCenter(:,i2)-obj.clusCenter(:,i1))'*obj.svdUpdateData(j).wRef;
%                             obj.preCompDistQuant(i1,i2).g{j} = 2*(obj.clusCenter(:,i2)-obj.clusCenter(:,i1))'*obj.phi0{j};
%                         end
%                     end
%                 end
%             end
%         end
        
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
        
        
        % 4/17/14 Rei: define wavelet basis here. (overwrite POD basis)
        function  [] = computePOD(obj,flag,svdOptions,rootFname)
            %Build nLocal local POD bases.
            
%             %Set defaults
%             if nargin < 3
%                 svdOptions=[];
%             end
%             if nargin < 4 || isempty(rootFname)
%                 rootFname = obj.fileRoot;
%             end
%             
%             if obj.nBases == 1
%                 switch lower(flag)                        
%                     case 'fromclustsnaps'
%                         %Read clustered snapshots, compute POD (using SVD
%                         %options specified; default is exact SVD)
%                         X = obj.readClusSnapFiles(1);
%                         [LSingVecs,SVs,RSingVecs]  = pod(X,[],svdOptions);
%                         obj.SingVals = diag(SVs); clear SVs;
%                         
%                         %Write SVD factors to files (not V because we don't
%                         %need it in the global case)
%                         fid = fopen([rootFname,'svdFactU1.bin'],'wb');
%                         fwrite(fid,size(LSingVecs)','double');
%                         fwrite(fid,LSingVecs,'double');
%                         fclose(fid);
%                         
%                         fid = fopen([rootFname,'svdFactS1.bin'],'wb');
%                         fwrite(fid,size(obj.SingVals)','double');
%                         fwrite(fid,obj.SingVals,'double');
%                         fclose(fid);
%                     case 'fromsvd'
%                         %If we want to compute POD directly from an SVD we
%                         %already have available, read the SVD files
%                         [LSingVecs,obj.SingVals,RSingVecs]  = obj.readSVDFiles(1,rootFname);
%                 end
%                 
%                 if obj.nYflag == 1 %(nYflag = 1 means we use specified nY)
%                     %If nY specified, just use standard POD and
%                     %return
%                     newNY = min(obj.nY,size(LSingVecs,2));
%                     if newNY ~= obj.nY
%                         warning(['nY changed from ',num2str(obj.nY),' to ',num2str(newNY),' due to available number of snapshots!']);
%                         obj.nY = newNY;
%                     end
%                     obj.phi = LSingVecs(:,1:obj.nY);
%                     
%                     
%                     % 4/17/14 Rei: use the Wavelet basis for nBases = 1
%                    
%                     W   = get_waveMatrix(obj.ndof,obj.waveLvl,obj.waveName);
%                     % extract low pass matrix
%                     L = W(1:obj.ndof/2^obj.waveLvl,:);
%                     PSI = L';
%                     % extract high pass matrix
%                     obj.highPass = W(obj.ndof/2^obj.waveLvl + 1:obj.ndof,:);
%                     str = sprintf('Using wavelet %s with level %d\n',obj.waveName, obj.waveLvl);
%                     fprintf(str);
%                     obj.phi = PSI;%


%                 for iLocBasis = 1:obj.nBases
%                     % Rei: 4/1/14 (test 3 wavelets as bases)
%                     % Define Wavelet basis
%                     
%                     if iLocBasis ==1
%                         W   = get_waveMatrix(128,obj.waveLvl,'rbio3.7');
%                         % extract low pass matrix
%                         L = W(1:128/2^obj.waveLvl,:);
%                         PSI = L';
%                         % extract high pass matrix
%                         H = W(128/2^obj.waveLvl + 1:128,:);
%                         fprintf('using wavelet basis 1\n');
%                         obj.phi0{iLocBasis,1} = PSI;%PSI(:,1:16);
%                         obj.highPass{iLocBasis} = H;
%                     elseif iLocBasis == 2
%                         W   = get_waveMatrix(128,obj.waveLvl,'rbio3.7');
%                         L = W(1:128/2^obj.waveLvl,:);
%                         PSI = L';
%                         % extract high pass matrix
%                         H = W(128/2^obj.waveLvl + 1:128,:);
%                         fprintf('using wavelet basis 2\n');
%                         obj.phi0{iLocBasis,1} = PSI;%PSI(:,1:16);
%                         obj.highPass{iLocBasis} = H;
%                     else
%                         W   = get_waveMatrix(128,obj.waveLvl,'haar');
%                         L = W(1:128/2^obj.waveLvl,:);
%                         PSI = L';
%                         % extract high pass matrix
%                         H = W(128/2^obj.waveLvl + 1:128,:);
%                         fprintf('using wavelet basis 3\n');
%                         obj.phi0{iLocBasis,1} = PSI;%PSI(:,1:16);  % 128 x 16
%                         obj.highPass{iLocBasis} = H;
%     
%                     end
%                 end


                obj.numswitch = 0;
                
                disp([' Generated local POD bases of size ',num2str(obj.nY(:)'),', respectively']);
                
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
            
%             obj.openCloseResJacFiles('open');
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
                
                        %Use Newton's Method to solve nonlinear system of equations
                        %and store: state vector, residual and jacobian snapshots,
                        %and number of newton iterations required
                        obj.NewtonRaphson();
                        if obj.killflag, warning('Encountered NaN at Time Step %i. Simulation Killed',obj.cTimeIter); return; end;
                       
                        %Check for steady convergence
                        if obj.time.steadyconverge
                            if norm(obj.sv(:,i_t+1)-obj.sv(:,i_t)) < obj.time.steadyconverge*norm(obj.sv(:,i_t+1))
%                                 tmp = obj.TimeScheme;
%                                 tmp2 = obj.newt.maxIter;
                                
%                                 obj.cTimeIter = i_t+1;
%                                 obj.TimeScheme=obj.lastTimeScheme;
%                                 obj.newt.maxIter=55;
%                                 
%                                 obj.NewtonRaphson();
%                                 
%                                 obj.TimeScheme = tmp;
%                                 obj.newt.maxIter = tmp2;
%                                 
                                obj.ontime = toc(tROMstart)+OldTime;
                                obj.time.nstep=i_t;
                                obj.time.T(2) = obj.time.T(1) + i_t*obj.time.dt;
                                obj.sv(:,i_t+2:end)=[];
%                                 obj.openCloseResJacFiles('close');
%                                 obj.newt.avgIter = mean(obj.newt.iter(1:obj.time.nstep));
                                
%                                 [R,J] = obj.prob.ResJac(obj.sv(:,end));
%                                 converge=norm(obj.phi'*(J'*R));
%                                 fprintf('Norm of phi''*J''*R = %f \n',converge);
%                                 if converge > 1e-2
%                                     obj.killflag = true;
%                                 end
                                
                                return;
                            end
                            
%                             if i_t == nstep-1
%                                 tmp = obj.TimeScheme;
%                                 tmp2 = obj.newt.maxIter;
%                                 
%                                 obj.cTimeIter = i_t+1;
%                                 obj.TimeScheme=obj.lastTimeScheme;
%                                 obj.newt.maxIter=55;
%                                 
%                                 obj.NewtonRaphson();
%                                 
%                                 obj.TimeScheme = tmp;
%                                 obj.newt.maxIter = tmp2;
%                                 
%                                 obj.ontime = toc(tROMstart)+OldTime; %Record simulation time
%                                 obj.openCloseResJacFiles('close');
%                                 obj.newt.avgIter = mean(obj.newt.iter);
%                                 fprintf('steady converge level = %f\n',norm(obj.sv(:,end)-obj.sv(:,end-1))/norm(obj.sv(:,end)));
%                                 [R,J] = obj.prob.ResJac(obj.sv(:,end));
%                                 converge=norm(obj.phi'*(J'*R));
%                                 fprintf('Norm of phi''*J''*R = %f \n',converge);
% %                                 if converge > 1e-2
% %                                     obj.killflag = true;
% %                                 end
%                                 
%                                 return;
%                             end
                        end
                    end
                    %[R,J] = obj.prob.ResJac(obj.sv(:,2)); % BUG: need
                    %second argument
                    %converge=norm(obj.phi'*(J'*R));
                    %fprintf('Norm of phi''*J''*R = %f \n',converge);
                    %fprintf('Number of Newton iterations = %i\n',obj.newt.iter);
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
%                         obj.LocBasisHist(obj.cTimeIter,1) = computeDistanceToCenters(obj,obj.UrLoc,obj.cTimeIter);
%                             case 'a_priori'
%                                 obj.LocBasisHist(obj.cTimeIter,1) = floor(obj.nBases*obj.cTimeIter/nstep)+1;
%                         end
%                         cbase=obj.LocBasisHist(obj.cTimeIter,1);
                        
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
%             obj.openCloseResJacFiles('close');
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
 
%                 obj.projerr = ColumnwiseNorm(vecs - obj.phi*(obj.phi'*vecs),normType);
                
                
                if normalizeFlag
                    obj.projerr = obj.projerr./ColumnwiseNorm(fomobj.sv,2);
                    fprintf('normalized projection error\n')
                end
                
                
                obj.projerr = ColumnwiseNorm(obj.highPass*fomobj.sv,normType);
%                 obj.projerr = ColumnwiseNorm(obj.highPass*vecs,normType);
                
                %Plot results
                if plotFlag
                    %figure;
                    plot(linspace(fomobj.time.T(1),fomobj.time.T(2),fomobj.time.nstep+1),obj.projerr,...
                        'k','linewidth',2);
                end
                xlabel('time','interpreter','latex'); ylabel('Projection Error','interpreter','latex');
                return;
            end
            
            %Compute projection error for multiple basis case
            for i = 1:obj.nBases
                vecs = bsxfun(@minus,fomobj.sv,obj.prob.ic);
                obj.projerr(:,i) = ColumnwiseNorm(vecs - obj.phi0{i}*(obj.phi0{i}'*vecs),normType);
%                 obj.projerr(:,i) = ColumnwiseNorm(obj.highPass{i}*vecs,normType);
                
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

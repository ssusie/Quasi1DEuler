classdef TPWL < handle
    
    properties (SetAccess = private, GetAccess = public)
        id;
        prob;
        nY;
        nLinPt; %scalar indicating the number of linearization points used for the current TPWL model
        nLinPtMax = inf; %Only applies for ResDist and TrajDist
        LinPtAlgo;
        LinPtAlgoParam;
        LinPtInt;
        basis;
        snapcoll = 'ref_init';
        estNextSv = 'prev'; %or linextrap
        beta = 25; %constant scalar for determining weight (value of 25 recommended in Rewienski Thesis)
        epsilon = 100e-2; %Establish parameter for time stepping
        %procedure that represents tolerance
        %for convergence of an initial guess of
        %state vector to solution (if tolerance
        %not met, need to decrease step size)
        
        FOMsize;
        time;
        
        ontime;
        LinPt;
        z0  %reduced initial condition
        zi; % matrix whose columns contain the project of the full state vector at the linearization points
        sv; %(reduced) state vector matrix
        phi;
        A_ri;  % nLinPt x 1 cell array whose entries contain the
               % reduced dynamics matrix of the linearized model at
               % each linearization point 
        B_hat; % nLinPt x 1 cell array whose entries contain the
               % reduced input matrix of the linearized model at
               % each linearization point 
        gamma; % ndof x nLinPt matrix containing the projected
               % "effective input" of LDS at each linearization point
        cTimeIter; %current time iteration number
    end
    
    properties (Hidden = true, SetAccess = private, GetAccess = public)
       GVARtxt; 
    end
    
    methods
        function [obj] = TPWL(ROMfile,TPWLid,probobj)
            %This is the constructor of the TPWL class.  It reads the
            %appropriate entries from the input file and stores them in the
            %class instance.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %ROMfile - string indicating the filename of the rom file to
            %          use
            %TPWLid  - id number of the TPWL object to use
            %proobj  - Problem object
            %
            %Outputs:
            %--------
            %obj     - instance of the TPWL object that was constructed
            %--------------------------------------------------------------
            
            if nargin == 0
                return;
            end
            
            %Store problem in model
            obj.prob = probobj;
            
            %Default time structure is the same as the FOM
            obj.time = obj.prob.config.time;
            %%%%%%%%%%%%%%%%%%%%Set default values%%%%%%%%%%%%%%%%%%%%
            dfname = [probobj.config.defaultFile,'.rom'];
            dVARtext = readInFile('VAR',dfname,1);
            %Extract the FOM text from the cfg file
            dTPWLtext = readInFile('TPWL',dfname,1);
            %Determine the TPWL properties based on the text in the ROM
            %section of the input file 
            determineTPWL(obj,dTPWLtext,dVARtext,1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %Extract the TPWL text from the rom file
            TPWLtext = readInFile('TPWL',ROMfile,1);
            %Extract parameters from ROM file
            obj.GVARtxt = probobj.config.GVARtxt;
            VARtext = readInFile('VAR',ROMfile,1);
            
            %Determine the TPWL properties based on the text in the TPWL section of the
            %input file
            determineTPWL(obj,TPWLtext,[obj.GVARtxt;VARtext],TPWLid);
        end
        
        function  [] = determineTPWL(obj,TPWLtext,VARtxt,id)
            %This function converts the text from the TPWL object into a structure
            %whose fields correspond to the TPWL properties
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj      - TPWL object
            %TPWLtxt  - a string containing the text from the TPWL object
            %id       - integer indicating the id number of the TPWL object
            %           of interest from the input file
            %
            %Outputs:
            %--------
            %There are no outputs.  Properties are stored in the TPWL
            %handle class.
            %--------------------------------------------------------------
            
            %Extract variables from VARtxt and evaluate them in the current
            %workspace.
            if ~isempty(VARtxt)
                for i = 1:length(VARtxt)
                    eval([VARtxt{i},';']);
                end
            end
            
%             %Determine if any properties are specified more than once
%             determineDuplicate(TPWLtext);
%             %Find all line break indicators in TPWLtext
%             TPWLdelim = findstr(TPWLtext,'//');
            
            %Set TPWL id number
            tmp = determinePropText(TPWLtext,'id');
            if isempty(tmp)
                error('In ROM file (TPWL section), the ID must be specified.  This controls the connection between different aspects of the program.');
            end
            idvec = eval(tmp);
            idvec = idvec(:);
            N = size(idvec,1);
            
            %Determine the entry in the input file that correspond to id
            ind = find(idvec == id);
            if isempty(ind)
               error('Can only initialize ROMs with ids defined in .rom file'); 
            end
            
            if N > 1
                obj.id = idvec(ind);
            elseif N == 1
                obj.id = idvec;
            else
                error('In ROM file (TPWL section), the ID must be specified with a scalar (double or int).  This controls the connection between different aspects of the program.');
            end
            
            %Initialize time structure
            %obj.time = struct('T',[],'dt',[],'nstep',[],'quiet',[]);
            %Determine the upper and lower bounds on the time interval
            obj.time.T = extractModelPropMultInput(TPWLtext,N,ind,'T',obj.time.T);
            %Determine the time step size
            obj.time.dt = extractModelPropMultInput(TPWLtext,N,ind,'dt',obj.time.dt);
            %Determine the number of time steps
            obj.time.nstep = extractModelPropMultInput(TPWLtext,N,ind,'nstep',obj.time.nstep);
            %Determine whether or not to time progress
            obj.time.quiet = extractModelPropMultInput(TPWLtext,N,ind,'timeQuiet',obj.time.quiet); 
            
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
            
            %Determine the requested TPWL model(s) size
            obj.nY = extractModelPropMultInput(TPWLtext,N,ind,'nY',obj.nY);
            %Determine the number of linearization points for the TPWL
            %model(s)
            obj.nLinPt = extractModelPropMultInput(TPWLtext,N,ind,'nLinPt',obj.nLinPt);
            %Determine maximum number of linearization points for the TPWL
            %model(s).  This will only be used in the TrajDist and ResDist
            %algorithms.
            obj.nLinPtMax = extractModelPropMultInput(TPWLtext,N,ind,'nLinPtMax',obj.nLinPtMax);
            %Determine the linearization point selection algorithm
            obj.LinPtAlgo = extractModelPropMultInput(TPWLtext,N,ind,'LinPtAlgo',obj.LinPtAlgo);
            %Determine the tunable parameter for the TPWL model(s)
            obj.LinPtAlgoParam = extractModelPropMultInput(TPWLtext,N,ind,'LinPtAlgoParam',obj.LinPtAlgoParam);
            %Determine the upper and lower bounds on the time interval where it is
            %valid to choose linearization points
            obj.LinPtInt = extractModelPropMultInput(TPWLtext,N,ind,'LinPtInt',obj.LinPtInt);
            %Determine the algorithm to use for selecting a reduced basis
            obj.basis = extractModelPropMultInput(TPWLtext,N,ind,'basis',obj.basis);
            %Determine the snapshot collection to use when constructing POD
            %reduced basis (only used if basis = 'pod')
            obj.snapcoll = extractModelPropMultInput(TPWLtext,N,ind,'snapcoll',obj.snapcoll);
            %Determine how to estimate the next reduced state vector in the
            %time integration scheme (either 'prev' or 'linextrap')
            obj.estNextSv = extractModelPropMultInput(TPWLtext,N,ind,'estNextSv',obj.estNextSv);
            %Establish the parameters beta and epsilon from Ch. 3 for
            %Reweinski thesis.  Recommended values are defaults.
            obj.beta = extractModelPropMultInput(TPWLtext,N,ind,'beta',obj.beta);
            obj.epsilon = extractModelPropMultInput(TPWLtext,N,ind,'epsilon',obj.epsilon);
            %Determine whether or not to time progress
            
%             %Determine the requested TPWL model(s) size
%             nYvec  = eval(determinePropText(TPWLtext,'nY',TPWLdelim));
%             if size(nYvec,1) == N
%                 obj.nY = nYvec(ind);
%             elseif size(nYvec,1) == 1
%                 obj.nY = nYvec;
%             elseif isempty(nYvec)
%                 obj.nY = [];
%             else
%                 error('In the TPWL field, nY must be either have 1 entry or the same number of entries as id');
%             end
%             
%             %Determine the number of linearization points for the TPWL model(s)
%             nLinPtVec = eval(determinePropText(TPWLtext,'nLinPt',TPWLdelim));
%             if size(nLinPtVec,1) == N
%                 obj.nLinPt = nLinPtVec(ind);
%             elseif size(nLinPtVec,1) == 1
%                 obj.nLinPt = nLinPtVec;
%             elseif isempty(nLinPtVec)
%                 obj.nLinPt = [];
%             else
%                 error('In the TPWL field, nLinPt must be either have 1 entry or the same number of entries as id');
%             end
%             
%             %Determine maximum number of linearization points for the TPWL
%             %model(s).  This will only be used in the TrajDist and ResDist
%             %algorithms.
%             nLinPtMaxVec = eval(determinePropText(TPWLtext,'nLinPtMax',TPWLdelim));
%             if size(nLinPtMaxVec,1) == N
%                 obj.nLinPtMax = nLinPtMaxVec(ind);
%             elseif size(nLinPtMaxVec,1) == 1
%                 obj.nLinPtMax = nLinPtMaxVec;
%             elseif isempty(nLinPtMaxVec)
%                 obj.nLinPtMax = [];
%             else
%                 error('In the TPWL field, nLinPtMax must be either have 1 entry or the same number of entries as id');
%             end
%             
%             %Determine the linearization point selection algorithm
%             LinPtAlgoVec = eval(determinePropText(TPWLtext,'LinPtAlgo',TPWLdelim));
%             if N > 1 && size(LinPtAlgoVec,1) == N
%                 obj.LinPtAlgo = LinPtAlgoVec{ind};
%             elseif size(LinPtAlgoVec,1) == 1
%                 if iscell(LinPtAlgoVec)
%                     obj.LinPtAlgo = LinPtAlgoVec{1};
%                 else
%                     obj.LinPtAlgo = LinPtAlgoVec;
%                 end
%             elseif isempty(LinPtAlgoVec)
%                 obj.LinPtAlgo = [];
%             else
%                 error('In the TPWL field, LinPtAlgo must be either have 1 entry or the same number of entries as id');
%             end
%             
%             %Determine the tunable parameter for the TPWL model(s)
%             LinPtAlgoParamVec = eval(determinePropText(TPWLtext,'LinPtAlgoParam',TPWLdelim));
%             if size(LinPtAlgoParamVec,1) == N
%                 obj.LinPtAlgoParam = LinPtAlgoParamVec(ind);
%             elseif size(LinPtAlgoParamVec,1) == 1
%                 obj.LinPtAlgoParam = LinPtAlgoParamVec;
%             elseif isempty(LinPtAlgoParamVec)
%                 obj.LinPtAlgoParam = [];
%             else
%                 error('In the TPWL field, LinPtAlgoParam must be either have 1 entry or the same number of entries as id');
%             end
%             
%             %Determine the upper and lower bounds on the time interval where it is
%             %valid to choose linearization points
%             LinPtIntVec = eval(determinePropText(TPWLtext,'LinPtInt',TPWLdelim));
%             if size(LinPtIntVec,1) == N
%                 obj.LinPtInt = LinPtIntVec(ind,:);
%             elseif size(LinPtIntVec,1) == 1
%                 obj.LinPtInt = LinPtIntVec;
%             elseif isempty(LinPtIntVec)
%                 obj.LinPtInt = [];
%             else
%                 error('In the TPWL field, LinPtInt must be either have 1 entry or the same number of entries as id');
%             end
%             
%             %Determine the algorithm to use for selecting a reduced basis
%             basisvec = eval(determinePropText(TPWLtext,'basis',TPWLdelim));
%             if N > 1 && size(basisvec,1) == N
%                 obj.basis = basisvec{ind};
%             elseif size(basisvec,1) == 1
%                 if iscell(basisvec)
%                     obj.basis = basisvec{1};
%                 else
%                     obj.basis = basisvec;
%                 end
%             elseif isempty(basisvec)
%                 obj.basis = [];
%             else
%                 error('In the TPWL field, basis must be either have 1 entry or the same number of entries as id');
%             end
% 
%             %Determine the snapshot collection to use when constructing POD
%             %reduced basis (only used if basis = 'pod')
%             snapcollvec = eval(determinePropText(TPWLtext,'snapcoll',TPWLdelim));
%             if N > 1 && size(snapcollvec,1) == N
%                 obj.snapcoll = snapcollvec{ind};
%             elseif size(snapcollvec,1) == 1
%                 if iscell(snapcollvec)
%                     obj.snapcoll = snapcollvec{1};
%                 else
%                     obj.snapcoll = snapcollvec;
%                 end
%             elseif isempty(snapcollvec)
%                 obj.snapcoll = [];
%             else
%                 error('In the TPWL field, snapcoll must be either have 1 entry or the same number of entries as id');
%             end
%             
%             %Determine how to estimate the next reduced state vector in the
%             %time integration scheme (either 'prev' or 'linextrap')
%             estNextSvVec = eval(determinePropText(TPWLtext,'estNextSv',TPWLdelim));
%             if N > 1 && size(estNextSvVec,1) == N
%                 obj.estNextSv = estNextSvVec(ind,:);
%             elseif size(estNextSvVec,1) == 1
%                 if iscell(estNextSvVec)
%                     obj.estNextSv = estNextSvVec{1};
%                 else
%                     obj.estNextSv = estNextSvVec;
%                 end
%             elseif isempty(estNextSvVec)
%                 obj.estNextSv = [];
%             else
%                 error('In the TPWL field, LinPtInt must be either have 1 entry or the same number of entries as id');
%             end
%             
%             %Establish the parameters beta and epsilon from Ch. 3 for
%             %Reweinski thesis.  Recommended values are defaults.
%             betavec = eval(determinePropText(TPWLtext,'beta',TPWLdelim));
%             if size(betavec,1) == N && size(betavec,1) > 1
%                 obj.beta = betavec(ind);
%             elseif size(betavec,1) == 1
%                 obj.beta = betavec;
%             elseif isempty(betavec)
%                 obj.beta = [];
%             else
%                 error('In the TPWL field, beta must be either have 1 entry or the same number of entries as id');
%             end
%             
%             epsilonvec = eval(determinePropText(TPWLtext,'epsilon',TPWLdelim));
%             if size(epsilonvec,1) == N && size(epsilonvec,1) > 1
%                 obj.epsilon = epsilonvec(ind);
%             elseif size(epsilonvec,1) == 1
%                 obj.epsilon = epsilonvec;
%             elseif isempty(epsilonvec)
%                 obj.epsilon = [];
%             else
%                 error('In the TPWL field, epsilon must be either have 1 entry or the same number of entries as id');
%             end
%             
%             %Determine whether or not to time progress
%             timequietvec = eval(determinePropText(TPWLtext,'timeQuiet',TPWLdelim));
%             if size(timequietvec,1) == N && size(timequietvec,1) > 1
%                 obj.time.quiet = timequietvec(ind);
%             elseif size(timequietvec,1) == 1
%                 obj.time.quiet = timequietvec;
%             elseif isempty(timequietvec)
%                 obj.time.quiet = [];
%             else
%                 error('In the TPWL field, timeQuiet must be either have 1 entry or the same number of entries as id');
%             end
        end
        
        function  [] = createTPWL(obj,X)
            %This method uses TPWL to reduce the FOM of the particular
            %problem.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - TPWL object
            %X       - matrix (nDOF_FOM x nstep+!) containing the solution
            %          of the FOM at each time step (stored columnwise).
            %          i.e. this is the training input.  The first column
            %          is the initial condition.
            %
            %Outputs:
            %-------
            %i_tstep  - vector containing the time step number that
            %            correspond to selected linearization points
            %t_vect    - vector containing the times where the reduced
            %            solution was computed during the special TPWL time
            %            stepping algorithm (in general, not equal to the
            %            vector of times time.T(1):time.dt:time.T(2))
            %--------------------------------------------------------------
            
            %Determine ndofs in FOM
            obj.FOMsize = obj.prob.config.ndof;
            N = obj.FOMsize;
            
            %Compute linearization points used specified algorithm
            [xi,R_I,Jf_I,~,s,obj.LinPt] = computeLinPoint(obj,X,obj.LinPtAlgo); 
            %When I include generalization to systems to the form:
            %d(g(U))/dt = f(U) + h(x) + B*u(t), the above ~ will be
            %replaced by Jg_I and the two lines of code below will be
            %removed.
            Jg_I = cell(s,1);
            [Jg_I{:}] = deal(eye(N));
            
            %Determine the actual number of linearization points used
            obj.nLinPt = s;
            
            %Compute the reduced basis via the model reduction method
            %specified in the input file (either POD, Krylov Moment
            %Matching, or balanced truncation)
            if obj.nY == N
                %If the reduced dimension is equal to the full dimension,
                %the reduced basis is simply the identity
                obj.phi = eye(N);
            else                
                obj.phi = computeTPWLbasis(obj,R_I,Jf_I,Jg_I,xi,N,X,obj.basis);
            end
            
            %Determine actual dimension of reduced model (after
            %linearization)...note, this may be different than the
            %specified dimension due to singularity issues when computing
            %reduced basis
            obj.nY = size(obj.phi,2);
            
            %Initialize a cell array where each entry will store the
            %reduced dynamics matrix (i.e. A matrix) at a linearization
            %point
            obj.A_ri = cell(s,1);
            
            %Initialize matrix to store projected "effective input" of LDS
            %at each linearization point
            obj.gamma = zeros(obj.nY,s);
            
            for ind = 1:s %Actual mathematical loop is from 0 to s-1, but we loop from
                %1 to s to adhere to MATLAB indexing
                %Project stored jacobians (the jacobian of the nonlinear
                %function serves as the dynamics matrix)
                obj.A_ri{ind,1} = obj.phi'*Jf_I{ind,1}*obj.phi;
                
                %Project RHS of the LDS equation that serves as an input matrix
                obj.gamma(:,ind) = obj.phi'*(R_I{ind,1} - Jf_I{ind,1}*xi(:,ind)); %-probobj.B*feval(probobj.config.inFunc,obj.time.T(1,1) + obj.time.dt*obj.LinPt(ind))
            end
            
            %Project the full-order input matrix to get the reduced input matrix
            obj.B_hat = obj.phi'*obj.prob.B;
            
            %Project the full-order state vector at linearization points
            %and the initial state
            obj.zi = obj.phi'*xi;
            obj.z0 = obj.phi'*obj.prob.ic;
        end
        
        function  [xi,R_I,Jf_I,Jg_I,s,i_tstep] = computeLinPoint(obj,X,flag)
            %This function computes the linearization points for the TPWL
            %model using the linearization point selection algorithm
            %specified in flag
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - TPWL object
            %X       - ndof x nstep state matrix, where the jth column
            %          corresponds to the solution of the FOM at time step
            %          j
            %flag    - string indicating the linearization point algorithm
            %          to use (valid values: 'All', 'TrajCurve',
            %          'TrajDist', 'ResDist')
            %
            %Outputs:
            %--------
            %xi      - ndof x s state matrix, where the jth column
            %          corresponds to the solution of the FOM at the jth
            %          linearization point
            %R_I     - s x 1 cell array the ith cell contains an nDofFOM x
            %          1 vector containing the residual at the ith
            %          linearization point
            %Jf_I     - s x 1 cell array the ith cell contains an nDofFOM x
            %          1 vector containing the jacobian (nl term) at the
            %          ith linearization point
            %Jg_I     - s x 1 cell array the ith cell contains an nDofFOM x
            %          1 vector containing the jacobian (nl time
            %          derivative)at the ith linearization point
            %s       - scalar specifying the number of linearization points
            %          selected
            %i_tstep - 1 x s array of indices indicating the time steps
            %          where linearization points were selected
            %--------------------------------------------------------------
            
            %Remove initial condition from state vector matrix
            X(:,1) = [];
            
            %Set up full time vector
            tvectF = obj.time.T(1):obj.time.dt:obj.time.T(2);

            %Set up a boolean vector that indicates whether a particular entry in
            %tvectF is out of bounds (1) or in bounds (0)
            OutOfBndsTbool = (tvectF(2:end) > obj.LinPtInt(2) | tvectF(2:end) < obj.LinPtInt(1));
            %Set up the time vector of only valid linearization point times
            tvectV = tvectF(~OutOfBndsTbool);
            %Determine the indices of the Full time vector that are inside the valid
            %linearization point range
            tvectVind = find(~OutOfBndsTbool);
            %Set the maximum number of time steps based on the
            %linearization point interval and time step size
            linptMAX = sum(~OutOfBndsTbool);
            %Remove state vectors that are not in the appropriate
            %time interval.
            X = X(:,tvectVind);
            
            if obj.nLinPt == linptMAX
                flag = 'All';
            elseif obj.nLinPt > linptMAX
                fprintf(['Warning...Requested too many linearization points for the interval specified.',...
                    'Will use all linearization points in specified interval (and the initial condition)']);
                flag = 'All';
            end
            
            switch lower(flag)
                case 'all' %Use all time steps (in valid lin. pt. selection interval)
                    %as linearization points
                    %Set the number of linearization points used
                    s = linptMAX+1;
                    %Initialize some of the outputs of appropriate type and size
                    xi = zeros(size(obj.prob.ic,1),s);
                    R_I = cell(s,1);
                    Jf_I = cell(s,1);
                    Jg_I = cell(s,1);
                    %Set output variable.  Always include the initial condition with
                    %this algorithm and then include all time steps in the appropriate
                    %range
                    i_tstep = [0,tvectVind];
                    
                    %Fill the residual and jacobian cell arrays with the IC information
                    [R_I{1,1},Jf_I{1,1}] = obj.prob.ResJacWithoutInput(obj.prob.ic,obj.time.T(1,1));
                    [~,Jg_I{1,1}] = obj.prob.ResJacTimeDer(obj.prob.ic,obj.time.T(1,1));
                    
                    %Set the first column of the state matrix as the IC
                    xi(:,1) = obj.prob.ic;
                    %Fill the remainder of the state matrix
                    xi(:,2:end) = X;
                    
                    %Loop over all valid time steps
                    for i = 1:linptMAX
                        %Record and store the residual and jacobian at the appropriate
                        %time steps
                        [R_I{i+1,1},Jf_I{i+1,1}] = obj.prob.ResJacWithoutInput(xi(:,i+1),tvectV(i));
                        [~,Jg_I{i+1,1}] = obj.prob.ResJacTimeDer(xi(:,i+1),tvectV(i));
                    end
                case 'uniform' %Use uniform distribution of time steps (in valid lin. pt. selection interval)
                               %as linearization points
                    %Set the number of linearization points used
                    s = obj.nLinPt+1;
                    %Initialize some of the outputs of appropriate type and size
                    xi = zeros(size(obj.prob.ic,1),s);
                    R_I = cell(s,1);
                    Jf_I = cell(s,1);
                    Jg_I = cell(s,1);
                    %Set output variable.  Always include the initial condition with
                    %this algorithm and then include all time steps in the appropriate
                    %range
                    ind = round(linspace(1,length(tvectVind),obj.nLinPt));
                    i_tstep = [0,tvectVind(ind)];
                    
                    %Fill the residual and jacobian cell arrays with the IC information
                    [R_I{1,1},Jf_I{1,1}] = obj.prob.ResJacWithoutInput(obj.prob.ic,obj.time.T(1,1));
                    [~,Jg_I{1,1}] = obj.prob.ResJacTimeDer(obj.prob.ic,obj.time.T(1,1));
                    
                    %Set the first column of the state matrix as the IC
                    xi(:,1) = obj.prob.ic;
                    %Fill the remainder of the state matrix
                    xi(:,2:end) = X(:,ind);
                    
                    %Loop over all valid time steps
                    for i = 1:obj.nLinPt
                        %Record and store the residual and jacobian at the appropriate
                        %time steps
                        [R_I{i+1,1},Jf_I{i+1,1}] = obj.prob.ResJacWithoutInput(xi(:,i+1),tvectV(i));
                        [~,Jg_I{i+1,1}] = obj.prob.ResJacTimeDer(xi(:,i+1),tvectV(i));
                    end
                    
                case 'trajcurve'
                    %Use a linearization point selection algorithm
                    %based on the curvature of the solution (in state
                    %space).  Developed by Zahr 2010.
                    
                    %Initialize some of the outputs of appropriate type and size
                    xi = zeros(size(obj.prob.ic,1),obj.nLinPt);
                    R_I = cell(obj.nLinPt,1);
                    Jf_I = cell(obj.nLinPt,1);
                    Jg_I = cell(obj.nLinPt,1);
                    i_tstep = zeros(1,obj.nLinPt);
                    %Fill the residual and jacobian cell arrays with the IC information
                    [R_I{1,1},Jf_I{1,1}] = obj.prob.ResJacWithoutInput(obj.prob.ic,obj.time.T(1,1));
                    [~,Jg_I{1,1}] = obj.prob.ResJacTimeDer(obj.prob.ic,obj.time.T(1,1));
                    %Set up the vector containing the pointwise curvature of the
                    %solution curve in state space.
                    rho = calcTimeCurv(X,obj.time.dt);
                    %Establish bins where a point will be selected from each bin
                    %bins = smartBins(s,linptMAX,1,rho);
                    bins = floor(linspace(1,linptMAX,obj.nLinPt))';
                    
                    %Set the first column of the state matrix as the IC
                    xi(:,1) = obj.prob.ic;
                    %Loop over s-1 linearization point (the first one has already been
                    %set as the IC)
                    for i = 1:obj.nLinPt-1
                        %Set the index of the current state solution
                        ind_x    = i+1;
                        %Set the indices of the current bin
                        ind_bins = [i,i+1];
                        
                        if i == obj.nLinPt-1
                            %If we are in the last bin, do not use the right edge as a
                            %bin bound because the curvature could not be computed at
                            %this point via finite differences (so just don't use it)
                            bin_rho = bins(ind_bins,1)-[0;1];
                        else
                            %Increment left bin bound to avoid including the same point
                            %twice (i.e. if this was not done, we might include a
                            %particular point at the right bound on iteration j and as
                            %the left bound on iteration j+1, which might cause this
                            %point to be selected twice)
                            bin_rho = bins(ind_bins,1)+[1;0];
                        end
                        
                        %Determine the index (or indices) of the maximum absolute
                        %curvature in the current bin
                        [~,I_rho]  = max(abs((rho(bin_rho(1):bin_rho(2),1))));
                        %If there are multiple indices with max curvature in this bin,
                        %take the one nearest the center
                        I_rho = I_rho(round(length(I_rho)/2));
                        %Compute the absolute index of the state matrix of where to
                        %grab the solution
                        I_X   = (bin_rho(1)-1) + I_rho;
                        %Set the appropriate vector in the state matrix at the current
                        %linearization point
                        xi(:,ind_x) = X(:,I_X);
                        %Record and store the residual and jacobian at the appropriate
                        %time steps
                        [R_I{i+1,1},Jf_I{i+1,1}] = obj.prob.ResJacWithoutInput(xi(:,i+1),tvectV(I_X));
                        [~,Jg_I{i+1,1}] = obj.prob.ResJacTimeDer(xi(:,i+1),tvectV(I_X));
                        %Set the indices of the linearization points selected
                        i_tstep(i+1) = tvectVind(I_X);
                    end
                case 'trajdist' %Use a linearization point selection algorithm
                    %based on the distance along the solution curve
                    %(in state space).  Developed by Rewienski.
                    %Initialize the outputs of appropriate type and size
                    xi = zeros(size(obj.prob.ic,1),obj.nLinPt);
                    R_I = cell(obj.nLinPt,1);
                    Jf_I = cell(obj.nLinPt,1);
                    Jg_I = cell(obj.nLinPt,1);
                    i_tstep = zeros(1,obj.nLinPt);
                    
                    %Set the first column of the state matrix as the IC
                    xi(:,1) = obj.prob.ic;
                    %Fill the residual and jacobian cell arrays with the IC information
                    [R_I{1,1},Jf_I{1,1}] = obj.prob.ResJacWithoutInput(obj.prob.ic,obj.time.T(1,1));
                    [~,Jg_I{1,1}] = obj.prob.ResJacTimeDer(obj.prob.ic,obj.time.T(1,1));
                    
                    %Initialize counter to ensure we do not select more linearization
                    %points than specified
                    i = 1; %This just ensures we don't go through the outer loop more
                    %than s-1 times
                    %Intialize counter to actually count the number of linearization
                    %points we select
                    inc = 1; %Already took initial condition as point
                    while i < obj.nLinPt %While we have not reached the specified number
                        %of linearization points
                        %Loop over all valid time steps
                        for k = 1:linptMAX
                            %Extract the appropriate state vector from the state matrix
                            x = X(:,k);
                            
                            %Set the denominator^2 of the distance function as the sum
                            %of the squares of the state vector at the selected
                            %linearization points
                            denom = sum(xi(:,1:i).^2);
                            
                            %Compute the min. relative distance from the current vector
                            %any of the state vectors at the selected linearization points
                            dist = sqrt(min(sum((repmat(x,1,i)-xi(:,1:i)).^2,1)./denom));
                            
                            %If we have moved del away from the closest vector (want to
                            %cover as much of the state space as possible), selected
                            %another linearization point
                            if dist >= obj.LinPtAlgoParam
                                %Record and store the residual and jacobian at the appropriate
                                %time steps
                                [R_I{i+1,1},Jf_I{i+1,1}] = obj.prob.ResJacWithoutInput(xi(:,i),tvectV(k+1));
                                [~,Jg_I{i+1,1}] = obj.prob.ResJacTimeDer(xi(:,i),tvectV(k+1));
                                inc = inc + 1; %Increment lin. pt. counter
                                break;
                            end
                        end
                        %Increment i
                        i = i+1;
                        %Store the appropriate snapshot in the output state matrix
                        %(solution at linearization points)
                        xi(:,i) = x; %need to use i+1 because i begins at 0, but MATLAB has no 0th entry
                        %Set the indices of the linearization points selected
                        i_tstep(i) = tvectVind(k);
                        %If we have taken the maximum number of points,
                        %stop
                        if inc >= obj.nLinPtMax
                            break;
                        end
                    end
                    %Update the number of linearization points
                    s = inc;
                    %Update all outputs
                    i_tstep(:,s+1:end) = [];
                    xi(:,s+1:end) = [];
                    temp1 = R_I;
                    temp2 = Jf_I;
                    temp3 = Jg_I;
                    R_I = cell(s,1); Jf_I = cell(s,1); Jg_I = cell(s,1);
                    [R_I{:,1}] = deal(temp1{1:s,1});
                    [Jf_I{:,1}] = deal(temp2{1:s,1});
                    [Jg_I{:,1}] = deal(temp3{1:s,1});
                case 'resdist'%Use a linearization point selection algorithm
                    %based on the distance along the residual.
                    %Initialize vector to store linearization points
                    linPT = [];
                    %Fill the residual and jacobian cell arrays with the IC information
                    [R1,A0] = obj.prob.ResJac(obj.prob.ic,obj.time.T(1,1));
                    %Determine the initial f vector
                    f0 = R1 - obj.prob.B*feval(obj.prob.config.inFunc,obj.time.T(1,1));
                    %Initialize a cell array to store information regarding the
                    %residual distance
                    g = cell(linptMAX,1);
                    %Loop over all valid time steps
                    for i = 1:linptMAX
                        %Determine the residual at the current time step
                        [fi,~] = obj.prob.ResJac(X(:,i),tvectV(i));
                        %Adjust fi vector
                        fi = fi - obj.prob.B*feval(obj.prob.config.inFunc,tvectV(i));
                        %Store the quantity that will be used to quantify the residual
                        %distance
                        g{i,1} = fi - f0 - A0*(X(:,i) - obj.prob.ic);
                        
                        %Determine the number of linearization points we have at this
                        %point
                        k = size(linPT,2);
                        %Initialize out distance variable to infinity (this ensures
                        %that the the first COMPUTED distance value will be taken on
                        %the first iteration in the loop below)
                        DELTA = inf;
                        
                        %Loop over all current linearization points
                        for j = 1:k
                            %Compute the distance between the current "residual" and
                            %all "residual" at linearization points
                            DELTAp = norm(g{i,1} - g{linPT(j,1),1});
                            %The distance we are interested in is the smallest of these
                            %distances (see next comment for the reason for this)
                            DELTA = min(DELTA,DELTAp);
                        end
                        
                        %If we have moved del away from the CLOSEST vector (want to
                        %cover as much of the space as possible), select another
                        %linearization point
                        if DELTA > obj.LinPtAlgoParam
                            %Store the linearization point
                            linPT = [linPT; i]; %#ok<AGROW>
                            if length(linPT) > obj.nLinPtMax
                                %If we have taken too many linearization points, break
                                %the loop
                                break;
                            end
                        end
                    end
                    %Store the solution vector at all selected linearization points
                    xi = [obj.prob.ic,X(:,linPT)];
                    %Determine the number of linearization points selected by the
                    %algorithm
                    s = size(xi,2);
                    %Initialize cell arrays of appropriate size to store the residual
                    %and jacobian at each linearization point
                    R_I = cell(s,1);
                    Jf_I = cell(s,1);
                    Jg_I = cell(s,1);
                    %Store the residual and jacobian corresponding to the initial
                    %condition
                    [R_I{1,1},Jf_I{1,1}] =  ResJacWithoutInput(obj.prob.ic,obj.time.T(1,1));
                    [~,Jg_I{1,1}] = obj.prob.ResJacTimeDer(obj.prob.ic,obj.time.T(1,1));
                    
                    %Loop over all linearization points selected
                    for i = 1:s-1
                        %Record and store the residual and jacobian at the appropriate
                        %time steps
                        [R_I{i+1,1},Jf_I{i+1,1}] = obj.prob.ResJacWithoutInput(xi(:,i+1),tvectV(linPT(i,1)));
                        [~,Jg_I{i+1,1}] = obj.prob.ResJacTimeDer(xi(:,i+1),tvectV(linPT(i,1)));
                    end
                    %Set the indices of the linearization points selected
                    i_tstep = [0,linPT'];
            end
            clear s;
            s = length(i_tstep);
        end
        
        function  [V] = computeTPWLbasis(obj,R_I,Jf_I,Jg_I,xi,N,X,flag)
            %This function computes the reduced basis for the TPWL model
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj  - TPWL object
            %R_I  - nLinPt x 1 cell array the ith cell contains an nDofFOM
            %       x 1 vector containing the residual at the ith
            %       linearization point
            %Jf_I  - nLinPt x 1 cell array the ith cell contains an nDofFOM
            %       x 1 vector containing the jacobian at the ith
            %       linearization point
            %Jf_I  - nLinPt x 1 cell array the ith cell contains an nDofFOM
            %       x 1 vector containing the jacobian at the ith
            %       linearization point
            %xi   - matrix whose columns represent the full state vector at
            %       the linearization point
            %N    - scalar indicating the size of the FOM
            %i_t  -
            %X    -  matrix (nDOF_FOM x nstep) containing the solution of
            %        the FOM at each time step (stored columnwise). i.e.
            %        this is the training input
            %flag - string indicating whether to use krylov moment
            %       matching, Arnoldi algorithm ('krylov') or POD ('pod')
            %       to construct a reduced order basis
            %
            %Outputs:
            %--------
            %V    - reduced order basis
            %--------------------------------------------------------------
            
            switch lower(flag)
                case 'krylov'
                    %This is the algorithm from pg. 51 of Reweinski thesis
                    Vagg = zeros(N,(2*obj.nY+1)*obj.nLinPt);
                    for i = 0:obj.nLinPt-1
                        ind = i+1;
                        
                        A_hat = Jf_I{ind,1};
                        G_hat = Jg_I{ind,1};
                        B1_hat = obj.prob.B;
                        B2_hat = R_I{ind,1} - Jf_I{ind,1}*xi(:,ind);
                        
                        [V1,~,~] = krylov(A_hat\G_hat,A_hat\B1_hat,[],obj.nY,2);
                        [V2,~,~] = krylov(A_hat\G_hat,A_hat\B2_hat,[],obj.nY,2);
                        
                        Vtilde = [V1, V2, xi(:,ind)];
                        Vagg(:,i*(2*obj.nY+1)+1:(i+1)*(2*obj.nY+1)) = Vtilde;
                    end
                    
                    V = svdTRUNC(Vagg,1e-10);
                case 'pod'
                    snap = snapshotCollection(X,obj.snapcoll,[]);
                    V = pod(snap,obj.nY);
                case 'baltrunc'
                    error('Not yet supported');
                otherwise
                    error('basis is TPWL module must be krylov, pod, or baltrunc');
            end
        end
        
        function  [w] = calcWeights(obj,ze)
            %This function calculates the weights for the TPWL model
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj  - TPWL object
            %ze   - vector containing the current projected state vector
            %
            %Outputs:
            %--------
            %w    - nLinPt x 1 vector containing the weights of the
            %       linearized model at each linearization point
            %--------------------------------------------------------------
            
            %Compute the distance (in state space) between the current
            %(reduced) state and the (reduced) state at all linearization
            %points
            di = zeros(1,obj.nLinPt);
            for i = 0:obj.nLinPt-1
                ind = i+1;
                di(1,ind) = round_nearest(norm(ze-obj.zi(:,ind),2),1e-13);
            end
            
            %Find the closest linearization point to the current state.
            [m,ind] = min(di);
            
            %Compute weights (pg. 41 of Reweinski Thesis)
            if m == 0
                w = zeros(1,obj.nLinPt);
                w(ind) = 1;
            else
                w = exp(-obj.beta*di/m);
                w = w/sum(w);
            end
        end
        
        function  [] = executeModel(obj)
            %This function performs the time stepping for the TPWL model
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - TPWL object
            %
            %Outputs:
            %--------
            %There are no outputs. 
            %--------------------------------------------------------------
            
            tTPWLstart = tic;
            
            %Initialize storage for reduced state vectors in time stepping
            %algorithm
            z_prev = obj.z0; %reduced state vector at previous time step
            ze = obj.z0; %estimate for current reduced state vector
            z = zeros(obj.nY,obj.time.nstep);
            
            %Initialize counters
            %Time stepping starts at DT
            t = obj.time.dt;
            ind = 1;
            %Initial time increment is DT
            dt = obj.time.dt;
            dtp = obj.time.dt;
            
            while t <= obj.time.T(2)
                %Calculate weights for the current reduced state vector
                w_e = calcWeights(obj,ze);
                
                %Compute dynamics matrix of TPWL model (i.e. stitch
                %together the linearized models with appropriate weights)
                A_hat = zeros(size(obj.A_ri{1,1}));
                for p = 1:obj.nLinPt
                    A_hat = A_hat + obj.A_ri{p,1}*w_e(1,p);
                end
                %Set up equations and solve for the reduced state vector at
                %the current time step
                LHS = (eye(obj.nY) - dt*A_hat);
                RHS = dt*(obj.gamma*w_e') + obj.B_hat*feval(obj.prob.config.inFunc,t)*dt + z_prev;
                z(:,ind) = LHS\RHS;
                
                %Record the time at which the current solution was computed
                t_vect(1,ind) = round_nearest(t,1e-10);  %#ok<AGROW>
                
                
                if norm(z(:,ind) - ze)/norm(z(:,ind)) > obj.epsilon
                    %If our estimate for the current time vector is not
                    %close enough to our computed solution, set back and
                    %try again
                    if ind == 1 && sum(obj.z0 ~= 0) == 0
                        if dt == obj.time.dt
                            t = t - dt;
                            ind = ind - 1;
                            dtp = dt;
                            dt = dt/4;
                        end
                    else
                        t = t - dt;
                        ind = ind - 1;
                        dtp = dt;
                        dt = dt/2;
                    end
                else
                    %Otherwise, update dt and display iteration progress
                    dtp = dt;
                    dt = obj.time.dt;
                    if ~obj.time.quiet && rem(ceil(100*t/obj.time.T(2)),10) == 0 && rem(ceil(100*(t+dt)/obj.time.T(2)),10) ~= 0
                        fprintf('Integration of %3.0f%% of the Time Domain Complete\n',ceil(100*t/obj.time.T(2)))
                    end
                end
                t = round_nearest(t,1e-10);
                dtp = round_nearest(dtp,1e-12);
                dt = round_nearest(dt,1e-12);
                
                %Compute next estimate using either the previous state
                %vector or linear extrapolation (depending on user input)
                if strcmpi(obj.estNextSv,'prev')
                    if ind == 0
                        ze = obj.z0;
                    else
                        ze = z(:,ind);
                    end
                elseif strcmpi(obj.estNextSv,'linextrap')
                    if ind == 0
                        ze = obj.z0;
                    elseif ind == 1
                        ze = z(:,ind) + (z(:,ind) - obj.z0)*((t+dt) - t)/(t-(t-dtp));
                    elseif ind == 2
                        ze = z(:,ind) + (z(:,ind) - obj.z0)*((t+dt) - t)/(t-(t-dtp));
                    else
                        ze = z(:,ind) + (z(:,ind) - z(:,ind-1))*((t+dt) - t)/(t-(t-dtp));
                    end
                end
                
                %Update dt, z_prev, and ind
                t = t + dt;
                if ind == 0
                    z_prev = obj.z0;
                else
                    z_prev = z(:,ind);
                end
                ind = ind + 1;
            end
            
            %If we won't simulate every time step specified in the input
            %file, perform linear extrapolation in each entry of the
            %reduced state vector vs. time
            if length(t_vect) ~= obj.time.nstep
                %If we won't simulate every time step specified in the input
                %file, perform linear extrapolation in each entry of the
                %reduced state vector vs. time
                
                %Initialize reduced state vector matrix
                obj.sv = zeros(size(z,1),obj.time.nstep);
                for i = 1:n
                    tVec = obj.T(1):obj.dt:obj.T(2);
                    %Loop over every entry of the reduced state vector and
                    %interpolate and extrapolate at all points in the time
                    %vector defined in the input file
                    obj.sv(i,:) = interp1(t_vect,z(i,:),tVec,'linear','extrap');
                end
            else
                %Otherwise, just store the reduced state vectors
                obj.sv = z;
            end
            
            obj.ontime = toc(tTPWLstart);
        end
        
        function  [svF] = reconstructFullState(obj,~)
            %This function reconstructs the full state vector from the
            %reduced state vector matrix and the reduced basis (phi)
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - TPWL object
            %
            %Outputs:
            %--------
            %svF     - full state vector matrix 
            %--------------------------------------------------------------
            svF = [obj.prob.ic,obj.phi*obj.sv];
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
               
            all2clear4PP = {'z0,zi','A_ri','B_ri','gamma'};
            
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
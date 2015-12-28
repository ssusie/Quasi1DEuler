classdef locGNATother < handle
    
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
        nGreed; %number of basis vectors to use in greed node selection
        nBases;
        addInd;
        
        %Compute quantities
        ontime;
        LocBasisHist;
        icHat; %initial condition
        S; %Local POD precomputed
        f; %Local POD precomputed
        
        sv; %State vector (reduced) nY x (nstep+1) matrix
        A; %Online matrix A from Carlberg et. al. 2011
        B; %Online matrix B from Carlberg et. al. 2011
        avgIter = 0;
        sampleNodes; %Vector containing the sample node numbers
        sampleInd; %Vector containing the sample indice numbers
        %(will be the same as sampleNodes if there is only 1 unknown per node)
        
        irstart; %vector of indices pointing to locations in jrow to indicate the start of a new indice
        jrow; %vector of indices indicating the indices of the full state vector where the state vector needs to be evaluated
        jdiag; %boolean vector the same size as jrow indicating whether or not that component corresponds to a diagonal entry
        phiYhat; %the partial reduced order basis (i.e. phiY(jrow,:))
        cTimeIter; %current time iteration number
    end
    
    properties (Hidden=true, SetAccess = private, GetAccess = public)
        %Temporary properties
        partialULoc; %partial state vector at the current time step
        partialUprevLoc; %partial state vector at the previous time step
        JhatTemp; %matrix of size nI x number of necessary state vectors used to store the jacobian
        phiYhatInd; %vector that will take phiYhat into its true partial form, i.e. phiYhat(phiYhatInd,:) is operator
        %Phi_hat from Carlberg et. al 2010 (online algorithm)
        reconstJhatInd; %vector that is used take the partial Jacobian from vector form into matrix form (JhatTemp),
        %i.e. JhatTemp(reconstJhatInd) = Jhat will take the vector Jhat into the appropriate matrix JhatTemp
        uniqueJROWind;
        jdiagHat;
        cLocBasis; %Local basis currently being used
        UrGlob;
        phiYhatUnion; %This is the unique union of the rows of all phiYhat
        phi;
    end
    
    methods
        function [obj] = locGNATother(ROMfile,romobj,id,oldobj)
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
            
            % Construct a new object based on a deep copy of an old object
            % of this class by copying properties over.
            if nargin == 4 && isa(oldobj,'locGNAT')
                copyProperties(obj,oldobj);
                return;
            end
            
            %Extract the GNAT text from the rom file
            GNATtext = readInFile('GNAT',ROMfile,1);
            
            %Copy the time structure and reduced order from the ROM object
            obj.time = romobj.time;
            obj.nY = romobj.nY;
            obj.nBases = romobj.nBases;
            obj.S = romobj.S;
            obj.f = romobj.f;
            obj.phi = romobj.phi;
            
            %Determine the GNAT properties based on the text in the ROM
            %section of the input file
            determineGNAT(obj,GNATtext,id);
        end
        
        function  [] = determineGNAT(obj,GNATtext,id)
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
            
            %Determine if any properties are specified more than once
            determineDuplicate(GNATtext);
            %Find all line break indicators in GNATtext
            GNATdelim = findstr(GNATtext,'//');
            
            %Determine and store configuration id
            idvec = eval(determinePropText(GNATtext,'id',GNATdelim));
            idvec = idvec(:);
            N = size(idvec,1);
            if N > 1
                obj.id = idvec(id);
            elseif N == 1
                obj.id = idvec;
            else
                obj.id = [];
            end
            
            %Determine the whether or not to add "free" indices.
            addindvec  = eval(determinePropText(GNATtext,'addInd',GNATdelim));
            if size(addindvec,1) == N && size(addindvec,1) > 1
                obj.addInd = addindvec{id,1};
            elseif size(addindvec,1) == 1
                obj.addInd = addindvec{1};
            elseif isempty(addindvec)
                obj.addInd = [];
            else
                error('In the GNAT field, addInd must be either have 1 entry or the same number of entries as id');
            end
            
            %Determine the number of entries of res/jac to use (may be a
            %vector for GNAT comparison plot)
            nIvec = eval(determinePropText(GNATtext,'nI',GNATdelim));
            if size(nIvec,1) == N && size(nIvec,1) > 1
                obj.nI = nIvec(id);
            elseif size(nIvec,1) == 1
                obj.nI = nIvec;
            elseif isempty(nIvec)
                obj.nI = [];
            else
                error('In the GNAT field, nI must be either have 1 entry or the same number of entries as id');
            end
            
            %Determine the size of subspace residual is constrained to lie in (may be a
            %vector for GNAT comparison plot)
            nRvec = eval(determinePropText(GNATtext,'nR',GNATdelim));
            if size(nRvec,1) == N && size(nRvec,1) > 1
                obj.nR = nRvec(id);
            elseif size(nRvec,1) == 1
                obj.nR = nRvec;
            elseif isempty(nRvec)
                obj.nR = [];
            else
                error('In the GNAT field, nR must be either have 1 entry or the same number of entries as id');
            end
            %Determine the size of subspace jacobian is constrained to lie in (may be a
            %vector for GNAT comparison plot)
            nJvec = eval(determinePropText(GNATtext,'nJ',GNATdelim));
            if size(nJvec,1) == N && size(nJvec,1) > 1
                obj.nJ = nJvec(id);
            elseif size(nJvec,1) == 1
                obj.nJ = nJvec;
            elseif isempty(nJvec)
                obj.nJ = [];
            else
                error('In the GNAT field, nJ must be either have 1 entry or the same number of entries as id');
            end
            
            %Determine the number of sample nodes to use and the number of
            %basis vectors to use in greedy selection of nodes
            nSampleVec = eval(determinePropText(GNATtext,'nSample',GNATdelim));
            if size(nSampleVec,1) == N && size(nSampleVec,1) > 1
                obj.nSample = nSampleVec(id);
            elseif size(nSampleVec,1) == 1
                obj.nSample = nSampleVec;
            elseif isempty(nSampleVec)
                obj.nSample = [];
            else
                error('In the GNAT field, nSample must be either have 1 entry or the same number of entries as id');
            end
            
            nGreedVec  = eval(determinePropText(GNATtext,'nGreed',GNATdelim));
            if size(nGreedVec,1) == N && size(nGreedVec,1) > 1
                obj.nGreed = nGreedVec(id);
            elseif size(nGreedVec,1) == 1
                obj.nGreed = nGreedVec;
            elseif isempty(nGreedVec)
                obj.nGreed = [];
            else
                error('In the GNAT field, nGreed must be either have 1 entry or the same number of entries as id');
            end
            
            %Determine the maximum number(s) of newton iterations for the GNAT(s)
            maxItervec  = eval(determinePropText(GNATtext,'maxIter',GNATdelim));
            if size(maxItervec,1) == N && size(maxItervec,1) > 1
                obj.newt.maxIter = maxItervec(id);
            elseif size(maxItervec,1) == 1
                obj.newt.maxIter = maxItervec;
            elseif isempty(maxItervec)
                obj.newt.maxIter = [];
            else
                error('In the GNAT field, maxIter must be either have 1 entry or the same number of entries as id');
            end
            
            %Determine the absolute tolerance(s) for newton convergence for the GNAT(s)
            epsvec  = eval(determinePropText(GNATtext,'eps',GNATdelim));
            if size(epsvec,1) == N && size(epsvec,1) > 1
                obj.newt.eps = epsvec(id);
            elseif size(epsvec,1) == 1
                obj.newt.eps = epsvec;
            elseif isempty(epsvec)
                obj.newt.eps = [];
            else
                error('In the GNAT field, eps must be either have 1 entry or the same number of entries as id');
            end
            
            %Determine whether or not to print Newton warnings
            newtquietvec = eval(determinePropText(GNATtext,'newtQuiet',GNATdelim));
            if size(newtquietvec,1) == N && size(newtquietvec,1) > 1
                obj.newt.quiet = newtquietvec(id);
            elseif size(newtquietvec,1) == 1
                obj.newt.quiet = newtquietvec;
            elseif isempty(newtquietvec)
                obj.newt.quiet = [];
            else
                error('In the GNAT field, newtQuiet must be either have 1 entry or the same number of entries as id');
            end
            
            %Determine whether or not to time progress
            timequietvec = eval(determinePropText(GNATtext,'timeQuiet',GNATdelim));
            if size(timequietvec,1) == N && size(timequietvec,1) > 1
                obj.time.quiet = timequietvec(id);
            elseif size(timequietvec,1) == 1
                obj.time.quiet = timequietvec;
            elseif isempty(timequietvec)
                obj.time.quiet = [];
            else
                error('In the GNAT field, timeQuiet must be either have 1 entry or the same number of entries as id');
            end
        end
        
        function  [] = createGNAT(obj,probobj,phiR,phiJ,phiY)
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
            
            %Determine the union of the masks of all bases
            jrowUnion = [];
            for kk = 1:obj.nBases
                jrowUnion = [jrowUnion; unique(obj.jrow{kk,1})'];
            end
            jrowUnion = unique(jrowUnion);
            
            %Determine the location in the mask union for each basis
            obj.phiYhatInd = cell(obj.nBases,1);
            for kk = 1:obj.nBases
                temp = unique(obj.jrow{kk,1})';
                [~,obj.phiYhatInd{kk,1},~] = intersect(jrowUnion,temp);
            end
            
            %Store handle to problem in instance
            for kk = 1:obj.nBases
                obj.icHat{kk,1} = probobj.ic(jrowUnion,1);
            end
            
            %Store the local basis over the union of all masks
            for kk = 1:obj.nBases
                obj.phiYhat{kk,1} = phiY(jrowUnion,:,kk);
            end
            
            %Determine nI based on the number of sample indices
            temp = length(obj.sampleInd{1,1});
            if obj.nBases == 1
                obj.nI = length(obj.sampleInd{1,1});
            else
                for kk = 2:obj.nBases
                    if temp ~= length(obj.sampleInd{kk,1}) && temp ~= 0 && ~isempty(obj.sampleInd{kk,1})
                        error('For Local GNAT, nI must be the same for all local bases');
                    end
                    temp = length(obj.sampleInd{kk,1});
                    
                    if ~isempty(obj.sampleInd{kk,1})
                        obj.nI = length(obj.sampleInd{kk,1});
                    end
                end
            end
            
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
            
            %Compute the online matrices A and B
            phiRtrunc = cell(obj.nBases,1);
            phiJtrunc = cell(obj.nBases,1);
            for kk = 1:obj.nBases
                if isempty(phiR{kk,1}) || isempty(phiJ{kk,1})
                    continue;
                end
                phiRtrunc{kk,1} = phiR{kk,1}(:,1:obj.nR);
                phiJtrunc{kk,1} = phiJ{kk,1}(:,1:obj.nJ);
            end
            computeOnlineMatrices(obj,phiRtrunc,phiJtrunc);
            
            %Initialize the state vector matrix (reduced coordinates)
            obj.sv = zeros(obj.nY,obj.time.nstep+1); %Reduced ic = 0
            
            %             %Set partialU and partialUprev to the initial condition
            %             %(reduced)
            %             setUpartial2ic(obj,probobj);
            
            %Setup UrGlob for local GNAT computations
            obj.UrGlob = zeros(obj.nY,obj.nBases);
            
            %Determine phiYhatInd.  Let X be a full state vector, then
            %Xpartial = X(phiYhatInd,:), where Xpartial has no repeats.
            obj.reconstJhatInd = cell(obj.nBases,1);
%             obj.phiYhatInd = cell(obj.nBases,1);
            obj.JhatTemp = cell(obj.nBases,1);
            
            obj.reconstJhatInd = cell(obj.nBases,1);
            for kk = 1:obj.nBases
                %Determine phiYhatInd.  Let X be a full state vector, then
                %Xpartial = X(phiYhatInd,:), where Xpartial has no repeats.
                %[~,obj.phiYhatInd{kk,1},J] = unique(obj.jrow{kk,1});
                [~,~,J] = unique(obj.jrow{kk,1});
                for i = 1:length(obj.irstart{kk,1})-1
                    obj.reconstJhatInd{kk,1} = [obj.reconstJhatInd{kk,1}, i + obj.nI*([J(obj.irstart{kk,1}(i):obj.irstart{kk,1}(i+1)-1)]-1)];%):J(obj.irstart(i+1)-1)]-1)];
                end
                obj.uniqueJROWind{kk,1} = J;
                
                %Initialize (sparse) JhatTemp
                obj.JhatTemp{kk,1} = spalloc(obj.nI,length(unique(obj.jrow{kk,1})),length(obj.reconstJhatInd{kk,1}));
                
                %Set up the vector containing the indices of the partial
                %vectors that correspond to diagonal entries
                obj.jdiagHat{kk,1} = obj.uniqueJROWind{kk,1}(logical(obj.jdiag{kk,1}));
            end
            %Store handle to problem in instance
            obj.probGNAT = probobj.createCopy4GNATloc(obj);
        end
        
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
            
            obj.partialULoc = obj.icHat{obj.cLocBasis,1}(obj.phiYhatInd{obj.cLocBasis},:);
            obj.partialUprevLoc = obj.icHat{obj.cLocBasis,1}(obj.phiYhatInd{obj.cLocBasis},:);
            
            %obj.partialULoc = obj.ic(unique(obj.jrow{obj.cLocBasis,1}),1);
            %obj.partialUprevLoc = obj.ic(unique(obj.jrow{obj.cLocBasis,1}),1);
        end
        
        function  [] = computeOnlineMatrices(obj,phiR,phiJ)
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
            
            obj.A = cell(obj.nBases,1);
            obj.B = cell(obj.nBases,1);
            
            for kk = 1:obj.nBases
                if isempty(phiR{kk,1}) || isempty(phiJ{kk,1})
                    continue;
                end
                obj.A{kk,1} = pinv(phiJ{kk,1}(obj.sampleInd{kk,1},:));
                obj.B{kk,1} = phiJ{kk,1}'*phiR{kk,1}*pinv(phiR{kk,1}(obj.sampleInd{kk,1},:));
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
            
            obj.sampleNodes = cell(obj.nBases,1);
            obj.sampleInd = cell(obj.nBases,1);
            tempNodes = cell(obj.nBases,1);
            tempInd = cell(obj.nBases,1);
            
            %Make sure that we have not specified too few sample nodes.
            if obj.nSample*probobj.ndim < obj.nGreed
                nSampleOLD = obj.nSample;
                obj.nSample = ceil(obj.nGreed/probobj.ndim);
                fprintf('Warning: not enough sample nodes! Increasing number of sample nodes from %d to %d\n',nSampleOLD,obj.nSample);
                clear nSampleOLD;
            end
            
            P = min(obj.nSample,obj.nGreed);
            Q = floor(obj.nGreed/P)*ones(P,1);
            ind = [1:P] < mod(obj.nGreed,P);
            Q(ind) = Q(ind) + 1;
            s = ceil(obj.nGreed/obj.nSample);
            SS = floor(obj.nSample*s/obj.nGreed)*ones(P,1);
            ind = [1:P] < mod(obj.nSample,obj.nGreed);
            SS(ind) = SS(ind)+1;
            
            for kk = 1:obj.nBases
                QQ = 0;
                
                if isempty(phiR{kk,1}) || isempty(phiJ{kk,1})
                    continue;
                end
                %Lines 6 and 7 from Algorithm 4 of Carlberg's thesis
                R = phiR{kk,1}(:,1:max(Q));
                J = phiJ{kk,1}(:,1:max(Q));
                
                for p = 1:P
                    QQ = QQ + Q(p);
                    for s = 1:SS(p)
                        %Determine the set of indices not already in our sample mesh
                        validNodes = setdiff(1:length(probobj.mesh.node),tempNodes{kk,1}');
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
                        tempNodes{kk,1} = [tempNodes{kk,1};n];
                        tempInd{kk,1} = [tempInd{kk,1};probobj.node2ind(n)];
                    end
                    %Lines 13-16 of Algorithm 4 of Carlberg thesis
                    for q = 1:Q(p)
                        a = phiR{kk,1}(tempInd{kk,1},1:QQ)\phiR{kk,1}(tempInd{kk,1},QQ+q);
                        b = phiJ{kk,1}( tempInd{kk,1},1:QQ)\phiJ{kk,1}(tempInd{kk,1},QQ+q);
                        
                        R(:,q) = phiR{kk,1}(:,QQ+q) - phiR{kk,1}(:,1:QQ)*a;
                        J(:,q) = phiJ{kk,1}(:,QQ+q) - phiJ{kk,1}(:,1:QQ)*b;
                    end
                end
                
                %Sort the sample nodes and indices
                obj.sampleInd{kk,1} = unique(tempInd{kk,1});
                obj.sampleNodes{kk,1} = unique(tempNodes{kk,1});
            end
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
            
            obj.jdiag = cell(obj.nBases,1);
            obj.jrow = cell(obj.nBases,1);
            obj.irstart = cell(obj.nBases,1);
            
            for kk = 1:obj.nBases
                if isempty(obj.sampleInd{kk,1})
                    continue;
                end
                obj.irstart{kk,1} = 1;
                for i = 1:length(obj.sampleInd{kk,1})
                    %Find the nonzero entries in the appropriate row of Jstruct
                    temp = find(probobj.Jstruct(obj.sampleInd{kk,1}(i,1),:) ~= 0);
                    %Store these indices in jrow
                    obj.jrow{kk,1} = [obj.jrow{kk,1},temp];
                    %Determine which of these values appended onto jrow
                    %corresponds to the diagonal.
                    obj.jdiag{kk,1} = [obj.jdiag{kk,1},temp == obj.sampleInd{kk,1}(i)];
                    
                    %Store the index of the start of the next node
                    obj.irstart{kk,1} = [obj.irstart{kk,1}, length(obj.jrow{kk,1})+1];
                end
                %Make jdiag a column vector
                obj.jdiag{kk,1} = obj.jdiag{kk,1}(:);
            end
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
            
            obj.sampleNodes = cell(obj.nBases,1);
            obj.sampleInd   = cell(obj.nBases,1);
            
            for kk = 1:obj.nBases
                if isempty(phiR{kk,1}) || isempty(phiJ{kk,1})
                    continue;
                end
                % Algorithm 5 of Carlberg's 2010 Paper
                R = phiR{kk,1}(:,1);
                J = phiJ{kk,1}(:,1);
                
                nbarI = 0;
                m = 1;
                
                while nbarI < obj.nI
                    %Determine the set of indices not already in our sample
                    %indices
                    validNodes = setdiff(1:length(R),obj.sampleInd{kk,1}');
                    %Add up the contribution of the residual and
                    %jacobian for each valid index
                    temp = R(validNodes,1).^2 + J(validNodes,1).^2;
                    %Determine the node with the max. value
                    [~,ID] = max(temp);
                    n = validNodes(ID);
                    %Add the computed indices to the set
                    obj.sampleInd{kk,1} = [obj.sampleInd{kk,1};n];
                    K = checkResEvaluatedForFree(obj,probobj,n,nbarI,kk);
                    obj.sampleInd{kk,1} = unique([obj.sampleInd{kk,1};K(:)]);
                    
                    nbarI = nbarI + 1 + length(K);
                    m = m+1;
                    pR = min(m-1,obj.nR); pJ = min(m-1,obj.nJ);
                    %Lines 10-11 of Algorithm 5 of Carlberg paper 2010
                    a = phiR{kk,1}(obj.sampleInd{kk,1},1:pR)\phiR{kk,1}(obj.sampleInd{kk,1},m);
                    b = phiJ{kk,1}( obj.sampleInd{kk,1},1:pJ)\phiJ{kk,1}(obj.sampleInd{kk,1},m);
                    
                    R = phiR{kk,1}(:,m) - phiR{kk,1}(:,1:pR)*a;
                    J = phiJ{kk,1}(:,m) - phiJ{kk,1}(:,1:pJ)*b;
                end
                %Sort the sample nodes and indices
                obj.sampleInd{kk,1} = unique(obj.sampleInd{kk,1});
            end
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
                            Jref = obj.jrow{iBases,1};
                    end
                    IndexComp = setdiff(1:probobj.config.ndof,obj.sampleInd{iBases,1});
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
            %Set initial guess as last time step (note that sv is the
            %partial state). Initial guess for the sv = 0 for GNAT, but we
            %do need to update the partial state initial guess.
            if obj.cTimeIter > 1
                if (obj.cLocBasis==obj.LocBasisHist(obj.cTimeIter-1))
                    obj.partialUprevLoc = obj.partialULoc;
                else
                    U = obj.icHat{obj.cLocBasis,1};
                    for kk=1:obj.nBases
                        if ~isempty(obj.phiYhat{kk,1})
                            U = U + obj.phiYhat{kk,1}*obj.UrGlob(:,kk);
                        end
                    end
                    obj.partialUprevLoc = U;%(unique(obj.jrow{obj.cLocBasis}));
                    obj.partialULoc = obj.partialUprevLoc;
                    %                     obj.partialUprevLoc = U(unique(obj.jrow{obj.cLocBasis}));
                end
            end
            %Determine the residual and jacobian based on initial guess
            [Rhat,JVhat] = BackwardEulerNLFuncGNATloc(obj);
            %Update local operators
            [RAbar,E] = updateLocOperators(obj,Rhat,JVhat);
            %Use relative tolerance
            tol = obj.newt.eps*norm(E);
            
            for i_N = 1:obj.newt.maxIter %loop until the maximum newton iterations are reached
                %Update the state vector (reduced) with search direction
                %(unit step length)
                obj.sv(:,itnump1) = obj.sv(:,itnump1) - RAbar\E;
                obj.partialULoc = obj.partialUprevLoc + obj.phiYhat{obj.cLocBasis}(obj.phiYhatInd{obj.cLocBasis},:)*obj.sv(:,itnump1);
                %Compute residual and jacobian with updated vector for next iteration
                [Rhat,JVhat] = BackwardEulerNLFuncGNATloc(obj);
                %Solve for the search direction
                [RAbar,E] = updateLocOperators(obj,Rhat,JVhat);
                
                %Stop iterating once residual is small enough
                if (norm(E)<tol)
                    break;
                end
            end
            
            obj.UrGlob(:,obj.cLocBasis) = obj.UrGlob(:,obj.cLocBasis) + obj.sv(:,itnump1);
            
            %Update average number of newton iterations
            obj.avgIter = obj.avgIter + i_N/obj.time.nstep;
            
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
            obj.JhatTemp{obj.cLocBasis}(obj.reconstJhatInd{obj.cLocBasis,1}) = Jhat;
            JVhat = obj.JhatTemp{obj.cLocBasis}*obj.phiYhat{obj.cLocBasis,1}(obj.phiYhatInd{obj.cLocBasis},:);
            
        end
        
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
            
            nstep = obj.time.nstep; %Extract the number of time steps from the object
            
            obj.LocBasisHist = zeros(nstep,1);
            tic;
            for i_t = 1:nstep %Loop over each time step to determine all state vectors
                % Compute reference state distances to bases centers
                distCenters = computeDistanceToCenters(obj,obj.UrGlob);
                [~,I] = sort(distCenters,'ascend');
                obj.cLocBasis = I(1);
                obj.LocBasisHist(i_t,1) = I(1);
                
                %                 obj.cLocBasis = romobj.LocBasisHist(i_t,1);
                %                 obj.LocBasisHist(i_t,1) = romobj.LocBasisHist(i_t,1);
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
                
                %Use Newton's Method to solve nonlinear system of equations
                %and store: state vector, residual and jacobian snapshots,
                %and number of newton iterations required
                obj.NewtonRaphsonLocal;
            end
            obj.ontime = toc; %Record simulation time
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
            
            svF = zeros(size(romobj.phi,1),romobj.time.nstep+1);
            %First column is the initial condition.
            svF(:,1) = obj.prob.ic;
            for i = 2:obj.time.nstep+1
                %The newton method that was implemented for the GNAT
                %computes the current reduced state vector relative to the
                %previous state vector.  Therefore, we reconstructed the
                %full state by looping through the time steps and adding
                %the reconstructed state difference to the previous state
                %vector.
                svF(:,i) = svF(:,i-1) + romobj.phi(:,:,obj.LocBasisHist(i-1))*obj.sv(:,i);
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
                'partialU','partialUprev','JhatTemp','phiYhatInd',...
                'reconstJhatInd','uniqueJROWind','jdiagHat','UrGlob',...
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
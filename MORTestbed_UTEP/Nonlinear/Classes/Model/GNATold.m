classdef GNATold < handle
    
    properties (SetAccess = private, GetAccess = public)
        %Quantities read from input file
        id;
        nY;
        nR;
        nJ;
        nI;
        newt;
        time;
        nSample; %number of sample nodes to use
        nGreed; %number of basis vectors to use in greed node selection
        
        %Compute quantities
        sv; %State vector (reduced) nY x (nstep+1) matrix
        A; %Online matrix A from Carlberg et. al. 2010
        B; %Online matrix B from Carlberg et. al. 2010
        avgIter = 0;
        sampleNodes; %Vector containing the sample node numbers
        sampleInd; %Vector containing the sample indice numbers 
                   %(will be the same as sampleNodes if there is only 1 unknown per node)
        
        irstart; %vector of indices pointing to locations in jrow to indicate the start of a new indice
        jrow; %vector of indices indicating the indices of the full state vector where the state vector needs to be evaluated
        jdiag; %boolean vector the same size as jrow indicating whether or not that component corresponds to a diagonal entry
        phiYhat; %the partial reduced order basis (i.e. phiY(jrow,:))
    end
    
    properties (Hidden=true, SetAccess = private, GetAccess = public)
        %Temporary properties
        partialU; %partial state vector at the current time step
        partialUprev; %partial state vector at the previous time step
        JhatTemp; %matrix of size nI x number of necessary state vectors used to store the jacobian
        phiYhatInd; %vector that will take phiYhat into its true partial form, i.e. phiYhat(phiYhatInd,:) is operator 
                    %Phi_hat from Carlberg et. al 2010 (online algorithm)
        reconstJhatInd; %vector that is used take the partial Jacobian from vector form into matrix form (JhatTemp), 
                        %i.e. JhatTemp(reconstJhatInd) = Jhat will take the vector Jhat into the appropriate matrix JhatTemp
    end
    
    properties (SetAccess = public, GetAccess = public)
        cTimeIter; %current time iteration number
    end
    
    methods
        function [obj] = GNATold(ROMfile,romobj)
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
            %Outputs:
            %--------
            %obj     - instance of the GNAT object that was constructed
            %--------------------------------------------------------------
            
            %Extract the GNAT text from the rom file
            GNATtext = readInFile('GNAT',ROMfile,1);
            
            %Copy the time structure and reduced order from the ROM object
            obj.time = romobj.time;
            obj.nY = romobj.nY;
            
            %Determine the GNAT properties based on the text in the ROM
            %section of the input file 
            determineGNAT(obj,GNATtext);
        end
        
        function  [] = determineGNAT(obj,GNATtext)
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
            
            %Determine the size of subspace residual is constrained to lie in (may be a
            %vector for GNAT comparison plot)
            obj.nR = eval(determinePropText(GNATtext,'nR',GNATdelim));
            %Determine the size of subspace jacobian is constrained to lie in (may be a
            %vector for GNAT comparison plot)
            obj.nJ = eval(determinePropText(GNATtext,'nJ',GNATdelim));
            %Determine the number of sample nodes to use and the number of
            %basis vectors to use in greedy selection of nodes
            obj.nSample = eval(determinePropText(GNATtext,'nSample',GNATdelim));
            obj.nGreed  = eval(determinePropText(GNATtext,'nGreed',GNATdelim));
            %Determine the maximum number of newton iterations for the GNAT model
            obj.newt.maxIter = eval(determinePropText(GNATtext,'maxIter',GNATdelim));
            %Determine the absolute tolerance for newton convergence for the GNAT model
            obj.newt.eps = eval(determinePropText(GNATtext,'eps',GNATdelim));
        end
        
        function  [] = createGNAT(obj,probobj,phiR,phiJ,phiY)
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
            
            %Compute the sample nodes based on a greedy node selection
            %algorithm
            computeSampleNodes(obj,probobj,phiR,phiJ);
            
            %Compute additional nodes required in the reduced mesh in order
            %to compute the residual and jacobian at the sample nodes.
            determineAdditionalNodes(obj,probobj);
            
            %Extract the phiYhat from the reduced order basis
            obj.phiYhat = phiY(obj.jrow,:);
            
            %Determine nI based on the number of sample nodes being sure to
            %account for the fact that we may have selected a node where 1
            %dof is a dbc and others are not.
%             rmdof = 0;
%             for i = 1:obj.nSample 
%                 if isempty(probobj.mesh.dbc)
%                     %If there are no DBCs in the problem, flag should
%                     %always be empty
%                     flag = [];
%                 else
%                     flag = find(probobj.mesh.dbc(:,1) == obj.sampleNodes(i));
%                 end
%                 if ~isempty(flag)
%                     for j = 1:probobj.ndim
%                         rmdof = rmdof + double(probobj.mesh.dbc(flag,j+1));
%                     end
%                 end
%             end
%             obj.nI = obj.nSample*probobj.ndim - rmdof;
            obj.nI = length(obj.sampleInd);

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
            computeOnlineMatrices(obj,phiR(:,1:obj.nR),phiJ(:,1:obj.nJ));
            
            %Initialize the state vector matrix (reduced coordinates)
            obj.sv = zeros(obj.nY,obj.time.nstep+1); %Reduced ic = 0
            
            %Set partialU and partialUprev to the initial condition
            %(reduced)
            setUpartial2ic(obj,probobj);
            
            %Determine phiYhatInd.  Let X be a full state vector, then
            %Xpartial = X(phiYhatInd,:), where Xpartial has no repeats.
            [~,obj.phiYhatInd,J] = unique(obj.jrow);
            obj.reconstJhatInd = [];
            for i = 1:length(obj.irstart)-1
                obj.reconstJhatInd = [obj.reconstJhatInd, i + obj.nI*([J(obj.irstart(i):obj.irstart(i+1)-1)]-1)];%):J(obj.irstart(i+1)-1)]-1)];
            end
            
            %Initialize (sparse) JhatTemp
            obj.JhatTemp = spalloc(obj.nI,length(unique(obj.jrow)),length(obj.reconstJhatInd));
        end
        
        function  [] = setUpartial2ic(obj,probobj)
            %This function sets partialU and partialUprev to the ic.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - instance of GNAT class
            %probobj - Problem object
            %
            %Outputs:
            %--------
            %No outputs.  partialU and partialUprev are overwritten in GNAT
            %handle class. 
            %--------------------------------------------------------------
            
            obj.partialU = probobj.ic(obj.jrow,1);
            obj.partialUprev = probobj.ic(obj.jrow,1);
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
           
            obj.A = pinv(phiJ(obj.sampleInd,:));
            obj.B = phiJ'*phiR*pinv(phiR(obj.sampleInd,:));
        
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
        
        function  [] = NewtonRaphson(obj,probobj)
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
            %Set initial guess as last time step (note that sv is the
            %partial state). Initial guess for the sv = 0 for GNAT, but we
            %do need to update the partial state initial guess.
            obj.partialUprev = obj.partialU;
            %Determine the residual and jacobian based on initial guess
            [Rhat,JVhat] = BackwardEulerNLFuncGNAT(obj,probobj);
            %Update local operators
            [RAbar,E] = updateLocOperators(obj,Rhat,JVhat);
            %Use relative tolerance
            tol = obj.newt.eps*norm(E);

            for i_N = 1:obj.newt.maxIter %loop until the maximum newton iterations are reached
                %Update the state vector (reduced) with search direction
                %(unit step length)
                obj.sv(:,itnump1) = obj.sv(:,itnump1) - RAbar\E;
                obj.partialU = obj.partialUprev + obj.phiYhat*obj.sv(:,itnump1);
                %Compute residual and jacobian with updated vector for next iteration
                [Rhat,JVhat] = BackwardEulerNLFuncGNAT(obj,probobj);
                %Solve for the search direction
                [RAbar,E] = updateLocOperators(obj,Rhat,JVhat);
                
                %Stop iterating once residual is small enough
                if (norm(E)<tol)
                    break;
                end
            end
            
            %Update average number of newton iterations
            obj.avgIter = obj.avgIter + i_N/obj.time.nstep;
            
%             %If the maxiumum newton iterations are reached, warn the user
%             if (i_N==obj.newt.maxIter)
%                 disp(['*** Warning; Newton solver reached max number of iterations before convergence : res = ',...
%                     num2str(norm(E,2)/tol),'  eps = ',num2str(obj.newt.eps), '***']);
%             end
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
            obj.JhatTemp(obj.reconstJhatInd) = Jhat;
            JVhat = obj.JhatTemp*obj.phiYhat(obj.phiYhatInd,:);
                        
        end
        
        function  [ctime] = executeModel(obj,probobj)
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
            %ctime   - online simulation time of GNAT
            %--------------------------------------------------------------
            
            nstep = obj.time.nstep; %Extract the number of time steps from the object
            
            %Make sure partialU and partialUprev are properly set whenever this
            %function is called
            modelobj.setUpartial2ic(probobj);
            
            tic;
            for i_t = 1:nstep %Loop over each time step to determine all state vectors
                if rem(i_t,round(nstep*0.1)) == 0
                    %Generate output so user can see progress
                    fprintf('%4.0f of %4.0f Timesteps Complete (%2.0f%%)\n',i_t,nstep,100*i_t/nstep)
                end
                
                %Set current time iteration (so model can keep track)
                obj.cTimeIter = i_t;
                
                %Use Newton's Method to solve nonlinear system of equations
                %and store: state vector, residual and jacobian snapshots,
                %and number of newton iterations required
                NewtonRaphson(obj,probobj);
            end
            ctime = toc; %Record simulation time
        end
        
        function  [svF] = reconstructFullState(obj,probobj,romobj)
            %This function reconstructs the full state vector from the
            %reduced state vector matrix and the reduced basis (phi)
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - GNAT object
            %probobj - Problem object
            %romobj  - ROM object
            %
            %Outputs:
            %--------
            %svF     - full state vector matrix 
            %--------------------------------------------------------------
            
            svF = zeros(probobj.config.ndof,romobj.time.nstep+1);
            %First column is the initial condition.
            svF(:,1) = probobj.ic;
            for i = 2:obj.time.nstep+1
                %The newton method that was implemented for the GNAT
                %computes the current reduced state vector relative to the
                %previous state vector.  Therefore, we reconstructed the
                %full state by looping through the time steps and adding
                %the reconstructed state difference to the previous state
                %vector.
                svF(:,i) = svF(:,i-1) + romobj.phi*obj.sv(:,i);
            end
            
        end
    end
end
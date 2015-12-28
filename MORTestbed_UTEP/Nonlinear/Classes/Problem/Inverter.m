classdef Inverter 
    
    properties
        config;
        staticFlag = false;
        A;
        B
        Vdd;
        Jstruct;
        ic;
    end
    
    properties (Hidden=true)
        caseNum;
        vInd;
        numInd;
    end
    
    methods
        function  [obj] = Inverter(cfgobj,oldobj)
            %This is the constructor for the Inverter Chain Circuit attractor probleme.
            %If there is 1 input, the constructor will create the class
            %instance from the variables defined in this functino.  If
            %there are 2 inputs, the first input will be unused and the
            %output will be a deep copy of the oldobj.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %cfgobj - CONFIG object
            %oldobj - Inverter Chain Circuit object whose properties will be copied
            %         into the new object.
            %
            %Outputs:
            %--------
            %obj    - Inverter Chain Circuit object
            %--------------------------------------------------------------
            
            if nargin == 0
                %Keep everything at defaults if no input arguments
                return;
            end
            
            % Construct a new object based on a deep copy of an old object
            % of this class by copying properties over.
            if nargin == 2 && isa(oldobj,'Inverter')
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
            
             %Store config object
            obj.config = cfgobj;
            
            [obj.A,obj.Vdd] = feval(obj.config.altFile);

            
            %Setup Jacobian Structure
            obj.Jstruct = JacobianStructure(obj);
            
           %Determine input matrix from the input matrix function
            if isa(cfgobj.BFunc,'function_handle')
                obj.B = feval(cfgobj.BFunc);
            else
                obj.B = feval(cfgobj.BFunc);
            end
            
            
            
            
            %Determine ic term vector from the ic function and the
            %coordinate vector
            if isa(cfgobj.icFunc,'function_handle')
                obj.ic = cfgobj.icFunc();
            else
                obj.ic = feval(cfgobj.icFunc);
            end
            
        end
                
        function  [R,J] =  ResJac(obj,v,t)
        %This function calcuates the residual and jacobian of the Inverter Chain Circuit problem
        %------------------------------------------------------------------
        %Inputs:
        %-------
        %obj - Inverter Chain Circuit object
        %v   - N x 1 vector containing the "solution" (or guess or
        %      approximation) at the current time
        %t   - scalar indicating the time
        %
        %Outputs:
        %--------
        %R   - N x 1 vector containing the nonlinear residual of the
        %      solution v at the current time t
        %J   - N x N vector containing the nonlinear jacobian of the
        %      solution v at the current time t
        %------------------------------------------------------------------
        N = length(v);

        [g,dg] = Inverter.nonlinTerm([0;v(1:end-1,1)],obj.A,obj.Vdd);
        R = - v + g + obj.B*feval(obj.config.inFunc,t);
           
        J = spdiags([-ones(N,1),[dg(2:end);0]],[0,-1],N,N);
       % J
        end
        
        function  [R,J] =  ResJacGNAT(obj,v,t)
            %This function calculates the residual and jacobian of the
            %inverter chain circuit
            %------------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - Inverter object
            %U       - length(jrow) x 1 vector containing the "solution" (or
            %          guess or approximation) at the current time
            %t       - scalar indicating the time
            %
            %Outputs:
            %--------
            %R   - indN x 1 vector containing the nonlinear residual of the
            %      solution U at the current time t
            %J   - vector containing the nonlinear jacobian of the solution
            %      U at the current time t at the rows indN
            %--------------------------------------------------------------
            
            switch obj.caseNum
                case 1 % row 1 only
                     R = [obj.sigma*(v(2,1)-v(1,1))];
                 
            end
        end
        
        function  [R,J] = ResJacTimeDer(obj,U,~)
            R = U;
            J = eye(size(U));
        end
        
        function  [R,J] = ResJacGNATTimeDer(obj,U,~,~,~,jrow,jdiag)
            R = U(jrow(logical(jdiag)));
            J = jdiag;
        end
        
                
        function  [ind] = node2ind(obj,nodenum)
            %This function takes a node in the mesh and maps it to an index
            %in the full state vector for the Inverter Chain Circuit attractor
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - Inverter Chain Circuit object
            %nodenum - scalar specifying node numbers
            %
            %Outputs:
            %--------
            %ind     - vector specifying the corresponding indices in the
            %          full state vector
            %--------------------------------------------------------------
            
            ind = nodenum;

        end
        
        function  [JacStruct] = JacobianStructure(obj)
            %This function forms a boolean matrix whose true entries
            %correspond to nonzero entries in the jacobian for the
            %Inverter Chain Circuit attractor
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj       - Inverter Chain Circuit object
            %
            %Outputs:
            %--------
            %JacStruct - ndof x ndof boolean matrix whose true values
            %            correspond to nonzero entries in the jacobian
            %--------------------------------------------------------------
            
            
            N = obj.config.ndof;
            JacStruct = spdiags(ones(N,2),-1:0,N,N); 
            
        end
        
        function  [probGNAT] = createCopy4GNAT(obj,gnat)
            %This function makes a copy of the current object and only
            %copies the source term and input matrix evaluated at indN
            %rows.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - Inverter Chain Circuit object
            %gnat    - GNAT object
            %
            %Outputs:
            %--------
            %probGNAT - new "reduced" Inverter Chain Circuit object
            %--------------------------------------------------------------
            
            probGNAT = Inverter([],obj);
            
            probGNAT.config = [];
            probGNAT.Jstruct = [];
            
            probGNAT.numInd = length(gnat.sampleInd);
            if (length(gnat.sampleInd)==1)
                if gnat.sampleInd(1) == 1
                   probGNAT.caseNum = 1;
                elseif gnat.sampleInd(1) == 2
                   probGNAT.caseNum = 2;
                elseif gnat.sampleInd(1) == 3
                   probGNAT.caseNum = 3;
                end
            elseif (length(gnat.sampleInd)==2)
                if gnat.sampleInd(1) == 1 && gnat.sampleInd(2) == 2
                    probGNAT.caseNum = 4;
                elseif gnat.sampleInd(1) == 1 && gnat.sampleInd(2) == 3
                    probGNAT.caseNum = 5;
                elseif gnat.sampleInd(1) == 2 && gnat.sampleInd(2) == 3
                    probGNAT.caseNum = 6;
                end
             elseif (length(gnat.sampleInd)==3)
                 probGNAT.caseNum = 7;
            end
            probGNAT.vInd = gnat.irstart;   
            if length(obj.B) > 1
                probGNAT.B = obj.B(gnat.sampleInd,1);
            end    

        end
        
        function  [probGNAT] = createCopy4GNATloc(obj,gnat)
            %This function makes a copy of the current object and only
            %copies the source term and input matrix evaluated at indN
            %rows.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - Inverter Chain Circuit object
            %gnat    - locGNAT object
            %
            %Outputs:
            %--------
            %probGNAT - new "reduced" Inverter Chain Circuit object
            %--------------------------------------------------------------
            
            for kk = 1:gnat.nBases
                probGNAT(kk) = Inverter([],obj);
                
                probGNAT(kk).config = [];
                probGNAT(kk).Jstruct = [];
                
                probGNAT(kk).numInd = length(gnat.sampleInd{kk,1});
                if (length(gnat.sampleInd{kk,1})==1)
                    if gnat.sampleInd{kk,1}(1) == 1
                        probGNAT(kk).caseNum = 1;
                    elseif gnat.sampleInd{kk,1}(1) == 2
                        probGNAT(kk).caseNum = 2;
                    elseif gnat.sampleInd{kk,1}(1) == 3
                        probGNAT(kk).caseNum = 3;
                    end
                elseif (length(gnat.sampleInd{kk,1})==2)
                    if gnat.sampleInd{kk,1}(1) == 1 && gnat.sampleInd{kk,1}(2) == 2
                        probGNAT(kk).caseNum = 4;
                    elseif gnat.sampleInd{kk,1}(1) == 1 && gnat.sampleInd{kk,1}(2) == 3
                        probGNAT(kk).caseNum = 5;
                    elseif gnat.sampleInd{kk,1}(1) == 2 && gnat.sampleInd{kk,1}(2) == 3
                        probGNAT(kk).caseNum = 6;
                    end
                elseif (length(gnat.sampleInd{kk,1})==3)
                    probGNAT(kk).caseNum = 7;
                end
                if length(obj.B) > 1
                    probGNAT(kk).B = obj.B(gnat.sampleInd{kk,1},1);
                end
                

            end
        end
        
        function   [axOUT] = problemPaperPlot(obj,modelobj,modelAuxobj,axIN)
            %This function generates the Inverter Chain Circuit  plot
            %
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj         - Problem object
            %modelobj    - FOM, ROM, GNAT, or TPWL object
            %modelAuxobj - ROM object if modelobj is of class GNAT or
            %              locGNAT (unused otherwise)
            %pstr        - cell array indicating plot object formatting
            %axIN        - axis handle to place new plots
            %
            %Outputs:
            %--------
            %axOUT       - axis handle that new plots were placed
            %--------------------------------------------------------------
            
            %Extract reconstructed state vectors (with DBCs)

            svF = obj.returnPosSV(modelobj,modelAuxobj,[],'sv');
            N =size(svF,1);
            nT = size(svF,2);
            
            %Set up axes
            if isempty(axIN)
                ax = axes;
                title(' ','fontsize',16,'interpreter','latex');
                xlabel('node','fontsize',14,'interpreter','latex');
                ylabel('Time (s)','fontsize',14,'interpreter','latex');
                zlabel('Voltage','fontsize',14,'interpreter','latex');
            else
                ax = axIN;
            end
            hold on
           surf(1:N,1:nT,svF','EdgeColor','none');
            axOUT = ax;
        end
        
        function  [svF] = returnPosSV(obj,modelobj,modelAuxobj,~,flag,basenum)
            %This function takes a state vector from a simulation and
            %returns the state vector.  A vector of positions (including those with DBCs)
            %is also returned.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj         - Problem object
            %modelobj    - FOM, ROM, GNAT, or TPWL object
            %modelAuxobj - ROM object if modelobj is of class GNAT, empty
            %              otherwise
            %flag        - string ('pod','res','jac','clusCenter','sv')
            %basenum     - basis number to use (only needed for phi)
            %
            %Outputs:
            %--------
            %X           - x coordinate positions (empty for this problem)
            %Y           - y coordinate positions (empty for this problem)
            %svDBC       - state vector matrix with DBCs
            %--------------------------------------------------------------
            
            %Extract data (state vectors, cluster centers,
            %residual/jacobians, pod vectors)
            switch flag
                case 'pod'
                    svF = modelobj.phi(:,:,basenum);
                case 'res'
                    svF = modelobj.res;
                case 'jac'
                    svF = modelobj.jac;
                case 'clusCenter'
                    svF = modelobj.clusCenter;
                case 'sv'
                    svF = modelobj.reconstructFullState(modelAuxobj);
            end

        end
        
        function   [] = problemAnimate(obj,modelobj,modelAuxobj,pstr,~,freq)
            %This function generates the animation for the Nonlinear
            %Transmission line.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj      - Problem object
            %modelobj - FOM, ROM, GNAT, or TPWL object
            %modelAuxobj - ROM object if modelobj is of class GNAT or
            %              locGNAT (unused otherwise)
            %pstr        - cell array indicating plot object formatting
            %freq        - scalar indicating how often to plot sv in
            %              animation
            %
            %Outputs:
            %--------
            %There are no outputs for this function.  It generates an
            %animiation.
            %--------------------------------------------------------------
            
            if isempty(pstr)
                pstr = {'k.','MarkerSize',6};
            end
            
            %Extract reconstructed state vectors (with DBCs)
            [X,Y,svF] = obj.returnPosSV(modelobj,modelAuxobj,[],'sv');
            
            
            %Set up axes
            set(gcf,'position',[439   566   714   540]);
            ax = axes;
            title(' ','fontsize',16,'fontname','georgia');
            xlabel('x','fontsize',14,'interpreter','latex');
            ylabel('y','fontsize',14,'fontname','georgia');
            zlabel('z','fontsize',14,'fontname','georgia');
            hold on
            %set(ax,'nextplot','replacechildren','ylim',[ymin,ymax]);

            for i = 1:freq:size(svF,2)
                plot3(svF(1,i),svF(2,i),svF(3,i),pstr{:});
            end
        end
    end
    
    methods (Static)
       function [g,dg] = nonlinTerm(U,A,Vdd)
            %This function defines the nonlinear term of the nonlinear ODE.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %U   - state vector
            %
            %Outputs:
            %--------
            %g   - vector containing the evaluation of the nonlinear term
            %      at each component of the state vector
            %dg  - vector containing the evaluation of the derivative of
            %      the nonlinear term at each component of the state vector
            %--------------------------------------------------------------
            
            g = Vdd*tanh(-A*U);
            dg = -A*Vdd*(1-(tanh(-A*U)).^2);
        end
    end
end


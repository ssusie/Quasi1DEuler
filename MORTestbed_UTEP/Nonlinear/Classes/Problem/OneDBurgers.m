classdef OneDBurgers < handle
    
    properties (SetAccess = private, GetAccess = public)
        config;
        staticFlag = false;
        mesh;
        ndim;
        Jstruct;
        ic;
        B;
        G;
        C;
        D;
        dx;
        p;
    end
    
    properties (Hidden=true)
        Uind;
        ind1;
        lengthJROW;
    end
    
    methods
        function  [obj] = OneDBurgers(cfgobj,oldobj)
            %This is the constructor for the 1D Burger's Equation.  If
            %there is 1 input, the constructor will create the class
            %instance from the variables defined in this functino.  If
            %there are 2 inputs, the first input will be unused and the
            %output will be a deep copy of the oldobj.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %cfgobj - CONFIG object
            %oldobj - OneDBurgers object whose properties will be copied
            %         into the new object.
            %
            %Outputs:
            %--------
            %obj - OneDBurgers object
            %--------------------------------------------------------------
            
            if nargin == 0
                %Keep everything at defaults if no input arguments
                return;
            end
            
            % Construct a new object based on a deep copy of an old object
            % of this class by copying properties over.
            if nargin == 2 && isa(oldobj,'OneDBurgers')
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
            if isempty(obj.config.altFile)
                uniformMesh = 1;
            else
                uniformMesh = 0;
            end
            
            %Store number of DOFs
            ndof = cfgobj.nNodes(1) - 1;
            
            %Determine grid spacing
            if uniformMesh
                obj.mesh.node = linspace(cfgobj.DLim(1,1),cfgobj.DLim(1,2),cfgobj.nNodes(1))';
                %obj.dx = (cfgobj.DLim(1,2) - cfgobj.DLim(1,1))/(ndof);
            else
                obj.mesh.node = feval(obj.config.altFile);
                %obj.dx = obj.mesh.node(2:end) - obj.mesh.node(1:end-1);
            end
            obj.dx = obj.mesh.node(2:end) - obj.mesh.node(1:end-1);
            
            %Determine mesh
            obj.mesh.conn = [1:cfgobj.nNodes(1)-1;2:cfgobj.nNodes(1)]';
            obj.mesh.dbc = [1,true];
            
            %Number of unknowns per node
            obj.ndim = 1;
            
            %Determine grid spacing
            %obj.dx = obj.mesh.node(2,1) - obj.mesh.node(1,1);
            coord = obj.mesh.node(2:end,1); %want coordinate vector to only include nodes that remain after removing dbc
            
            %Setup Jacobian Structure
            JacobianStructure(obj);
            
            %Determine source term vector from the source function and the
            %coordinate vector
            if isa(cfgobj.GFunc,'function_handle')
                obj.G = cfgobj.GFunc(coord);
            else
                obj.G = feval(cfgobj.GFunc,coord);
            end
            
            %Determine ic term vector from the ic function and the
            %coordinate vector
            if isa(cfgobj.icFunc,'function_handle')
                obj.ic = cfgobj.icFunc(coord);
            else
                obj.ic = feval(cfgobj.icFunc,coord);
            end
            
            %Determine input matrix from the input matrix function
            switch obj.config.form
                case 'non_descriptor'
		    disp('non_descriptor')
                    if isa(cfgobj.BFunc,'function_handle')
                        obj.B = (0.5./obj.dx).*feval(cfgobj.BFunc);
                    else
                        obj.B = (0.5./obj.dx).*feval(cfgobj.BFunc);
                    end
                case 'descriptor'
                    disp('descriptor')
		    if isa(cfgobj.BFunc,'function_handle')
                        obj.B = 0.5*feval(cfgobj.BFunc);
                    else
                        obj.B = 0.5*feval(cfgobj.BFunc);
                    end
                case 'hybrid'
		    disp('hybrid')
                    if isa(cfgobj.BFunc,'function_handle')
                        obj.B = (0.5./sqrt(obj.dx)).*feval(cfgobj.BFunc);
                    else
                        obj.B = (0.5./sqrt(obj.dx)).*feval(cfgobj.BFunc);
                    end
            end
            
            %Determine output matrix from the output matrix function
            if isa(cfgobj.CFunc,'function_handle')
                obj.C = feval(cfgobj.CFunc);
            else
                obj.C = feval(cfgobj.CFunc);
            end
            
            %Determine feedthrough matrix from the feedthrough matrix function
            if isa(cfgobj.DFunc,'function_handle')
                obj.D = feval(cfgobj.DFunc);
            else
                obj.D = feval(cfgobj.DFunc);
            end
        end %Done
        
        function  [R,J] = ResJac(obj,U,t)
            %This function calcuates the residual and jacobian of Burger's Equation
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj - OneDBurgers object
            %U   - N x 1 vector containing the "solution" (or guess or
            %      approximation) at the current time
            %t   - scalar indicating the time
            %
            %Outputs:
            %--------
            %R   - N x 1 vector containing the nonlinear residual of the solution
            %      U at the current time t
            %J   - N x N vector containing the nonlinear jacobian of the solution
            %      U at the current time t
            %--------------------------------------------------------------
            
            %Determine the number of DOFs in the FOM
            N = obj.config.ndof;
            
            %Fill the entries of the residual and Jacobian appropriately
            switch obj.config.form
                case 'non_descriptor'
                    R = (1./obj.dx).*[-0.5*U(1,1)^2; 0.5*(U(1:N-1,1).^2-U(2:N,1).^2)] + ...
                        obj.G + obj.B*feval(obj.config.inFunc,t);
%                     keyboard 
                    dummy = obj.dx(end);
                    J = spdiags([(1./[obj.dx(2:end);dummy]).*U,-(1./obj.dx).*U],[-1,0],N,N);
                case 'descriptor'
                    R = [-0.5*U(1,1)^2; 0.5*(U(1:N-1,1).^2-U(2:N,1).^2)] + obj.dx.*obj.G + obj.B*feval(obj.config.inFunc,t);
                    J = spdiags([U,-U],[-1,0],N,N);
                case 'hybrid'
                    R = (1./sqrt(obj.dx)).*[-0.5*U(1,1)^2; 0.5*(U(1:N-1,1).^2-U(2:N,1).^2)] + ...
                        sqrt(obj.dx).*obj.G + obj.B*feval(obj.config.inFunc,t);
                    
                    dummy = sqrt(obj.dx(end));
                    J = spdiags([(1./[sqrt(obj.dx(2:end));dummy]).*U,-(1./sqrt(obj.dx)).*U],[-1,0],N,N);
            end
        end %Done
        
        function  [R,J] = ResJacGNAT(obj,U,t)
            %This function calcuates the residual and jacobian of Burger's
            %Equation 
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - OneDBurgers object
            %U       - N x 1 vector containing the "solution" (or guess or
            %          approximation) at the current time
            %t       - scalar indicating the time
            %indN    - vector containing the sample indicies
            %irstart - irstart vector from GNAT class
            %
            %Outputs:
            %--------
            %R   - indN x 1 vector containing the nonlinear residual of the
            %      solution U at the current time t
            %J   - vector containing the nonlinear jacobian of the solution
            %      U at the current time t at the rows indN
            %--------------------------------------------------------------
            
            switch obj.config.form
                case 'non_descriptor'
                    if obj.ind1 == 1
                        %This is done first so the entire array is initialized
                        R = [-U(1,1)*U(1,1)/(2*obj.dx(1));...
                            (U(obj.Uind,1).*U(obj.Uind,1)-U(obj.Uind+1,1).*U(obj.Uind+1,1))./(2*obj.dx(2:end))] + ...
                            obj.G + obj.B*feval(obj.config.inFunc,t);
                        J = [-U(1,1)./obj.dx(1); reshape([(U(obj.Uind,1)./(obj.dx(2:end)))'; -(U(obj.Uind+1,1)./(obj.dx(2:end)))'],obj.lengthJROW-1,1)];
                    else
                        R = (U(obj.Uind,1).*U(obj.Uind,1)-U(obj.Uind+1,1).*U(obj.Uind+1,1))./(2*obj.dx) + ...
                            obj.G + obj.B*feval(obj.config.inFunc,t);
                        J = reshape([(U(obj.Uind,1)./(obj.dx))'; -(U(obj.Uind+1,1)./(obj.dx))'],obj.lengthJROW,1);
                    end
                case 'descriptor'
                    if obj.ind1 == 1
                        %This is done first so the entire array is initialized
                        R = [-U(1,1)*U(1,1)/2;...
                            (U(obj.Uind,1).*U(obj.Uind,1)-U(obj.Uind+1,1).*U(obj.Uind+1,1))./(2)] + ...
                            obj.dx.*obj.G + obj.B*feval(obj.config.inFunc,t);
                        J = [-U(1,1); reshape([(U(obj.Uind,1))'; -(U(obj.Uind+1,1))'],obj.lengthJROW-1,1)];
                    else
                        R = (U(obj.Uind,1).*U(obj.Uind,1)-U(obj.Uind+1,1).*U(obj.Uind+1,1))./(2) + ...
                            obj.dx.*obj.G + obj.B*feval(obj.config.inFunc,t);
                        J = reshape([U(obj.Uind,1)'; -U(obj.Uind+1,1)'],obj.lengthJROW,1);
                    end
                    
                case 'hybrid'
                    if obj.ind1 == 1
                        %This is done first so the entire array is initialized
                        R = [-U(1,1)*U(1,1)/(2*sqrt(obj.dx(1)));...
                            (U(obj.Uind,1).*U(obj.Uind,1)-U(obj.Uind+1,1).*U(obj.Uind+1,1))./(2*sqrt(obj.dx(2:end)))] + ...
                            sqrt(obj.dx).*obj.G + obj.B*feval(obj.config.inFunc,t);
                        J = [-U(1,1)./sqrt(obj.dx(1)); reshape([(U(obj.Uind,1)./(sqrt(obj.dx(2:end))))'; -(U(obj.Uind+1,1)./(sqrt(obj.dx(2:end))))'],obj.lengthJROW-1,1)];
                    else
                        R = (U(obj.Uind,1).*U(obj.Uind,1)-U(obj.Uind+1,1).*U(obj.Uind+1,1))./(2*sqrt(obj.dx)) + ...
                            sqrt(obj.dx).*obj.G + obj.B*feval(obj.config.inFunc,t);
                        J = reshape([(U(obj.Uind,1)./sqrt(obj.dx))'; -(U(obj.Uind+1,1)./sqrt(obj.dx))'],obj.lengthJROW,1);
                    end
            end           
            
        end %Done
        
        function  [R,J] = ResJacTimeDer(obj,U,~)
           switch obj.config.form
                case 'non_descriptor'
                    R = U;
                    J = speye(size(U,1));
                case 'descriptor'
                    N = obj.config.ndof;
                    R = obj.dx.*U;
                    J = spdiags(obj.dx,0,N,N);
                case 'hybrid'
                    N = obj.config.ndof;
                    R = sqrt(obj.dx).*U;
                    J = spdiags(sqrt(obj.dx),0,N,N);                
            end
        end %Done
        
        function  [R,J] = ResJacGNATTimeDer(obj,U,~,~,~,~,jdiag)
            switch obj.config.form
                case 'non_descriptor'
                    R = U;
                    J = jdiag;
                    
                case 'descriptor'
                    R = obj.dx.*U;
                    if obj.ind1 == 1
                        J = [obj.dx(1);  reshape([zeros(size(obj.dx(2:end)'));obj.dx(2:end)'],2*size(obj.dx(2:end),1),1)];
                    else
                        J = reshape([zeros(size(obj.dx'));obj.dx'],2*size(obj.dx,1),1);
                    end
                case 'hybrid'
                    R = sqrt(obj.dx).*U;
                    if obj.ind1 == 1
                        J = [sqrt(obj.dx(1));  reshape([zeros(size(obj.dx(2:end)'));sqrt(obj.dx(2:end))'],2*size(obj.dx(2:end),1),1)];
                    else
                        J = reshape([zeros(size(obj.dx'));sqrt(obj.dx)'],2*size(obj.dx,1),1);
                    end
            end
        end %Done
        
        function  [R,J] = ResJacWithoutInput(obj,U,t)
            %This function returns the residual and jacobian without the
            %input term contributions.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj - OneDBurgers object
            %
            %Outputs:
            %--------
            %Radj - residual without input term contribution
            %Jadj - jacobian without input term contribution
            %--------------------------------------------------------------
            
            [R,J] = obj.ResJac(U,t);
            R = R - obj.B*feval(obj.config.inFunc,t);
            
        end %Done
        
        function  [ind] = node2ind(obj,nodenum) 
            %This function takes a node in the mesh and maps it to an index
            %in the full state vector
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - OneDBurgers object
            %nodenum - scalar specifying node numbers
            %
            %Outputs:
            %--------
            %ind     - vector specifying the corresponding indices in the
            %          full state vector
            %--------------------------------------------------------------
            
            if nodenum == 1
                ind = nan;
            else
                ind = nodenum - 1;
            end
        end %Done
    
        function  [] = JacobianStructure(obj)
            %This function forms a boolean matrix whose true entries
            %correspond to nonzero entries in the jacobian
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - OneDBurgers object
            %
            %Outputs:
            %--------
            %There are no outputs.  It store a ndof x ndof boolean matrix
            %whose true values correspond to nonzero entries in the
            %jacobian in the class instance.
            %--------------------------------------------------------------
            
            N = obj.config.ndof;
            obj.Jstruct = spdiags([ones(N,1),ones(N,1)],[-1,0],N,N);            
        end %Done
                
        function  [probGNAT] = createCopy4GNAT(obj,gnat)
            %This function makes a copy of the current object and only
            %copies the source term and input matrix evaluated at indN
            %rows.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - OneDBurgers object
            %gnat    - gnat object
            %
            %Outputs:
            %--------
            %probGNAT - new "reduced" OneDBurgers object
            %--------------------------------------------------------------
            
            probGNAT = OneDBurgers([],obj);
            
            temp = probGNAT.config.inFunc;
            temp2 = probGNAT.config.form;
            probGNAT.config = [];
            probGNAT.config.inFunc = temp;
            probGNAT.config.form = temp2;
            probGNAT.mesh   = [];
            probGNAT.ndim = [];
            probGNAT.Jstruct = [];
        
            probGNAT.ind1 = gnat.sampleInd(1);
            if probGNAT.ind1 == 1
                probGNAT.Uind = gnat.irstart(2:end-1);
            else
                probGNAT.Uind = gnat.irstart(1:end-1);
            end
            probGNAT.lengthJROW = gnat.irstart(end)-1;
            
            if length(obj.B) > 1
                probGNAT.B = obj.B(gnat.sampleInd,1);
            end
            probGNAT.dx = obj.dx(gnat.sampleInd);
            
            if length(obj.G) > 1
                probGNAT.G = obj.G(gnat.sampleInd,1);
            end
            
            if length(obj.ic) > 1
                probGNAT.ic = obj.ic(unique(gnat.jrow),1);
            end
        end %Done
        
        function  [probGNAT] = createCopy4GNATloc(obj,gnat)
            %This function makes a copy of the current object and only
            %copies the source term and input matrix evaluated at indN
            %rows.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - OneDBurgers object
            %gnat    - gnat object
            %
            %Outputs:
            %--------
            %probGNAT - array of new "reduced" OneDBurgers object
            %--------------------------------------------------------------
           
            for kk = 1:gnat.nBases
                if isempty(gnat.sampleInd{kk,1})
                    continue;
                end
                probGNAT(kk) = OneDBurgers([],obj);
                
                temp = probGNAT(kk).config.inFunc;
                temp2 = probGNAT(kk).config.form;
                
                probGNAT(kk).config = [];
                probGNAT(kk).config.inFunc = temp;
                probGNAT(kk).config.form = temp2;
                probGNAT(kk).mesh   = [];
                probGNAT(kk).ndim = [];
                probGNAT(kk).Jstruct = [];
                
                probGNAT(kk).dx = obj.dx(gnat.sampleInd{kk,1});

                probGNAT(kk).ind1 = (gnat.sampleInd{kk,1}(1) == 1);
                if probGNAT(kk).ind1 == 1
                    probGNAT(kk).Uind = gnat.irstart{kk,1}(2:end-1);
                else
                    probGNAT(kk).Uind = gnat.irstart{kk,1}(1:end-1);
                end
                probGNAT(kk).lengthJROW = gnat.irstart{kk,1}(end)-1;
                
                if length(obj.B) > 1
                    probGNAT(kk).B = obj.B(gnat.sampleInd{kk,1},1);
                end
                
                if length(obj.G) > 1
                    probGNAT(kk).G = obj.G(gnat.sampleInd{kk,1},1);
                end
                
                if length(obj.ic) > 1
                    probGNAT(kk).ic = obj.ic(unique(gnat.jrow{kk,1}),1);
                end
            end
        end %Done
        
        function   [axOUT,lineobj] = problemPaperPlot(obj,modelobj,modelAuxobj,Tplots,pstr,axIN)
            %This function generates the 1D Burger's Equation plot from the
            %Rewienski thesis.
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
            %axOUT    - axis handle that new plots were placed
            %--------------------------------------------------------------
            
            %Store the times instances to plot the solution (based on
            %Reweinski thesis)
%             Tplots = 0.35;
%             Tplots = [0.30,2.5,7,10,20,30,40,50];
            
            %Extract reconstructed state vectors (with DBCs)
            [~,~,svF] = obj.returnPosSV(modelobj,modelAuxobj,[],'sv');
            
            %Extract mesh spacing
            x = obj.mesh.node;
            n = length(Tplots);
            tvec = modelobj.time.T(1):modelobj.time.dt:modelobj.time.T(2);
            ind = zeros(1,2);
            
            %Set up axes
            if isempty(axIN)
                ax = axes;
                title('Shock Propagation Through Fluid','fontsize',16,'interpreter','latex');
                xlabel('Position','fontsize',14,'interpreter','latex');
                ylabel('Conserved Quantity','fontsize',14,'interpreter','latex');
            else
                ax = axIN;
            end
            set(ax,'nextplot','add');
            
            for i = 1:n
                %Determine indices in svF that correspond to entries in
                %Tplots
                if sum(abs(Tplots(i) - tvec) < 1e-8) == 0
                    [~,temp] = sort(abs(tvec - Tplots(i)));
                    ind = temp(1:2);
                else
                    ind(1,:) = find(abs(Tplots(i) - tvec) < 1e-8);
                end
                U = 0.5*(svF(:,ind(1,1)) + svF(:,ind(1,2)));
                
                lineobj(i) = plot(ax,x,U,pstr{:});
            end
            axOUT = ax;
        end %Done
        
        function  [X,Y,svDBC] = returnPosSV(obj,modelobj,modelAuxobj,~,flag,basenum)
            %This function takes a state vector from a simulation and
            %returns the state vector with the Dirichlet boundary
            %conditions.  A vector of positions (including those with DBCs)
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
            %X           - x coordinate positions
            %Y           - y coordinate positions (empty for this problem)
            %svDBC       - state vector matrix with DBCs
            %--------------------------------------------------------------
            
            %Extract data (state vectors, cluster centers,
            %residual/jacobians, pod vectors)
            switch flag
                case 'pod'
                    svF = modelobj.phi{basenum,1};
                    %svF = modelobj.phi(:,:,basenum);
                    
                    u = zeros(1,size(svF,2));
                case 'res'
                    svF = modelobj.res;
                    
                    u = zeros(1,size(svF,2));
                case 'jac'
                    svF = modelobj.jac;
                    
                    u = zeros(1,size(svF,2));
                case 'clusCenter'
                    svF = modelobj.clusCenter;
                    
                    u = zeros(1,size(svF,2));
                case 'sv'
                    svF = modelobj.reconstructFullState(modelAuxobj);
                    
                    t = modelobj.time.T(1):modelobj.time.dt:modelobj.time.T(2);
                    u = feval(obj.config.inFunc,t);
            end
            n = size(svF,2);

            svDBC = zeros(obj.config.ndof+1,n);
            if ~isempty(svF)
                svDBC(1,:) = sqrt(u);
                svDBC(2:end,:) = svF;
            else
                error('Requested empty entry (most likely clusCenters or pod).');
            end
            
            X = obj.mesh.node;
            Y = [];
            
        end %Done
        
        function   [] = problemAnimate(obj,modelobj,modelAuxobj,pstr,~,freq)
            %This function generates a problem specific 1D Burger's
            %Equation plot.
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
            
            %Extract reconstructed state vectors (with DBCs) and mesh
            %spacing
            [x,~,svF] = obj.returnPosSV(modelobj,modelAuxobj);

            tvec = modelobj.time.T(1):modelobj.time.dt:modelobj.time.T(2);
            
            %Determine ylimits
            ymin = min(min(svF)); ymax = max(max(svF));
            
            set(gcf,'position',[439   566   714   540]);
            ax = axes;
            title(ax,'Shock Propagation Through Fluid','fontsize',16,'fontname','georgia');
            xlabel(ax,'Position','fontsize',14,'fontname','georgia');
            ylabel(ax,'Conserved Quantity','fontsize',14,'fontname','georgia');
            set(ax,'nextplot','replacechildren','ylim',[ymin,ymax]);
            for i = 1:freq:length(tvec)
                plot(ax,x,svF(:,i),pstr{:}); pause(0.001);
            end
        end %Done

        function sqrtMass = returnSqrtMass(obj)
            switch obj.config.form
                case 'non_descriptor'
                    sqrtMass = speye(length(obj.dx));
                case 'descriptor'
                    sqrtMass = spdiags(obj.dx,0,length(obj.dx),length(obj.dx));
                case 'hybrid'
                    sqrtMass = spdiags(sqrt(obj.dx),0,length(obj.dx),length(obj.dx));
            end
        end
        
        
        
        
        function [] = setPrevSV(obj,iprevsv)

        end
    end
end

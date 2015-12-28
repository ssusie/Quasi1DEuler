classdef OneDFKPP < handle
    
    properties (SetAccess = private, GetAccess = public)
        config;
        staticFlag = false;
        mesh;
        ndim;
        Jstruct;
        ic;
        dx;
        Delta;
        Mass;
        sqrtMass;
        alpha;
        discretization;
    end
    
    properties (Hidden=true)
        Uind;
        ind1;
        lengthJROW;
    end
    
    methods
        function  [obj] = OneDFKPP(cfgobj,oldobj)
            %This is the constructor for the 1D FKPP Equation.  If
            %there is 1 input, the constructor will create the class
            %instance from the variables defined in this function.  If
            %there are 2 inputs, the first input will be unused and the
            %output will be a deep copy of the oldobj.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %cfgobj - CONFIG object
            %oldobj - OneDFKPP object whose properties will be copied
            %         into the new object.
            %
            %Outputs:
            %--------
            %obj - OneDFKPP object
            %--------------------------------------------------------------
            
            if nargin == 0
                %Keep everything at defaults if no input arguments
                return;
            end
            
            % Construct a new object based on a deep copy of an old object
            % of this class by copying properties over.
            if nargin == 2 && isa(oldobj,'OneDFKPP')
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

            obj.discretization = 'FEM';
            % obj.discretization = 'FD';
            
            %Store number of DOFs
            ndof = cfgobj.nNodes(1) - 2; % Dirichlet BC with FD
            
            %Determine grid spacing
            [obj.mesh.node,obj.alpha] = feval(obj.config.altFile);
            obj.dx = obj.mesh.node(2:end) - obj.mesh.node(1:end-1);


            %Determine mesh
            obj.mesh.conn = [1:cfgobj.nNodes(1)-1;2:cfgobj.nNodes(1)]';
            obj.mesh.dbc = [1,true;...
                cfgobj.nNodes(1), true];
            
            %Number of unknowns per node
            obj.ndim = 1;
            
            %Determine grid spacing
            coord = obj.mesh.node(2:end-1,1); %want coordinate vector to only include nodes that remain after removing dbc
            
            %Setup Jacobian Structure
            JacobianStructure(obj);
            
            switch obj.discretization
                case 'FD'
                    BuildLaplacian_FD(obj);
                case 'FEM'
                    BuildStiffness_FEM(obj);
                    %BuildMass_FEM(obj);   
                    BuildMass_FEM_Lumped(obj); 
                    if (strcmp(obj.config.form,'hybrid'))
                        obj.sqrtMass = sparse(sqrtm(full(obj.Mass)));
                    end
            end
            %Determine ic term vector from the ic function and the
            %coordinate vector
            if isa(cfgobj.icFunc,'function_handle')
                obj.ic = cfgobj.icFunc(coord);
            else
                obj.ic = feval(cfgobj.icFunc,coord);
            end
            

        end
        
        function  [R,J] = ResJac(obj,U,t)
            %This function calcuates the residual and jacobian of the FKPP's Equation
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj - OneDFKPP object
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
            %Fill the entries of the residual appropriately

            switch obj.discretization
                case 'FD'
                    R = obj.Delta*U + obj.alpha*U.*(ones(size(U))-U);
                    J = obj.Delta + spdiags(obj.alpha*(ones(size(U))-2*U),0,N,N);
                    
                case 'FEM'
                    [F_NL,J_NL] = nonLinearTerm_FEM(obj,U);
                    switch obj.config.form
                        
                        case 'non_descriptor'
                            
                            R = obj.Mass\(-obj.Delta*U - obj.alpha*F_NL) + obj.alpha*U;
                            J = obj.Mass\(-obj.Delta -obj.alpha*J_NL) + obj.alpha*speye(size(U,1));
                            
                        case 'descriptor'
                            R = -obj.Delta*U - obj.alpha*F_NL + obj.alpha* obj.Mass*U;
                            J = -obj.Delta -obj.alpha*J_NL + obj.alpha* obj.Mass;
                        case 'hybrid'
                            R = obj.sqrtMass\(-obj.Delta*U - obj.alpha*F_NL) + obj.alpha*obj.sqrtMass*U;
                            J = obj.sqrtMass\(-obj.Delta -obj.alpha*J_NL) + obj.alpha*obj.sqrtMass;
                    end
            end
        end
        
        function  [R,J] = ResJacGNAT(obj,U,t)
            %This function calculates the residual and jacobian of the FKPP's
            %Equation 
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - OneDFKPP object
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
                            R = 0
                            J = 0
                        else
                            R = 0
                            J = 0
                        end       
                case 'descriptor'
                        if obj.ind1 == 1
                            %This is done first so the entire array is initialized
                            R = 0
                            J = 0
                        else
                            R = 0
                            J = 0
                        end                     
                    
                case 'hybrid'
                        if obj.ind1 == 1
                            %This is done first so the entire array is initialized
                            R = 0
                            J = 0
                        else
                            R = 0
                            J = 0
                        end                      
                    
            end


        end
        
        function  [R,J] = ResJacTimeDer(obj,U,~)
            
            switch obj.discretization
                case 'FD'
                    R = U;
                    J = speye(size(U,1));
                case 'FEM'
                    

                    switch obj.config.form
                        case 'non_descriptor'
                            R = U; 
                            J = speye(size(U,1)); 
                        case 'descriptor'
                            
                            R = obj.Mass*U;
                            J = obj.Mass;
                        case 'hybrid'
                            R = obj.sqrtMass*U;
                            J = obj.sqrtMass;                

                    end
            end
        end
        
        function  [R,J] = ResJacGNATTimeDer(obj,U,~,~,~,~,jdiag)
           switch obj.config.form
                case 'non_descriptor'
                    R = U;
                    J = jdiag;
        
               case 'descriptor'
                    R = 0;% obj.dx.*U;
%                     if obj.ind1 == 1
%                         J = [obj.dx(1);  reshape([zeros(size(obj.dx(2:end)'));obj.dx(2:end)'],2*size(obj.dx(2:end),1),1)];
%                     else
%                         J = reshape([zeros(size(obj.dx'));obj.dx'],2*size(obj.dx,1),1);
%                     end             
               case 'hybrid'
                    R = 0;%sqrt(obj.dx).*U;
%                     if obj.ind1 == 1
%                         J = [sqrt(obj.dx(1));  reshape([zeros(size(obj.dx(2:end)'));sqrt(obj.dx(2:end))'],2*size(obj.dx(2:end),1),1)];
%                     else
%                         J = reshape([zeros(size(obj.dx'));sqrt(obj.dx)'],2*size(obj.dx,1),1);
%                     end                    
           end
        end
        
        function  [R,J] = ResJacWithoutInput(obj,U,t)
            %This function returns the residual and jacobian without the
            %input term contributions.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj - OneDFKPP object
            %
            %Outputs:
            %--------
            %Radj - residual without input term contribution
            %Jadj - jacobian without input term contribution
            %--------------------------------------------------------------
            
            [R,J] = obj.ResJac(U,t);
            R = R;
            
        end
        
        function  [ind] = node2ind(obj,nodenum) 
            %This function takes a node in the mesh and maps it to an index
            %in the full state vector
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - OneDFKPP object
            %nodenum - scalar specifying node numbers
            %
            %Outputs:
            %--------
            %ind     - vector specifying the corresponding indices in the
            %          full state vector
            %--------------------------------------------------------------
            
            if nodenum == 1 || nodenum == obj.config.nNodes(1)
                ind = nan;
            else
                ind = nodenum - 1;
            end
        end
    
        function  [] = JacobianStructure(obj)
            %This function forms a boolean matrix whose true entries
            %correspond to nonzero entries in the jacobian
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - OneDFKPPs object
            %
            %Outputs:
            %--------
            %There are no outputs.  It store a ndof x ndof boolean matrix
            %whose true values correspond to nonzero entries in the
            %jacobian in the class instance.
            %--------------------------------------------------------------
            
            N = obj.config.ndof;
            obj.Jstruct = spdiags([ones(N,1),ones(N,1),ones(N,1)],[-1,0,1],N,N);
%             obj.Jstruct = sparse(temp);
            
        end 
           
        function [] = BuildLaplacian_FD(obj)
            N = obj.config.ndof;

            dxm = obj.dx(1:end-1,1);   
            dxp = obj.dx(2:end,1);
            size(dxm)
            size(dxp)
            lowD = 2./(dxm(2:end,1).*(dxm(2:end,1)+dxp(2:end,1)));
            D = -2*(1./dxm+1./dxp)./(dxm+dxp);
            upD =  2./(dxp(1:end-1,1).*(dxm(1:end-1,1)+dxp(1:end-1,1)));

           obj.Delta = spdiags([[lowD;0],D,[0;upD]],[-1 0 1],N,N);

           % obj.Delta = 1/(obj.dx(1)^2)*spdiags([ones(N,1),-2*ones(N,1),ones(N,1)],[-1 0 1],N,N);

        end
        
        function [] = BuildStiffness_FEM(obj)
            N = obj.config.ndof;

            %Elementary matrix
            K_e = [1 -1;-1 1];
            
            obj.Delta = sparse(N,N);
            obj.Delta(1,1) = obj.Delta(1,1) + K_e(2,2)/obj.dx(1,1);
            for i=2:size(obj.dx,1)-1
                obj.Delta(i-1:i,i-1:i) = obj.Delta(i-1:i,i-1:i) + K_e/obj.dx(i,1);  
            end
            obj.Delta(N,N) = obj.Delta(N,N) + K_e(1,1)/obj.dx(end,1);
            
        end
        
        function [] = BuildMass_FEM(obj)
             N = obj.config.ndof;

            %Elementary matrix
            M_e = [1/3 1/6;1/6 1/3];       
            obj.Mass = sparse(N,N);
            obj.Mass(1,1) = obj.Mass(1,1) + M_e(2,2)*obj.dx(1,1);
            for i=2:size(obj.dx,1)-1
                obj.Mass(i-1:i,i-1:i) = obj.Mass(i-1:i,i-1:i) + M_e*obj.dx(i,1);  
            end
            obj.Mass(N,N) = obj.Mass(N,N) + M_e(1,1)*obj.dx(end,1);
             
        end
        
        function [] = BuildMass_FEM_Lumped(obj)
            N = obj.config.ndof;
            
            %Elementary matrix
            M_el = [1/2 0;0 1/2];
            M_e = [1/3 1/6;1/6 1/3]; 
            obj.Mass = sparse(N,N);
            obj.Mass(1,1) = obj.Mass(1,1) + M_e(2,2)*obj.dx(1,1);
            for i=2:size(obj.dx,1)-1
                obj.Mass(i-1:i,i-1:i) = obj.Mass(i-1:i,i-1:i) + M_el*obj.dx(i,1);
            end
            obj.Mass(N,N) = obj.Mass(N,N) + M_e(1,1)*obj.dx(end,1);
            
        end
        
        function [F,J] = nonLinearTerm_FEM(obj,U)
            N = obj.config.ndof;

            %Elementary matrices
            M_e1 = [1/4 1/12;1/12 1/12];
            M_e2 = [1/12 1/12;1/12 1/4];
            
            F = zeros(N,1);
            F(1,1) = F(1,1) + obj.dx(1,1)*M_e2(2,2)*U(1,1)^2;
             for i=2:size(obj.dx,1)-1
                F(i-1:i,1) = F(i-1:i,1) + [U(i-1:i,1)'*M_e1*U(i-1:i,1);U(i-1:i,1)'*M_e2*U(i-1:i,1)]*obj.dx(i,1) ;  
            end           
            F(N,1) = F(N,1) + obj.dx(end,1)*M_e1(1,1)*U(N,1)^2;
            
            J = sparse(N,N);
            J(1,1) = J(1,1) + 2*obj.dx(1,1)*M_e2(2,2)*U(1,1);
            for i=2:size(obj.dx,1)-1
                J(i-1:i,i-1:i) = J(i-1:i,i-1:i) + [U(i-1:i,1)'*M_e1;U(i-1:i,1)'*M_e2]*obj.dx(i,1) ; 
            end
            J(N,N) = J(N,N) + 2*obj.dx(end,1)*M_e1(1,1)*U(N,1);       
            
            
        end
        
        
        function  [probGNAT] = createCopy4GNAT(obj,gnat)
            %This function makes a copy of the current object and only
            %copies the source term and input matrix evaluated at indN
            %rows.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - OneDFKPP object
            %gnat    - gnat object
            %
            %Outputs:
            %--------
            %probGNAT - new "reduced" OneDFKPP object
            %--------------------------------------------------------------
            
            probGNAT = OneDFKPP([],obj);
            
            temp = probGNAT.config.inFunc;
            temp2 = probGNAT.config.form;
            probGNAT.config = [];
            probGNAT.config.inFunc = temp;
            probGNAT.config.form = temp2;
            probGNAT.mesh   = [];
            probGNAT.ndim = [];
            probGNAT.Jstruct = [];
            
           disp('TO DO start');
           pause
            probGNAT.ind1 = gnat.sampleInd(1);
            if probGNAT.ind1 == 1
                probGNAT.Uind = gnat.irstart(2:end-1);
            else
                probGNAT.Uind = gnat.irstart(1:end-1);
            end
            probGNAT.lengthJROW = gnat.irstart(end)-1;
            disp('TO DO end');

            probGNAT.dx = obj.dx(gnat.sampleInd);
                    
            if length(obj.ic) > 1
                probGNAT.ic = obj.ic(unique(gnat.jrow),1);
            end
        end
        
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
%                 probGNAT(kk).ic = obj.ic(jrowUnion);
                probGNAT(kk).config = [];
                probGNAT(kk).config.inFunc = temp;
                probGNAT(kk).mesh   = [];
                probGNAT(kk).ndim = [];
                probGNAT(kk).Jstruct = [];
   
                            
               disp('TO DO start');
               pause
                probGNAT(kk).ind1 = gnat.sampleInd{kk,1}(1);
                if probGNAT(kk).ind1 == 1
                    probGNAT(kk).Uind = gnat.irstart{kk,1}(2:end-1);
                else
                    probGNAT(kk).Uind = gnat.irstart{kk,1}(1:end-1);
                end
                probGNAT(kk).lengthJROW = gnat.irstart{kk,1}(end)-1;
                disp('TO DO end');
                
                if length(obj.ic) > 1
                    probGNAT(kk).ic = obj.ic(unique(gnat.jrow{kk,1}),1);
                end
            end
        end
        
        function   [axOUT,lineobj] = problemPaperPlot(obj,modelobj,modelAuxobj,Tplots,pstr,axIN)
            %This function generates the 1D FKPP Equation plot from Jean Frederic Gerbeau's paper.
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
                xlabel('Position','fontsize',14,'interpreter','latex');
                ylabel('U','fontsize',14,'interpreter','latex');
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
                U = svF(:,ind(1,1));
                
                lineobj(i) = plot(ax,x,U,pstr{:});
            end
            axOUT = ax;
        end
        
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

            svDBC = zeros(obj.config.ndof+2,n);
            if ~isempty(svF)
                svDBC(1,:) = zeros(1,n);
                svDBC(2:end-1,:) = svF;
                svDBC(end,:) = zeros(1,n);
            else
                error('Requested empty entry (most likely clusCenters or pod).');
            end
            
            X = obj.mesh.node;
            Y = [];
            
        end
        
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
        end
    end
end
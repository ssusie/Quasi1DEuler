classdef OneDKdV < handle
    
    properties (SetAccess = private, GetAccess = public)
        config;
        staticFlag = false;
        mesh;
        ndim;
        Jstruct;
        ic;
        dx;
        LinMat;
        eta;
        dm2;
        dm1;
        d0;
        dp1;
        dp2;
    end
    
    properties (Hidden=true)
        Uind;
        ind1;
        ind2;
        indendm1;
        indend;
        lengthJROW;
        ndof;
        D0;
    end
    
    methods
        function  [obj] = OneDKdV(cfgobj,oldobj)
            %This is the constructor for the 1D KdV Equation.  If
            %there is 1 input, the constructor will create the class
            %instance from the variables defined in this function.  If
            %there are 2 inputs, the first input will be unused and the
            %output will be a deep copy of the oldobj.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %cfgobj - CONFIG object
            %oldobj - OneDKdV object whose properties will be copied
            %         into the new object.
            %
            %Outputs:
            %--------
            %obj - OneDKdV object
            %--------------------------------------------------------------
            
            if nargin == 0
                %Keep everything at defaults if no input arguments
                return;
            end
            
            % Construct a new object based on a deep copy of an old object
            % of this class by copying properties over.
            if nargin == 2 && isa(oldobj,'OneDKdV')
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

            
            %Store number of DOFs
            obj.ndof = cfgobj.nNodes(1) - 1; % Periodic BC
            
            %Determine grid spacing
            [obj.mesh.node,obj.eta] = feval(obj.config.altFile); 
            dxtemp = obj.mesh.node(2:end) - obj.mesh.node(1:end-1);
            if (norm((dxtemp-dxtemp(1)*ones(size(dxtemp)))/dxtemp(1)) > 10^(-10))
                disp('Non-uniform mesh not supported in KdV model');
                return;
            end
            obj.dx=dxtemp(1);
            
            %Determine mesh
            obj.mesh.conn = [1:cfgobj.nNodes(1)-1;2:cfgobj.nNodes(1)]';
            obj.mesh.pbc = [1,true]; 
            
            %Number of unknowns per node
            obj.ndim = 1;
            
            %Determine grid spacing
            coord = obj.mesh.node(1:end-1,1); %want coordinate vector to only include nodes that remain after removing pbc
            
            %Setup Jacobian Structure
            JacobianStructure(obj); 
            
            %Determine ic term vector from the ic function and the
            %coordinate vector 
            if isa(cfgobj.icFunc,'function_handle')
                obj.ic = cfgobj.icFunc(coord);
            else
                obj.ic = feval(cfgobj.icFunc,coord);
            end

            
            % set up linear operators
            N = obj.config.ndof;
            e = ones(N,1);
           % obj.dxm = [obj.dx(end);obj.dx(1:end-1)];
            Dm = spdiags([-1/obj.dx*e 1/obj.dx*e], -1:0, N, N);
            Dm(1,N) = -1/obj.dx;
            Dp = spdiags([-(1/obj.dx)*e (1/obj.dx)*e], -0:1, N, N);
            Dp(N,1) = 1/obj.dx;
            obj.D0 = spdiags([-1/(2*obj.dx)*e 0*e 1/(2*obj.dx)*e], -1:1, N, N);
            obj.D0(1,N) = -1/(2*obj.dx(1));
            obj.D0(N,1) = 1/(2*obj.dx(end));
            % artificial viscosity eta
            %DeltaMat = bsxfun(@times,obj.dx,obj.eta*(Dp*Dm*Dp*Dm));
            DeltaMat = obj.dx*obj.eta*(Dp*Dm*Dp*Dm);
            obj.LinMat = Dp*obj.D0*Dm + DeltaMat;

            obj.dm2 = obj.LinMat(3,1);
            obj.dm1 = obj.LinMat(3,2);
            obj.d0 = obj.LinMat(3,3);
            obj.dp1 = obj.LinMat(3,4);
            obj.dp2 = obj.LinMat(3,5);

        end
        
        function  [R,J] = ResJac(obj,U,t)
            %This function calcuates the residual and jacobian of the KdV's Equation
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj - OneDKdV object
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


            %Fill the entries of the residual appropriately
            N = length(U);

            R = -obj.LinMat*U -2* (U.*(obj.D0*U) + obj.D0*(U.*U));
          %  R2 = -(obj.dm2*[U(end-1:end,1);U(1:end-2,1)] + obj.dm1*[U(end,1);U(1:end-1,1)] + obj.d0*U + obj.dp1*[U(2:end,1);U(1,1)] + obj.dp2*[U(3:end,1);U(1:2,1)])  -2* (1./(obj.dx+obj.dxm).*(U+[U(2:end,1);U(1,1)]+[U(end,1);U(1:end-1,1)]).*([U(2:end,1);U(1,1)]-[U(end,1);U(1:end-1,1)]));

            J = -obj.LinMat - 2*( bsxfun(@times,U,obj.D0)+2*bsxfun(@times,obj.D0',U)');
%                   J2 =   sparse(N,N);
%                   % first row
%                   J2(1,end-1) = -obj.dm2;
%                   J2(1,end) = -obj.dm1 + 2*(2*U(end,1)+U(1,1))/(2*obj.dx);
%                   J2(1,1) = -obj.d0 - 2*(U(2,1)-U(end,1))/(2*obj.dx);
%                   J2(1,2) = -obj.dp1 - 2*(U(1,1)+2*U(2,1))/(2*obj.dx);
%                   J2(1,3) = -obj.dp2;
%                   
%                   % second row
%                   J2(2,end) = -obj.dm2;
%                   J2(2,1) = -obj.dm1 + 2*(2*U(1,1)+U(2,1))/(obj.dx(2)+obj.dxm(2));
%                   J2(2,2) = -obj.d0 - 2*(U(3,1)-U(1,1))/(obj.dx(2)+obj.dxm(2));
%                   J2(2,3) = -obj.dp1 - 2*(U(2,1)+2*U(3,1))/(obj.dx(2)+obj.dxm(2));
%                   J2(2,4) = -obj.dp2;
%                   % middle
%                   e = ones(size(J,1)-4,1);
%                   size(e)
%                   size(2*(2*U(2:end-3,1)+U(3:end-2,1))./(obj.dx(3:end-2,1)+obj.dxm(3:end-2,1)))
%                   size(2*(U(4:end-1,1)-U(2:end-3,1))./(obj.dx(3:end-2,1)+obj.dxm(3:end-2,1)))
%                   size(2*(U(2:end-3,1)+2*U(3:end-2,1))./(obj.dx(3:end-2,1)+obj.dxm(3:end-2,1)))
%                   size(J2)
%                   size(J2(3:N-2,3:N-2))
%                   J2(3:N-2,3:N-2) = spdiags([-obj.dm2*e,-obj.dm1*e+2*(2*U(2:end-3,1)+U(3:end-2,1))./(obj.dx(3:end-2,1)+obj.dxm(3:end-2,1)),-obj.d0*e-2*(U(4:end-1,1)-U(2:end-3,1))./(obj.dx(3:end-2,1)+obj.dxm(3:end-2,1)), -obj.dp1*e-2*(U(2:end-3,1)+2*U(3:end-2,1))./(obj.dx(3:end-2,1)+obj.dxm(3:end-2,1)),-obj.dp2*e],-2:2,size(J,1)-4,size(J,1)-4);
%                   % second to last row
%                   J2(end-1,end-3) = -obj.dm2;
%                   J2(end-1,end-2) = -obj.dm1 + 2*(2*U(end-2,1)+U(end-1,1))/(obj.dx(end-1)+obj.dxm(end-1));
%                   J2(end-1,end-1) = -obj.d0 - 2*(U(end,1)-U(end-2,1))/(obj.dx(end-1)+obj.dxm(end-1));
%                   J2(end-1,end) = -obj.dp1 - 2*(U(end-1,1)+2*U(end,1))/(obj.dx(end-1)+obj.dxm(end-1));
%                   J2(end-1,1) = -obj.dp2;         
%                   % last row
%                   J2(end,end-2) = -obj.dm2;
%                   J2(end,end-1) = -obj.dm1 + 2*(2*U(end-1,1)+U(end,1))/(obj.dx(end)+obj.dxm(end));
%                   J2(end,end) = -obj.d0 - 2*(U(1,1)-U(end-1,1))/(obj.dx(end)+obj.dxm(end));
%                   J2(end,1) = -obj.dp1 - 2*(U(end,1)+2*U(1,1))/(obj.dx(end)+obj.dxm(end));
%                   J2(end,2) = -obj.dp2;                   
%    normest(J-J2,2)/normest(J,2) 
%    pause
        end
        
        function  [R,J] = ResJacGNAT(obj,U,t)
            %This function calculates the residual and jacobian of the KdV's
            %Equation 
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - OneDKdV object
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
            Nind = length(U)/5;
            N = obj.ndof;
            R = zeros(Nind,1);
            J = zeros(5*Nind,1);
            v = [obj.dm2;obj.dm1;obj.d0;obj.dp1;obj.dp2];
            count = 0;
            if (obj.ind1==1)
                R(count+1,1) = -v'*U([4,5,1,2,3],1)-2*(U(2,1)-U(5,1))*(U(5,1)+U(1,1)+U(2,1))/(2*obj.dx);
                J(5*count+1:5*count+5,1) = [-obj.d0-(U(2,1)-U(5,1))/obj.dx;...
                                                               -obj.dp1-(U(1,1)+U(2,1))/obj.dx;...
                                                               -obj.dp2;...
                                                               -obj.dm2;...
                                                               -obj.dm1+(2*U(5,1)+U(1,1))/obj.dx];                                                   
               count = count + 1;
                                                           
            end
            if (obj.ind1==2 || obj.ind2==2)
                locInd = [5,1,2,3,4]+5*count*ones(1,5);
                R(count+1,1) = -v'*U(locInd,1)-2*(U(locInd(4),1)-U(locInd(2),1))*(U(locInd(2),1)+U(locInd(3),1)+U(locInd(4),1))/(2*obj.dx);
                J(5*count+1:5*count+5,1) = [-obj.dm1+(2*U(locInd(2),1)+U(locInd(3),1))/obj.dx;...
                                                              -obj.d0-(U(locInd(4),1)-U(locInd(2),1))/obj.dx;...
                                                               -obj.dp1-(U(locInd(3),1)+U(locInd(4),1))/obj.dx;...
                                                               -obj.dp2;...
                                                               -obj.dm2];
                count = count + 1;
            end
                     
            % middle block
            R(count+1:count+length(obj.Uind),1) = -obj.dm2*U(obj.Uind,1)-obj.dm1*U(obj.Uind+1,1)-obj.d0*U(obj.Uind+2,1)-obj.dp1*U(obj.Uind+3,1)-obj.dp2*U(obj.Uind+4,1)...
                                                                              -2*(U(obj.Uind+3,1)-U(obj.Uind+1,1)).*(U(obj.Uind+1,1)+U(obj.Uind+2,1)+U(obj.Uind+3,1))/(2*obj.dx);
            J(5*count+1:5*count+5*length(obj.Uind),1) = reshape([(-obj.dm2*ones(size(U(obj.Uind,1))))';...
                                                                                                     (-obj.dm1*ones(size(U(obj.Uind,1)))+(2*U(obj.Uind+1,1)+U(obj.Uind+2,1))/obj.dx)';...
                                                                                                     (-obj.d0*ones(size(U(obj.Uind,1)))-(U(obj.Uind+3,1)-U(obj.Uind+1,1))/obj.dx)';...
                                                                                                     (-obj.dp1*ones(size(U(obj.Uind,1)))-(U(obj.Uind+2,1)+2*U(obj.Uind+3,1))/obj.dx)';...
                                                                                                     (-obj.dp2*ones(size(U(obj.Uind,1))))'],5*length(obj.Uind),1);                                                                
                             
            count = count + length(obj.Uind);

            if (obj.indendm1 == N-1 || obj.indend == N-1)
                locInd = [2,3,4,5,1]+5*count*ones(1,5);  
                R(count+1,1) = -v'*U(locInd,1)-2*(U(locInd(4),1)-U(locInd(2),1))*(U(locInd(2),1)+U(locInd(3),1)+U(locInd(4),1))/(2*obj.dx);
                J(5*count+1:5*count+5,1) = [-obj.dp2;...
                                                              -obj.dm2;...
                                                              -obj.dm1+(2*U(locInd(2),1)+U(locInd(3),1))/obj.dx;...
                                                              -obj.d0-(U(locInd(4),1)-U(locInd(2),1))/obj.dx;...
                                                               -obj.dp1-(U(locInd(3),1)+U(locInd(4),1))/obj.dx];
                count = count + 1;
            end
            
            if (obj.indend == N)
                locInd = [3,4,5,1,2]+5*count*ones(1,5);  
                R(count+1,1) = -v'*U(locInd,1)-2*(U(locInd(4),1)-U(locInd(2),1))*(U(locInd(2),1)+U(locInd(3),1)+U(locInd(4),1))/(2*obj.dx);
                J(5*count+1:5*count+5,1) = [-obj.dp2;...
                                                               -obj.dp1-(U(locInd(3),1)+U(locInd(4),1))/obj.dx;...
                                                               -obj.dm2;...
                                                              -obj.dm1+(2*U(locInd(2),1)+U(locInd(3),1))/obj.dx;...
                                                              -obj.d0-(U(locInd(4),1)-U(locInd(2),1))/obj.dx];

            end




        end
        
        function  [R,J] = ResJacTimeDer(obj,U,~)
            

            R = U;
            J = speye(size(U,1));
  
        end
        
        function  [R,J] = ResJacGNATTimeDer(obj,U,~,~,~,~,jdiag)

            R = U;
            J = jdiag;
        
  
        end
        
        function  [R,J] = ResJacWithoutInput(obj,U,t)
            %This function returns the residual and jacobian without the
            %input term contributions.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj - OneDKdV object
            %
            %Outputs:
            %--------
            %Radj - residual without input term contribution
            %Jadj - jacobian without input term contribution
            %--------------------------------------------------------------
            
            [R,J] = obj.ResJac(U,t);
            
        end
        
        function  [ind] = node2ind(obj,nodenum) 
            %This function takes a node in the mesh and maps it to an index
            %in the full state vector
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - OneDKdV object
            %nodenum - scalar specifying node numbers
            %
            %Outputs:
            %--------
            %ind     - vector specifying the corresponding indices in the
            %          full state vector
            %--------------------------------------------------------------
            
            if  nodenum == obj.config.nNodes(1)
                ind = nan;
            else
                ind = nodenum ;
            end
        end
    
        function  [] = JacobianStructure(obj)
            %This function forms a boolean matrix whose true entries
            %correspond to nonzero entries in the jacobian
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - OneDKdVs object
            %
            %Outputs:
            %--------
            %There are no outputs.  It store a ndof x ndof boolean matrix
            %whose true values correspond to nonzero entries in the
            %jacobian in the class instance.
            %--------------------------------------------------------------
            
            N = obj.config.ndof;
            obj.Jstruct = spdiags(ones(N,5),-2:2,N,N);
            obj.Jstruct(1,N-1) = 1;
            obj.Jstruct(1,N) = 1;
            obj.Jstruct(2,N) = 1;
            obj.Jstruct(N-1,1) = 1;
            obj.Jstruct(N,1) = 1;
            obj.Jstruct(N,2) = 1;
        end 

        function  [probGNAT] = createCopy4GNAT(obj,gnat)
            %This function makes a copy of the current object and only
            %copies the source term and input matrix evaluated at indN
            %rows.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - OneDKdV object
            %gnat    - gnat object
            %
            %Outputs:
            %--------
            %probGNAT - new "reduced" OneDKdV object
            %--------------------------------------------------------------
            
            probGNAT = OneDKdV([],obj);
            
            temp = probGNAT.config.inFunc;
            temp2 = probGNAT.config.form;
            tempdm2 = probGNAT.dm2;
            tempdm1 = probGNAT.dm1;
            tempd0 = probGNAT.d0;
            tempdp1 = probGNAT.dp1;
            tempdp2 = probGNAT.dp2;
            N = probGNAT.config.ndof;
            probGNAT.config = [];
            probGNAT.config.inFunc = temp;
            probGNAT.config.form = temp2;
            probGNAT.mesh   = [];
            probGNAT.ndim = [];
            probGNAT.Jstruct = [];
            probGNAT.LinMat = [];
            probGNAT.D0 = [];
            probGNAT.ndof = N;
            probGNAT.dm2 = tempdm2;
            probGNAT.dm1 = tempdm1;
            probGNAT.d0 = tempd0;
            probGNAT.dp1 = tempdp1;
            probGNAT.dp2 = tempdp2;
                
           probGNAT.ind1 = gnat.sampleInd(1);
           probGNAT.ind2 = gnat.sampleInd(2);
           probGNAT.indendm1 = gnat.sampleInd(end-1);
           probGNAT.indend = gnat.sampleInd(end);
           if (probGNAT.ind1 == 1 && probGNAT.ind2 == 2)
               irstart_start = 3;
           elseif ((probGNAT.ind1 == 1 && probGNAT.ind2 ~= 2) || probGNAT.ind1 == 2)
                irstart_start = 2;
           else
              irstart_start = 1;
           end 
           if (probGNAT.indendm1 == N-1 && probGNAT.indend == N)
               irstart_end = length(gnat.irstart)-3;
           elseif ((probGNAT.indendm1 ~= N-1 && probGNAT.indend == N) || probGNAT.indend == N-1)
               irstart_end = length(gnat.irstart)-2;
           else
               irstart_end = length(gnat.irstart)-1;
           end
           probGNAT.Uind = gnat.irstart(irstart_start:irstart_end);
                     
           probGNAT.lengthJROW = gnat.irstart(end)-1;
         %  disp('TO DO end');

           probGNAT.dx = obj.dx;
                    
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
            %obj     - OneDKdV object
            %gnat    - gnat object
            %
            %Outputs:
            %--------
            %probGNAT - array of new "reduced" OneDKdV  object
            %--------------------------------------------------------------
           
            for kk = 1:gnat.nBases
                if isempty(gnat.sampleInd{kk,1})
                    continue;
                end
                probGNAT(kk) = OneDKdV([],obj);
                
                temp = probGNAT(kk).config.inFunc;
                temp2 = probGNAT(kk).config.form;
                tempdm2 = probGNAT(kk).dm2;
                tempdm1 = probGNAT(kk).dm1;
                tempd0 = probGNAT(kk).d0;
                tempdp1 = probGNAT(kk).dp1;
                tempdp2 = probGNAT(kk).dp2;
                N = probGNAT(kk).config.ndof;
                probGNAT(kk).config = [];
                probGNAT(kk).config.inFunc = temp;
                probGNAT(kk).config.form = temp2;
                probGNAT(kk).mesh   = [];
                probGNAT(kk).ndim = [];
                probGNAT(kk).Jstruct = [];
                probGNAT(kk).LinMat = [];
                probGNAT(kk).D0 = [];
                probGNAT(kk).ndof = N;
                probGNAT(kk).dm2 = tempdm2;
                probGNAT(kk).dm1 = tempdm1;
                probGNAT(kk).d0 = tempd0;
                probGNAT(kk).dp1 = tempdp1;
                probGNAT(kk).dp2 = tempdp2;                
                
                probGNAT(kk).ind1 = gnat.sampleInd{kk,1}(1);
                probGNAT(kk).ind2 = gnat.sampleInd{kk,1}(2);
                probGNAT(kk).indendm1 = gnat.sampleInd{kk,1}(end-1);
                probGNAT(kk).indend = gnat.sampleInd{kk,1}(end);
                if (probGNAT(kk).ind1 == 1 && probGNAT(kk).ind2 == 2)
                    irstart_start = 3;
                elseif ((probGNAT(kk).ind1 == 1 && probGNAT(kk).ind2 ~= 2) || probGNAT(kk).ind1 == 2)
                    irstart_start = 2;
                else
                    irstart_start = 1;
                end
                if (probGNAT(kk).indendm1 == N-1 && probGNAT(kk).indend == N)
                    irstart_end = length(gnat.irstart{kk,1})-3;
                elseif ((probGNAT(kk).indendm1 ~= N-1 && probGNAT(kk).indend == N) || probGNAT(kk).indend == N-1)
                    irstart_end = length(gnat.irstart{kk,1})-2;
                else
                    irstart_end = length(gnat.irstart{kk,1})-1;
                end
                probGNAT(kk).Uind = gnat.irstart{kk,1}(irstart_start:irstart_end);
                
                probGNAT(kk).lengthJROW = gnat.irstart{kk,1}(end)-1;

                probGNAT(kk).dx = obj.dx;            
                
                if length(obj.ic) > 1
                    probGNAT(kk).ic = obj.ic(unique(gnat.jrow{kk,1}),1);
                end
            end
        end
        
        
        function  [probGNAT] = updateIcCopy4GNATloc(obj,gnat,probGNAT)
            for kk = 1:gnat.nBases
                if length(obj.ic) > 1
                    probGNAT(kk).ic = obj.ic(unique(gnat.jrow{kk,1}),1);
                end
            end
        end
        
        function [] = modifyInitialCondition(obj,probobj)
            
            obj.ic = probobj.ic;            
            
        end
        
        function   [axOUT,lineobj] = problemPaperPlot(obj,modelobj,modelAuxobj,Tplots,pstr,axIN)
            %This function generates the 1D KdV Equation plot from Jean Frederic Gerbeau's paper.
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
        
        function  [X,Y,svPBC] = returnPosSV(obj,modelobj,modelAuxobj,~,flag,basenum)
            %This function takes a state vector from a simulation and
            %returns the state vector with the Dirichlet boundary
            %conditions.  A vector of positions (including those with PBCs)
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
            %svDBC       - state vector matrix with PBCs
            %--------------------------------------------------------------
            
            %Extract data (state vectors, cluster centers,
            %residual/jacobians, pod vectors)
            switch flag
                case 'pod'
                    svF = modelobj.phi{basenum,1};

                case 'res'
                    svF = modelobj.res;
                    

                case 'jac'
                    svF = modelobj.jac;

                case 'clusCenter'
                    svF = modelobj.clusCenter;

                case 'sv'
                    svF = modelobj.reconstructFullState(modelAuxobj);
                    
                    t = modelobj.time.T(1):modelobj.time.dt:modelobj.time.T(2);

            end
            n = size(svF,2);

            svPBC = zeros(obj.config.ndof+1,n);
            if ~isempty(svF)
                svPBC(1:end-1,:) = svF;
                svPBC(end,:) = svF(1,1);
            else
                error('Requested empty entry (most likely clusCenters or pod).');
            end
            
            X = obj.mesh.node;
            Y = [];
            
        end
        
    end
end
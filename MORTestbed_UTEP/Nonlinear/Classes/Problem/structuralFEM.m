classdef  structuralFEM < handle
    
    properties (SetAccess = private, GetAccess = public)
        config;
        
        %FEM structures
        msh;
        phys;
        integration;
        dbc;
        
        %Sparsity structures for Jacobian
        nnzeros; %number of nonzeros in Jacobian
        irow_crs; %row sparsity structure in compressed row format
        jcol_crs; %column sparsity structure in compressed row format
                  %(not usually necessary but needed in C++ evaluation of df)
        irow_coor;%row sparsity structure in coordinate format (for Matlab)
        jcol_coor;%column sparsity structure in coordinate format (for Matlab)
        
        ndof;
        
        prevsv;
        ic;
        ic_v;
        M;
        V;   %nel x 1 vector of volumes
        V_r; %reduced vector of volumes
        
        fExtNodal0;
        fExtNodal;
        fExtNodal_red;
        
        udbc; %DBC displacements
        adbc; %DBC accelerations
        
        b0;
        bdbc; %DBC body forces
        b;    %body forces (everywhere except DBCs)
        
        PrecompROM;
        p;
    end
    
    methods
        function  [obj] = structuralFEM(cfgobj,oldobj)
            %This is the constructor for structural FEM.  If
            %there is 1 input, the constructor will create the class
            %instance from the variables defined in this function.  If
            %there are 2 inputs, the first input will be unused and the
            %output will be a deep copy of the oldobj.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %cfgobj - CONFIG object
            %oldobj - structuralFEM object whose properties will be copied
            %         into the new object.
            %
            %Outputs:
            %--------
            %obj - structuralFEM object
            %--------------------------------------------------------------
            
            if nargin == 0
                %Keep everything at defaults if no input arguments
                return;
            end
            
            % Construct a new object based on a deep copy of an old object
            % of this class by copying properties over.
            if nargin == 2 && isa(oldobj,'structuralFEM')
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
            obj.config=cfgobj;
            
            [obj.ic,obj.ic_v,obj.fExtNodal0,obj.fExtNodal,...
               obj.b0,obj.b,obj.bdbc,obj.udbc,obj.adbc,obj.dbc,...
                  obj.msh,obj.phys,obj.integration] = feval(obj.config.altFile,obj.config.param{:});
              
            obj.config.setProperties('ndof',obj.msh.ndof);  
            %Create sparsity structure of tangent stiffness and mass matrices
            [obj.nnzeros,obj.irow_crs,obj.irow_coor,obj.jcol_coor] = createSparsityStructure(obj.msh);
            
            %Compute initial forces
            obj.setBCs(0);
            
            %Volume vector
            obj.volumeVector();
%             if obj.config.staticFlag
%                 return;
%             end
            
            %Compute Mass matrix
            Msp1d=zFEM_mxMass(obj.msh,obj.phys,obj.integration,obj.irow_crs,obj.jcol_coor,obj.nnzeros);
            obj.M=sparse(double(obj.irow_coor),double(obj.jcol_coor),Msp1d,double(obj.msh.ndof),double(obj.msh.ndof));
            
            obj.prevsv.v=obj.ic_v;
            fext = obj.fExtNodal(obj.fExtNodal0,0); %Nodal forces
            fext = fext + zFEM_mxBodyForces(obj.b(0),0,obj.msh,obj.phys,obj.integration,obj.dbc); 
            fext = fext - zFEM_mxMassVec(obj.msh,obj.phys,obj.integration,obj.dbc,obj.adbc(0));
            fint = zFEM_mxIntForces(obj.ic,0,obj.msh,obj.phys,obj.integration,obj.dbc);
            obj.prevsv.a = obj.M\( fext - fint );
        end
       
        function  [] = setICs(obj)
            obj.prevsv.v=zeros(size(obj.M,1),1);
            obj.prevsv.a=zeros(size(obj.M,1),1);
            
            disp('Debug statement in structuralFEM.setICs');
            return;
            obj.prevsv.v=obj.ic_v;
            fext = obj.fExtNodal(obj.fExtNodal0,0); %Nodal forces
            fext = fext + zFEM_mxBodyForces(obj.b(0),0,obj.msh,obj.phys,obj.integration,obj.dbc);
            fext = fext - zFEM_mxMassVec(obj.msh,obj.phys,obj.integration,obj.dbc,obj.adbc(0));
            fint = zFEM_mxIntForces(obj.ic,0,obj.msh,obj.phys,obj.integration,obj.dbc);
            obj.prevsv.a = obj.M\( fext - fint );
        end
        
        function  [] = setBCs(obj,t)
           obj.dbc.udbc = obj.udbc(t);
           obj.dbc.bdbc = obj.bdbc(t);
        end
        
        function  [] = setPrevSV(obj,iprev)
            obj.prevsv = iprev;
        end
        
        function  [R,J] =  ResJac(obj,U,t)
            % %Form jacobian from sparsity structure
            obj.setBCs(t);

            if nargout == 1
                R = zFEM_mxIntForces(U,t,obj.msh,obj.phys,obj.integration,obj.dbc);
                J=[];
            else
                [R,dfsp]=zFEM_mxIntForces(U,t,obj.msh,obj.phys,obj.integration,...
                    obj.dbc,obj.irow_crs,obj.jcol_coor,obj.nnzeros);
                J = sparse(double(obj.irow_coor),double(obj.jcol_coor),dfsp,double(obj.msh.ndof),double(obj.msh.ndof));
            end
            
            R = R - obj.fExtNodal(obj.fExtNodal0,t); %Nodal forces
            R = R - zFEM_mxBodyForces(obj.b(t),t,obj.msh,obj.phys,obj.integration,...
                                         obj.dbc,obj.irow_crs,obj.jcol_coor,obj.nnzeros);
            R = R + zFEM_mxMassVec(obj.msh,obj.phys,obj.integration,obj.dbc,obj.adbc(t));
        end
        
        function  [r,j] =  ResJacPrecomp(obj,y,t,basenum)
            
            [r,j] = obj.reconstRedForce(y,basenum);
            r = r - obj.fExtNodal(obj.PrecompROM(basenum).fExtNodal_red,t);
            r = r - obj.PrecompROM.fExtBody*t;
        end
        
        function  [r,drdy,drdp] =  ResJacPrecompSens(obj,y,t,basenum)
            n1=size(obj.PrecompROM.factors_lam,1);
            n2=size(obj.PrecompROM.factors_mu,1);
            n3=n2;
            
            fext = obj.fExtNodal(obj.PrecompROM(basenum).fExtNodal_red,t) ...
                 - obj.PrecompROM.fExtBody*t;
            
            [r,drdy] = obj.reconstRedForce(y,basenum);
            r = r - fext;
            drdp = zeros(size(r,1),n1+n2+n3);
            oldfactors=obj.PrecompROM.factors;
            for j = 1:n1
                obj.PrecompROM.factors{1}=obj.PrecompROM.factors_lam{j,1};
                obj.PrecompROM.factors{2}=obj.PrecompROM.factors_lam{j,2};
                obj.PrecompROM.factors{3}=obj.PrecompROM.factors_lam{j,3};
                obj.PrecompROM.factors{4}=obj.PrecompROM.factors_lam{j,4};
            
                drdp(:,j)=obj.reconstRedForce(y,basenum);
            end
            
            for j = 1:n2
                obj.PrecompROM.factors{1}=obj.PrecompROM.factors_mu{j,1};
                obj.PrecompROM.factors{2}=obj.PrecompROM.factors_mu{j,2};
                obj.PrecompROM.factors{3}=obj.PrecompROM.factors_mu{j,3};
                obj.PrecompROM.factors{4}=obj.PrecompROM.factors_mu{j,4};
                
                drdp(:,j+n1)=obj.reconstRedForce(y,basenum);
            end
            obj.PrecompROM.factors=oldfactors;
            drdp(:,n1+n2+1:end)=obj.PrecompROM.fExtBody_base*t;
        end
        
        function  [] = precomputeROM(obj,romobj,precompFlag,phi_lam,phi_mu,phi_rho)
            
            %Determine number of bases
            if iscell(romobj.phi)
                nBase=length(romobj.phi);
                phi = romobj.phi;
            else
                nBase=1;
                phi={romobj.phi};
            end
            
            if nBase == 1 && ~isempty(obj.prevsv)
                obj.prevsv.v = romobj.phi'*obj.prevsv.v;
                obj.prevsv.a = romobj.phi'*obj.prevsv.a;
            end
            
            %Loop over each basis and pre-compute terms
            for i = 1:nBase
                nY = size(phi{i},2);
                
                obj.PrecompROM(i).fExtNodal_red = phi{i}'*romobj.prob.fExtNodal0;
                
                if strcmpi(precompFlag,'nonparam')
%                     phi_lam = zeros(obj.msh.nel,1); phi_mu = zeros(obj.msh.nel,1);
%                     for jj = 1:obj.msh.nummat
%                        [phi_lam(obj.msh.mat==jj-1,1),phi_mu(obj.msh.mat==jj-1,1)] = Region.matphys2vec(jj-1,obj.phys); 
%                     end
%                     
%                     [matI,physI] = Region.vec2matphys(phi_lam,phi_mu);
%                     oldmsh = obj.msh;
%                     oldphys= obj.phys;
%                     physI(end).rho0=0;
%                     
%                     obj.updateMaterial(matI,physI);
                        
                    tStart=tic;
%                     [obj.PrecompROM(i).factors{1},obj.PrecompROM(i).factors{2},...
%                         obj.PrecompROM(i).factors{3},obj.PrecompROM(i).factors{4}] = ...
%                            zFEM_mxPrecompRom(phi{i}',zeros(romobj.prob.msh.ndof,1),nY,romobj.prob.msh,...
%                               romobj.prob.phys,romobj.prob.integration,romobj.prob.dbc);
                    
                    %Permanent change: use finite differences to compute 
                    [~,J] = romobj.prob.ResJac(zeros(romobj.prob.msh.ndof,1),1);
                    
                    obj.PrecompROM(i).factors{1} = zeros(nY,1);
                    obj.PrecompROM(i).factors{2} = phi{i}'*J*phi{i};
                    obj.PrecompROM(i).factors{3} = 0.5*red_hess_fd(romobj.prob,nY,phi{i},zeros(nY,1),1);
                    obj.PrecompROM(i).factors{4} = (1/6)*red_der3_fd(romobj.prob,nY,phi{i},zeros(nY,1),1);
                    obj.PrecompROM(i).time = toc(tStart);
                    
                    bF = romobj.prob.b(1);%obj.b0(~obj.dbc.loc);
                    oldbdbc = romobj.prob.dbc.bdbc;
                    obj.dbc.bdbc = romobj.prob.bdbc(1);%obj.b0(obj.dbc.loc);
                    t=1;
                    obj.PrecompROM(i).fExtBody = phi{i}'*zFEM_mxBodyForces(bF,t,romobj.prob.msh,...
                                                          romobj.prob.phys,romobj.prob.integration,...
                                                          romobj.prob.dbc,romobj.prob.irow_crs,romobj.prob.jcol_coor,...
                                                          romobj.prob.nnzeros);
                    obj.dbc.bdbc=oldbdbc;
%                     obj.updateMaterial(oldmsh.mat,oldphys);
                elseif strcmpi(precompFlag,'param')
                    %Reduce volume vector
                    obj.V = romobj.prob.V;
                    obj.reduceVolumeVector(phi_rho);
                    obj.V=[];
            
                    n_lam = size(phi_lam,2);
                    n_mu  = size(phi_mu ,2);
                    n_rho = size(phi_rho,2);
                    
                    if ~(n_lam == n_mu && n_mu ==n_rho)
                        error('Material bases must have same number of columns');
                    end
                    n = n_lam; % = n_mu = n_rho
                    
                    obj.PrecompROM(i).n_lam = n_lam;
                    obj.PrecompROM(i).n_mu  = n_mu;
                    obj.PrecompROM(i).n_rho = n_rho;
%                     n_min = min(n_lam,n_mu);
%                     [n_max,ind] = max([n_lam,n_mu]);
                    
                    obj.PrecompROM(i).factors_lam = cell(n,4);
                    obj.PrecompROM(i).factors_mu = cell(n,4);
                    
                    obj.PrecompROM(i).fExtBody_base = zeros(nY,n);
                    
                    oldmat = romobj.prob.msh.mat;
                    oldphys= romobj.prob.phys;
                    oldbdbc = romobj.prob.dbc.bdbc;

                    tStart=tic;
                    for j = 1:n
                        [matI,physI] = Region.vec2matphys(phi_lam(:,j),phi_mu(:,j),phi_rho(:,j));
                        
                        romobj.prob.updateMaterial(matI,physI);
                        
%                         [obj.PrecompROM(i).factors_lam{j,1},obj.PrecompROM(i).factors_lam{j,2},...
%                          obj.PrecompROM(i).factors_lam{j,3},obj.PrecompROM(i).factors_lam{j,4},...
%                          obj.PrecompROM(i).factors_mu{j,1}, obj.PrecompROM(i).factors_mu{j,2},...
%                          obj.PrecompROM(i).factors_mu{j,3}, obj.PrecompROM(i).factors_mu{j,4}] = ...
%                             zFEM_mxPrecompRom(phi{i}',zeros(obj.msh.ndof,1),nY,obj.msh,...
%                             physI,obj.integration,obj.dbc);

                        %disp('DEBUG STATEMENTs in precomputeROM!');
                        obj.PrecompROM(i).factors_lam{j,1} = zeros(nY,1);
                        [~,J] = romobj.prob.ResJac(zeros(romobj.prob.msh.ndof,1),1);
                        obj.PrecompROM(i).factors_lam{j,2} = phi{i}'*J*phi{i}; clear J;
                        obj.PrecompROM(i).factors_lam{j,3} = 0.5*red_hess_fd(romobj.prob,nY,phi{i},zeros(nY,1),1);
                        obj.PrecompROM(i).factors_lam{j,4} = (1/6)*red_der3_fd(romobj.prob,nY,phi{i},zeros(nY,1),1);

%                         disp('DEBUG STATEMENTs in precomputeROM!');
                        obj.PrecompROM(i).factors_mu{j,1} = zeros(nY,1);
                        obj.PrecompROM(i).factors_mu{j,2} = zeros(nY,nY);
                        obj.PrecompROM(i).factors_mu{j,3} = zeros(nY,nY,nY);
                        obj.PrecompROM(i).factors_mu{j,4} = zeros(nY,nY,nY,nY);
                        
                        bF = romobj.prob.b(1);%obj.b0(~obj.dbc.loc);
                        romobj.prob.dbc.bdbc = romobj.prob.bdbc(1);%obj.b0(obj.dbc.loc);
                        t=1;
                        obj.PrecompROM(i).fExtBody_base(:,j) = phi{i}'*zFEM_mxBodyForces(bF,t,romobj.prob.msh,...
                                                          romobj.prob.phys,romobj.prob.integration,...
                                                          romobj.prob.dbc,romobj.prob.irow_crs,romobj.prob.jcol_coor,...
                                                          romobj.prob.nnzeros);
                        
                    end
                    romobj.prob.updateMaterial(oldmat,oldphys);
                    romobj.prob.dbc.bdbc=oldbdbc;
                    disp('Hard coded body force treatment for precomputed ROM...cannot handle time dependent body force');
                    
%                     for j = n_min+1:n_max
%                         if ind == 1
%                             [matI,physI] = Region.vec2matphys(phi_lam(:,j),zeros(obj.msh.nel,1));
%                             oldmsh = obj.msh;
%                             oldphys= obj.phys;
%                             
%                             obj.updateMaterial(matI,physI);
%                             
%                             [obj.PrecompROM(i).factors_lam{j,1},obj.PrecompROM(i).factors_lam{j,2},...
%                                 obj.PrecompROM(i).factors_lam{j,3},obj.PrecompROM(i).factors_lam{j,4}] = ...
%                                 zFEM_mxPrecompRom(phi{i}',zeros(obj.msh.ndof,1),nY,obj.msh,...
%                                 obj.phys,obj.integration,obj.dbc,ind);
%                             
%                             obj.updateMaterial(oldmsh.mat,oldphys);
%                         elseif ind == 2
%                             [matI,physI] = Region.vec2matphys(zeros(obj.msh.nel,1),phi_mu(:,j));
%                             oldmsh = obj.msh;
%                             oldphys= obj.phys;
%                             
%                             obj.updateMaterial(matI,physI);
%                             
%                             [obj.PrecompROM(i).factors_mu{j,1},obj.PrecompROM(i).factors_mu{j,2},...
%                                 obj.PrecompROM(i).factors_mu{j,3},obj.PrecompROM(i).factors_mu{j,4}] = ...
%                                 zFEM_mxPrecompRom(phi{i}',zeros(obj.msh.ndof,1),nY,obj.msh,...
%                                 obj.phys,obj.integration,obj.dbc,ind);
%                             
%                             obj.updateMaterial(oldmsh.mat,oldphys);
%                         end
%                     end
                    obj.PrecompROM(i).time = toc(tStart);
                end
            end
        end
        
        function  [f,df] = reconstRedForce(obj,y,basenum)
            
            nY = length(y);
            %Reduced internal force
            f  = obj.PrecompROM(basenum).factors{1} + obj.PrecompROM(basenum).factors{2}'*y;
            for r1 = 1:nY
                f(r1) = f(r1) + y'*obj.PrecompROM(basenum).factors{3}(:,:,r1)*y;
                for r2 = 1:nY
                    f(r1) = f(r1) + y(r2)*(y'*obj.PrecompROM(basenum).factors{4}(:,:,r2,r1)*y);
                end
            end
            
            if nargout == 2
                %disp('DEBUG STATEMENT in reconstRedForce!');
                df = obj.PrecompROM(basenum).factors{2}';
                for r3 = 1:nY
                    df = df + 2*y(r3)*squeeze(obj.PrecompROM(basenum).factors{3}(r3,:,:));
                    for r4 = 1:nY
                        df = df + 3*(squeeze(obj.PrecompROM(basenum).factors{4}(r4,r3,:,:)))*y(r3)*y(r4);
                    end
                end

%                 %Reduced tangent stiffness
%                 df = obj.PrecompROM(basenum).factors{2}';
%                 for r3 = 1:nY
%                     df = df + y(r3)*(squeeze(obj.PrecompROM(basenum).factors{3}(r3,:,:))' + squeeze(obj.PrecompROM(basenum).factors{3}(:,r3,:))');
%                     for r4 = 1:nY
%                         df = df + (squeeze(obj.PrecompROM(basenum).factors{4}(r4,r3,:,:))' + ...
%                                    squeeze(obj.PrecompROM(basenum).factors{4}(r4,:,r3,:))' + ...
%                                    squeeze(obj.PrecompROM(basenum).factors{4}(:,r4,r3,:))')*y(r3)*y(r4);
%                     end
%                 end
            end
        end
                
        function  [] = updateMaterial(obj,mat,phys1)
            obj.msh.mat=uint32(mat);
            obj.msh.nummat=uint32(length(unique(mat)));
            obj.phys=phys1;
        end
        
        function  [] = PrecompFromRedMatPrecomp(obj,lam_r,mu_r,rho_r)
            
            obj.p = [lam_r;mu_r;rho_r];
            
            obj.PrecompROM.factors{1}=obj.PrecompROM.factors_lam{1,1}*lam_r(1);
            obj.PrecompROM.factors{2}=obj.PrecompROM.factors_lam{1,2}*lam_r(1);
            obj.PrecompROM.factors{3}=obj.PrecompROM.factors_lam{1,3}*lam_r(1);
            obj.PrecompROM.factors{4}=obj.PrecompROM.factors_lam{1,4}*lam_r(1);
            
            for i = 2:size(obj.PrecompROM.factors_lam,1)
                obj.PrecompROM.factors{1}=obj.PrecompROM.factors{1}+obj.PrecompROM.factors_lam{i,1}*lam_r(i);
                obj.PrecompROM.factors{2}=obj.PrecompROM.factors{2}+obj.PrecompROM.factors_lam{i,2}*lam_r(i);
                obj.PrecompROM.factors{3}=obj.PrecompROM.factors{3}+obj.PrecompROM.factors_lam{i,3}*lam_r(i);
                obj.PrecompROM.factors{4}=obj.PrecompROM.factors{4}+obj.PrecompROM.factors_lam{i,4}*lam_r(i);
            end
            
            for i = 1:size(obj.PrecompROM.factors_mu,1)
                obj.PrecompROM.factors{1}=obj.PrecompROM.factors{1}+obj.PrecompROM.factors_mu{i,1}*mu_r(i);
                obj.PrecompROM.factors{2}=obj.PrecompROM.factors{2}+obj.PrecompROM.factors_mu{i,2}*mu_r(i);
                obj.PrecompROM.factors{3}=obj.PrecompROM.factors{3}+obj.PrecompROM.factors_mu{i,3}*mu_r(i);
                obj.PrecompROM.factors{4}=obj.PrecompROM.factors{4}+obj.PrecompROM.factors_mu{i,4}*mu_r(i);
            end
            
            obj.PrecompROM.fExtBody = obj.PrecompROM.fExtBody_base*rho_r;
        end
        
        function  [] = clearMatPrecompFact(obj)
           obj.PrecompROM.factors_lam=[];
           obj.PrecompROM.factors_mu =[];
        end
        
        function  [] = volumeVector(obj)
            
            obj.V = zeros(obj.msh.nel,1);
            for j=1:obj.msh.nel
                MM = [obj.msh.X(:,obj.msh.IEN(1,j)+1)',1;...
                     obj.msh.X(:,obj.msh.IEN(2,j)+1)',1;...
                     obj.msh.X(:,obj.msh.IEN(3,j)+1)',1;...
                     obj.msh.X(:,obj.msh.IEN(4,j)+1)',1];
                obj.V(j) = (1/6)*abs(det(MM));
            end
        end
        
        function  [] = reduceVolumeVector(obj,phi_rho)
            obj.PrecompROM.V_r = phi_rho'*obj.V;
        end
        
        %Postprocessing
        function  [] = export2AeroS(obj,fname,fomobj)
            %This assumes tetrahedra!!!!
            %Cannot handle time dependent loads and forces
            %Assumes STATIC!
            
            fid = fopen(fname,'w+');
            
            ind=regexp(fname,'/');
            newfname=fname;
            if ~isempty(ind), newfname=fname(ind(end)+1:end); end;
            
            fprintf(fid,'STATIC\n');
            fprintf(fid,'sparse\n');
            fprintf(fid,'*\n');
            fprintf(fid,'NONLINEAR\n');
            fprintf(fid,'maxit %i\n',fomobj.newt.maxIter);
            fprintf(fid,'nltol %e\n',fomobj.newt.eps(1));
            fprintf(fid,'dlambda %f %f\n',fomobj.time.dt,fomobj.time.T(2));
            fprintf(fid,'*\n');
            fprintf(fid,'OUTPUT\n');
            fprintf(fid,'gdisplac 20 16 "%s_disp" %i\n',newfname,1);
            fprintf(fid,'stressvm 20 16 "%s_stressvm" %i\n',newfname,1);
            fprintf(fid,'stressp1 20 16 "%s_stressp1" %i\n',newfname,1);
            fprintf(fid,'*\n');
            
            fclose(fid);
            
            distmesh2aeros(fname,obj.msh,obj.fExtNodal0,...
                obj.dbc.loc,obj.udbc(0),obj.phys,'a');
        end
        
        function  [] = AeroSTop2zFEM(obj,fname)
           %This function should accept an AeroS top file and convert it into 
           %zFEM structures: msh.X, msh.IEN, msh.fExtNodal0, msh.udbc
           
           %Also need to make a separate function that does the same thing,
           %except it takes AeroS input files (not top files)
           
%             fid = fopen(fname,'r+');
%             
%             while 1
%                tline = fgetl(fid);
%                if ~ischar(tline), break, end;
%                
%                if ~isempty(strfind(tline),'Nodes nodeset')
%                    while 1
%                        tline = fgetl(fid);
%                        if ~ischar(tline), break, end;
%                        if ~isempty(strfind(tline),'Elements elemset using nodeset'), break, end;
%                        
%                        tline=deblank(tline);
%                        
%                        ind = regexp(tline,'\s');
%                        ind = [1,ind,length(ind)];
%                        
%                        nodenum = num2str(ind(1):ind(2));
%                        for i = 2:length(ind)-1
%                            node(i-1) = num2str(ind(i):ind(i+1));
%                        end
%                    end
%                end
%                
%                
%             end
%             
%             fclose(fid);
            
        end
        
        function  [] = ml2xpost_disp(obj,fname,model)
            
            fid = fopen(fname,'w');
            fprintf(fid,'Vector DISP under NLStatic for nodeset\n');
            fprintf(fid,'%i\n',obj.msh.np);
            
            u=zeros(obj.msh.nsd*obj.msh.np,1);
            for i = 1:model.time.nstep+1
                t=model.time.T(1)+(i-1)*model.time.dt;
                
                fprintf(fid,'   %5.4E\n',t);
                
                u(~obj.dbc.loc(:)) = model.sv(:,i);
                u(obj.dbc.loc(:))  = obj.udbc(t);
                u = reshape(u,obj.msh.nsd,obj.msh.np);
                for j = 1:obj.msh.np
                    fprintf(fid,' %5.4E   %5.4E   %5.4E\n',u(1,j),u(2,j),u(3,j));
                end
                u = reshape(u,obj.msh.np,obj.msh.nsd);
            end
            fclose(fid);
        end
        
        function  [] = ml2paraview_disp(obj,fname,model)
            
            for step = 1:size(model.sv,2)
                fid = fopen([fname,num2str(step),'.vtk'],'w');
                fprintf(fid,'# vtk DataFile Version 3.1\n');
                fprintf(fid,'zFEM structural displacement output\n');
                fprintf(fid,'ASCII\n\n');
                
                fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');
                fprintf(fid,'POINTS \t %i \t FLOAT\n',obj.msh.np);
                for i = 1:obj.msh.np
                    fprintf(fid,'%7.6f\t%7.6f\t%7.6f\n',obj.msh.X(1,i),obj.msh.X(2,i),obj.msh.X(3,i));
                end
                fprintf(fid,'\n');
                
                fprintf(fid,'CELLS \t %i \t %i\n',obj.msh.nel,5*obj.msh.nel);
                for i = 1:obj.msh.nel
                    fprintf(fid,'4\t%i\t%i\t%i\t%i\n',obj.msh.IEN(1,i),obj.msh.IEN(2,i),obj.msh.IEN(3,i),obj.msh.IEN(4,i));
                end
                fprintf(fid,'\n');
                
                fprintf(fid,'CELL_TYPES \t %i\n',obj.msh.nel);
                for i = 1:obj.msh.nel
                    fprintf(fid,'10\n');
                    
                end
                fprintf(fid,'\n');
                
                fprintf(fid,'POINT_DATA %i\n',obj.msh.np);
                u=zeros(obj.msh.nsd*obj.msh.np,1);
                t=model.time.T(1)+(step-1)*model.time.dt;

                fprintf(fid,'VECTORS DISP FLOAT\n');
                u(~obj.dbc.loc(:)) = model.sv(:,step);
                u(obj.dbc.loc(:))  = obj.udbc(t);
                u = reshape(u,obj.msh.nsd,obj.msh.np);
                for j = 1:obj.msh.np
                    fprintf(fid,'%20.16f\t%20.16f\t%20.16f\n',u(1,j),u(2,j),u(3,j));
                end
                fprintf(fid,'\n');
                
                fclose(fid);
            end
        end

        function  [] = ml2paraview_mat(obj,fname,X_mat,U)
            
            for step = 1:size(X_mat,2)
                if size(X_mat,2) == 1
                    fid = fopen([fname,'.vtk'],'w');
                else
                    fid = fopen([fname,num2str(step),'.vtk'],'w');
                end
                fprintf(fid,'# vtk DataFile Version 3.1\n');
                fprintf(fid,'zFEM structural material output\n');
                fprintf(fid,'ASCII\n\n');
                
                fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');
                fprintf(fid,'POINTS \t %i \t FLOAT\n',obj.msh.np);
                for i = 1:obj.msh.np
                    fprintf(fid,'%7.6f\t%7.6f\t%7.6f\n',obj.msh.X(1,i),obj.msh.X(2,i),obj.msh.X(3,i));
                end
                fprintf(fid,'\n');
                
                numcells=sum(X_mat(:,step)~=0);
                fprintf(fid,'CELLS \t %i \t %i\n',numcells,5*numcells);
                for i = 1:obj.msh.nel
                    if X_mat(i,step) ~= 0
                        fprintf(fid,'4\t%i\t%i\t%i\t%i\n',obj.msh.IEN(1,i),obj.msh.IEN(2,i),obj.msh.IEN(3,i),obj.msh.IEN(4,i));
                    end
                end
                fprintf(fid,'\n');
                
                fprintf(fid,'CELL_TYPES \t %i\n',numcells);
                for i = 1:numcells
                    fprintf(fid,'10\n');
                end
                fprintf(fid,'\n');
                
                if nargin == 4 && ~isempty(U)
                    fprintf(fid,'POINT_DATA %i\n',obj.msh.np);
                    u=zeros(obj.msh.nsd*obj.msh.np,1);
                    t=1;
                    
                    fprintf(fid,'VECTORS DISP FLOAT\n');
                    u(~obj.dbc.loc(:)) = U;
                    u(obj.dbc.loc(:))  = obj.udbc(t);
                    u = reshape(u,obj.msh.nsd,obj.msh.np);
                    for j = 1:obj.msh.np
                        fprintf(fid,'%20.16f\t%20.16f\t%20.16f\n',u(1,j),u(2,j),u(3,j));
                    end
                    fprintf(fid,'\n');
                end
                
                fprintf(fid,'CELL_DATA %i\n',numcells);
                fprintf(fid,'SCALARS MAT FLOAT\n');
                fprintf(fid,'LOOKUP_TABLE default\n');
                for j = 1:obj.msh.nel
                    if X_mat(j,step) ~= 0
                        fprintf(fid,'%20.16f\n',X_mat(j,step));
                    end
                end
                fprintf(fid,'\n');
                fclose(fid);
            end
        end

        function  [] = ml2paraview_matdist(obj,fname,X_mat)
            
            for step = 1:size(X_mat,2)
                fid = fopen([fname,num2str(step),'.vtk'],'w');
                fprintf(fid,'# vtk DataFile Version 3.1\n');
                fprintf(fid,'zFEM structural material distribution output\n');
                fprintf(fid,'ASCII\n\n');
                
                fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');
                fprintf(fid,'POINTS \t %i \t FLOAT\n',obj.msh.np);
                for i = 1:obj.msh.np
                    fprintf(fid,'%7.6f\t%7.6f\t%7.6f\n',obj.msh.X(1,i),obj.msh.X(2,i),obj.msh.X(3,i));
                end
                fprintf(fid,'\n');
                
                fprintf(fid,'CELLS \t %i \t %i\n',obj.msh.nel,5*obj.msh.nel);
                for i = 1:obj.msh.nel
                    fprintf(fid,'4\t%i\t%i\t%i\t%i\n',obj.msh.IEN(1,i),obj.msh.IEN(2,i),obj.msh.IEN(3,i),obj.msh.IEN(4,i));
                end
                fprintf(fid,'\n');
                
                fprintf(fid,'CELL_TYPES \t %i\n',obj.msh.nel);
                for i = 1:obj.msh.nel
                    fprintf(fid,'10\n');
                end
                fprintf(fid,'\n');
                
                fprintf(fid,'CELL_DATA %i\n',obj.msh.nel);
                fprintf(fid,'SCALARS MAT FLOAT\n');
                fprintf(fid,'LOOKUP_TABLE default\n');
                for j = 1:obj.msh.nel
                    fprintf(fid,'%20.16f\n',X_mat(j,step));
                end
                fprintf(fid,'\n');
                fclose(fid);
            end
        end

        function  [] = ml2xpost_stress(obj,fname,model)
            
            fid = fopen(fname,'w');
            fprintf(fid,'Scalar VONMISES under NLStatic for nodeset\n');
            fprintf(fid,'%i\n',obj.msh.np);
            
            for i = 1:model.time.nstep+1
                t=model.time.T(1)+(i-1)*model.time.dt;
                
                fprintf(fid,'   %5.4E\n',t);
                
                sigma = zFEM_mxStresses(model.sv(:,i),t,obj.msh,obj.phys,obj.dbc);
                for j = 1:obj.msh.np
                    fprintf(fid,'  %5.4E\n',sigma(j));
                end
            end
            fclose(fid);
        end
        
        function  [] = createtop(obj,fname)
           
            fid = fopen(fname,'w');
            fprintf(fid,'Nodes nodeset\n');
            for i = 1:obj.msh.np
                fprintf(fid,'%i\t%7.6f\t%7.6f\t%7.6f\n',i,obj.msh.X(1,i),obj.msh.X(2,i),obj.msh.X(3,i));
            end
            fprintf(fid,'Elements elemset using nodeset\n');
            for i = 1:obj.msh.nel
                fprintf(fid,'\t%i\t 123\t%i\t%i\t%i\t%i\n',i,obj.msh.IEN(1,i)+1,obj.msh.IEN(2,i)+1,obj.msh.IEN(3,i)+1,obj.msh.IEN(4,i)+1);
            end
            fprintf(fid,'Pattern default using elemset\n');
            fprintf(fid,'SDBoundary bcondset using nodeset\n');
            cnt=0; u=obj.udbc(0);
            for i = 1:obj.msh.np
                for j = 1:obj.msh.nsd
                    if obj.dbc.loc(j,i)
                        cnt=cnt+1;
                        fprintf(fid,'%i \t %i \t %7.6f\n',i,j,u(cnt));
                    end
                end
            end
            fprintf(fid,'SFBoundary bcondset using nodeset\n');
            cnt=0; f=obj.fExtNodal(obj.fExtNodal0,0);
            for i = 1:obj.msh.np
                for j = 1:obj.msh.nsd
                    if ~obj.dbc.loc(j,i)
                        cnt=cnt+1;
                        if f(cnt)~=0, fprintf(fid,'%i \t %i \t %7.6f\n',i,j,f(cnt)); end;
                    end
                end
            end
        end
        
        function  [probROM] = createReducedCopy(obj,rom,precompflag,philam,phimu,phirho)
            %This function makes a copy of the current object and only
            %copies the source term and input matrix evaluated at indN
            %rows.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj    - structuralFEM object
            %rom    - ROM object
            %
            %Outputs:
            %--------
            %probROM - new "reduced" structuralFEM object
            %--------------------------------------------------------------
            
            probROM = structuralFEM();
            
            probROM.ic = rom.phi'*rom.prob.ic;
            %probROM.ic_v=rom.phi'*rom.prob.ic_v;
            probROM.M  = rom.phi'*rom.prob.M*rom.phi;
            probROM.udbc = rom.prob.udbc;
            probROM.bdbc = rom.prob.bdbc;
            probROM.fExtNodal = rom.prob.fExtNodal;
            
            probROM.precomputeROM(rom,precompflag,philam,phimu,phirho);
            
            probROM.msh = [];
            probROM.integration = [];
            probROM.dbc = [];
            
            probROM.irow_coor = [];
            probROM.jcol_coor = [];
            probROM.irow_crs = [];
            probROM.nnzeros = [];
            
            
            probROM.fExtNodal0=[];
            %probROM.udbc = [];
            probROM.adbc = [];
            %probROM.bdbc = [];
            %probROM.b = [];
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

        %Optimization
        function  [f,dfdy,dfdp] = RedCompliance(obj,y,p)
            n1 = obj.PrecompROM.n_lam; n2=obj.PrecompROM.n_mu; n3=obj.PrecompROM.n_rho;
            rho_r = p(n1+n2+1:end);
            
            f    = (obj.PrecompROM.fExtNodal_red + obj.PrecompROM.fExtBody_base*rho_r)'*y;
            dfdy = (obj.PrecompROM.fExtNodal_red + obj.PrecompROM.fExtBody_base*rho_r);
            dfdp = [zeros(n1+n2,1);...
                    obj.PrecompROM.fExtBody_base'*y];
        end
        
        function  [f,dfdy,dfdp] = RedWeight(obj,y,p)
            n1 = obj.PrecompROM.n_lam; n2=obj.PrecompROM.n_mu; n3=obj.PrecompROM.n_rho;
            rho_r = p(n1+n2+1:end);
            
            f    = 9.81*(obj.PrecompROM.V_r'*rho_r);
            dfdy = zeros(size(y));
            dfdp = [zeros(n1+n2,1);9.81*obj.PrecompROM.V_r];
        end
        
        function  [c,ceq,dc,dceq] = ExactHole(obj,p)
            n1 = obj.PrecompROM.n_lam; n2=obj.PrecompROM.n_mu; n3=obj.PrecompROM.n_rho;
            lam_r = p(1:n1); mu_r = p(n1+1:n1+n2); rho_r = p(n1+n2+1:end);
            c=[]; dc=[];
            
            ceq = lam_r(2:end).*(lam_r(2:end) - lam_r(1));
            ceq = [ceq; mu_r(2:end).*(mu_r(2:end) - mu_r(1))];
            ceq = [ceq; rho_r(2:end).*(rho_r(2:end) - rho_r(1))];
            
            n=length(p);
            
            vals_lam = zeros(2*n1-2,1);
            vals_lam(1:n1-1)=-lam_r(2:end);
            vals_lam(n1:end)=2*lam_r(2:end)-lam_r(1);
            
            vals_mu = zeros(2*n2-2,1);
            vals_mu(1:n2-1)=-mu_r(2:end);
            vals_mu(n2:end)=2*mu_r(2:end)-mu_r(1);
            
            vals_rho = zeros(2*n3-2,1);
            vals_rho(1:n3-1)=-rho_r(2:end);
            vals_rho(n3:end)=2*rho_r(2:end)-rho_r(1);
            
            rows = [1:n1-1,1:n1-1]'; cols = [ones(1,n1-1),2:n1]';
            rows=[rows;rows+n1-1;rows+n1+n2-2]; cols = [cols;cols+n1;cols+n1+n2];
            vals = [vals_lam;vals_mu;vals_rho];
            
            dceq = sparse(rows,cols,vals)';
        end
        
        function  [A,b] = MatchMatParams(obj)
            %This makes sure the lam, mu, and rho
            %distributions are the same
            n1 = obj.PrecompROM.n_lam; n2=obj.PrecompROM.n_mu; n3=obj.PrecompROM.n_rho;

            irows = [1:n1,n1+1:n1+n2,1:n1,n1+1:n1+n2]';
            icols = [1:n1,1:n1,n1+1:n1+n2,n1+n2+1:n1+n2+n3]';
            vals  = [ones(2*n1,1);-ones(2*n2,1)];
            A = sparse(irows,icols,vals);
            b=zeros(2*n1,1);
        end
        
        function  [A,b] = PositiveMat(obj)
            n1 = obj.PrecompROM.n_lam; n2=obj.PrecompROM.n_mu; n3=obj.PrecompROM.n_rho;

            vals_lam = zeros(2*n1-2,1);
            vals_lam(1:n1-1)=-1;
            vals_lam(n1:end)=1;
            
            vals_mu = zeros(2*n2-2,1);
            vals_mu(1:n2-1)=-1;
            vals_mu(n2:end)=1;
            
            vals_rho = zeros(2*n3-2,1);
            vals_rho(1:n3-1)=-1;
            vals_rho(n3:end)=1;
            
            rows = [1:n1-1,1:n1-1]'; cols = [ones(1,n1-1),2:n1]';
            rows=[rows;rows+n1-1;rows+n1+n2-2]; cols = [cols;cols+n1;cols+n1+n2];
            vals = [vals_lam;vals_mu;vals_rho];
            A = sparse(rows,cols,vals);
            b=zeros(3*(n1-1),1);
        end
        
        function  [c,ceq,dc,dceq] = RedMaxDisplacement(obj,p,rom,eps,phi_rowsum,flag)
            n1 = obj.PrecompROM.n_lam; n2=obj.PrecompROM.n_mu; n3=obj.PrecompROM.n_rho;
            obj.PrecompFromRedMatPrecomp(p(1:n1),p(n1+1:n1+n2),p(n1+n2+1:end));
            rom.executeModel();
            
            y = rom.sv(:,end);
            ceq=[]; dceq=[];
            [R,dRdy,dRdp] = obj.ResJacPrecompSens(y,1,1);
            
            dydp = -dRdy\dRdp;
            
            c1 = [y - eps*phi_rowsum;-eps*phi_rowsum-y];
            dc1= [dydp;...
                 -dydp];
            dc1=dc1';
            if nargin < 6 || isempty(flag)
                c = c1;
                dc = dc1;
                return;
            elseif strcmpi(flag,'sizes')
                c = length(c1);
                ceq=0;
            elseif strcmpi(flag,'constraint')
                c=c1;
            elseif strcmpi(flag,'jacobian')
                c = sparse(dc1);
            elseif strcmpi(flag,'jacobianstructure')
                c = sparse(ones(length(c1),length(p)));
            end
        end
        
        function  [c,ceq,dc,dceq] = MaxDisplacement(obj,p,rom,ind,eps,phi,flag)
           n1 = obj.PrecompROM.n_lam; n2=obj.PrecompROM.n_mu;
           obj.PrecompFromRedMatPrecomp(p(1:n1),p(n1+1:n1+n2),p(n1+n2+1:end));
           rom.executeModel();
           
           y = rom.sv(:,end);
           ceq=[]; dceq=[];
           [R,dRdy,dRdp] = obj.ResJacPrecompSens(y,1,1);
           
           dydp = -dRdy\dRdp;
           
           c1 = [phi(ind,:)*y - eps;-eps-phi(ind,:)*y];
           dc1= [phi(ind,:)*dydp;...
               -phi(ind,:)*dydp];
           dc1=dc1';
           if nargin < 7 || isempty(flag)
               c = c1;
               dc = dc1;
               return;
           elseif strcmpi(flag,'sizes')
               c = length(c1);
               ceq=0;
           elseif strcmpi(flag,'constraint')
               c=c1;
           elseif strcmpi(flag,'jacobian')
               c = sparse(dc1);
           elseif strcmpi(flag,'jacobianstructure')
               c = sparse(ones(length(c1),length(p)));
           end
        end
    end
        
    methods (Static)
        function  [] = plotDeformed(model,U,geom)
            
            p=model.prob.msh.X';
            
            if nargin < 3 || isempty(geom)
                t=model.prob.msh.IEN'+1;
            else
                if strcmpi(class(geom),'Region')
                    t = zeros(0,4);
                    for i = 1:length(geom)
                        t = union(t,geom(i).elem,'rows','stable');
                    end
                else
                    t = geom;
                end
            end
            
            p_def = p';
            p_def(~model.prob.dbc.loc) = p_def(~model.prob.dbc.loc) + U;
            p_def=p_def';

            simpplot(p_def,t);
        end
        
        function  [u,uDOF] = readGDISPLAC(fname,steps,dbc)
            
            if nargin < 2 || isempty(steps), steps='all'; end;
            
            text=fileread(fname);
            lines = regexp(text,'\n','split');
            np=str2double(lines{2});
            lines(:,[1:2,length(lines)]) = [];
            lines(:,1:(np+1):end) = [];
            
            totsteps = length(lines)/np;
            if strcmpi(steps,'all'), steps=1:totsteps; end;
            if strcmpi(steps,'last'), steps=totsteps; end;
            
            lines=lines';
            u=zeros(3,np,length(steps));
            uDOF=[];
            cnt=0;
            for i=steps
                cnt=cnt+1;
                u(:,:,cnt)=cell2mat(cellfun(@str2num,lines(np*(i-1)+1:np*i),'UniformOutput',false))';
            end
            clear lines;
            
            if nargin == 3
                uDOF=u(repmat(~dbc.loc,[1,1,length(steps)]));
                uDOF=reshape(uDOF,sum(sum(~dbc.loc)),length(steps));
            end
            u=reshape(u,3*np,length(steps));
        end
    end
end
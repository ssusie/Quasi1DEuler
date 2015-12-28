classdef Region < handle
    
    properties (SetAccess=private, GetAccess=public)
        name;
        
        elem;
        numelem;
        
        filters;
        nfilters;
        
        plotexpr;
    end
    
    properties (Hidden=true, SetAccess=private, GetAccess=public)
        xlim;
        ylim;
        zlim;
        xlimX;
        ylimX;
        zlimX;
        
        hFig;
    end
    
    methods
        function [obj] = Region(regionname,p,t,varargin)
            
            if nargin == 0, return; end
            
            obj.name = regionname;
            
            %Extract filters
            obj.nfilters = length(varargin);
            for i = 1:obj.nfilters
                curr_filt = varargin{i};
                
                obj.filters(i).type = curr_filt{1};
                
                switch obj.filters(i).type
                    case 'plane'
                        nplanes = (length(curr_filt)-1)/3;
                        obj.filters(i).relate = cell(nplanes,1);
                        obj.filters(i).n = zeros(3,nplanes);
                        obj.filters(i).val = zeros(nplanes,1);
                        
                        for j = 1:nplanes
                            obj.filters(i).relate{j} = curr_filt{2+3*(j-1)};
                            obj.filters(i).n(:,j)    = curr_filt{3+3*(j-1)};
                            obj.filters(i).val(j)    = curr_filt{4+3*(j-1)};
                        end
                    case 'distmesh'
                        obj.filters(i).func = curr_filt{2};
                    case 'bbox'
                        nbbox = (length(curr_filt)-1)/2;
                        
                        for j = 1:nbbox
                            obj.filters(i).relate{j} = curr_filt{2+2*(j-1)};
                            obj.filters(i).val{j} = curr_filt{3+2*(j-1)};
                        end
                end
            end
            
            %Extended plot limits
            obj.xlim = [min(min(p(:,1))),max(max(p(:,1)))];
            obj.ylim = [min(min(p(:,2))),max(max(p(:,2)))];
            obj.zlim = [min(min(p(:,3))),max(max(p(:,3)))];
            dx=obj.xlim(2)-obj.xlim(1); dy=obj.ylim(2)-obj.ylim(1); dz=obj.zlim(2)-obj.zlim(1);
            obj.xlimX = obj.xlim + [-0.25,0.25]*dx;
            obj.ylimX = obj.ylim + [-0.25,0.25]*dy;
            obj.zlimX = obj.zlim + [-0.25,0.25]*dz;
            
            %Apply filters
            obj.applyFilters(p,t);
        end
        
        function  [] = setProperty(obj,prop,val,varargin)
            if iscell(obj.(prop))
                obj.(prop){varargin{:}}=val;
            else
                obj.(prop)=val;
            end
        end
        
        function  [] = generateAllPlots(obj,p,t)
            
            obj.hFig(1) = figure('pos',[316   677   560   420]); axes;
            obj.plotRegion(obj.hFig(1),p,t);
            
            obj.hFig(2) = figure('pos',[886   677   560   420]); axes;
            obj.plotRegion(obj.hFig(2),p,t);
            obj.plotPlanes(obj.hFig(2));
            obj.plotBBox(obj.hFig(2));
            
            obj.hFig(3) = figure('pos',[583   174   560   420]); axes;
            obj.plotRegion(obj.hFig(3),p,obj.elem);
        end
        
        function  [h]   = plotRegion(obj,hF,p,t,col)
            if nargin < 5 || isempty(col), col=[]; end;
            
            %Determine bounding box for view
            bbox(1,1) = min(min(p(:,1))); bbox(1,2) = max(max(p(:,1)));
            bbox(2,1) = min(min(p(:,2))); bbox(2,2) = max(max(p(:,2)));
            bbox(3,1) = min(min(p(:,3))); bbox(3,2) = max(max(p(:,3)));
            
            figure(hF);
            h=simpplot(p,t,[],col);
            
            hAxes=get(hF,'Children'); hAxes=hAxes(end);
            view(hAxes,3); grid(hAxes,'on');
            axis equal;
            set(hAxes,'xlim',bbox(1,:),'ylim',bbox(2,:),'zlim',bbox(3,:));
            hold on;
            
        end
        
        function  [surf2] = hardCopy(obj)
            surf2 = Region();
            props=properties(obj);
            for i=1:length(props)
                surf2.(props{i})=obj.(props{i});
            end
        end
        
        %Filters
        function  []  = applyFilters(obj,p,t)
            tF=t;
            for i = 1:obj.nfilters
                switch obj.filters(i).type
                    case 'plane'
                        tF = obj.filterWithPlane(p,tF,i);
                    case 'distmesh'
                        tF = obj.filterWithDMesh(p,tF,i);
                    case 'bbox'
                        tF = obj.filterWithBBox(p,tF,i);
                end
            end
            obj.elem=tF;
        end
        
        function  []  = plotPlanes(obj,hF)
            hAx = get(hF,'Children'); hAx=hAx(end);
            for i = 1:obj.nfilters
                if ~strcmpi(obj.filters(i).type,'plane')
                    continue;
                end
                
                nplanes = size(obj.filters(i).n,2);
                x=obj.xlimX; y=obj.ylimX; z=obj.zlimX;
                for j = 1:nplanes
                    if (obj.filters(i).n(1,j) ~= 0) && (obj.filters(i).n(2,j) == 0) && (obj.filters(i).n(3,j) == 0)
                        X = (obj.filters(i).val(j)/obj.filters(i).n(1,j))*ones(4,1);
                        Y = [y(1);y(2);y(2);y(1)];
                        Z = [z(1);z(1);z(2);z(2)];
                    elseif (obj.filters(i).n(1,j) == 0) && (obj.filters(i).n(2,j) ~= 0) && (obj.filters(i).n(3,j) == 0)
                        Y = (obj.filters(i).val(j)/obj.filters(i).n(2,j))*ones(4,1);
                        X = [x(1);x(2);x(2);x(1)];
                        Z = [z(1);z(1);z(2);z(2)];
                    elseif (obj.filters(i).n(1,j) == 0) && (obj.filters(i).n(2,j) == 0) && (obj.filters(i).n(3,j) ~= 0)
                        Z = (obj.filters(i).val(j)/obj.filters(i).n(3,j))*ones(4,1);
                        X = [x(1);x(2);x(2);x(1)];
                        Y = [y(1);y(1);y(2);y(2)];
                    elseif (obj.filters(i).n(3,j) == 0)
                        Y = [y(1);y(2);y(2);y(1)];
                        X = (1/obj.filters(i).n(1,j))*(obj.filters(i).val(j)-obj.filters(i).n(2,j)*Y);
                        Z = [z(1);z(1);z(2);z(2)];
                    else
                        X = [x(1);x(2);x(2);x(1)];
                        Y = [y(1);y(1);y(2);y(2)];
                        Z = (1/obj.filters(i).n(3,j))*(obj.filters(i).val(j) - obj.filters(i).n(1,j)*X - obj.filters(i).n(2,j)*Y);
                    end
                    switch obj.filters(i).relate{j}
                        case 'lt'
                            patch('parent',hAx,'faces',1:4,'vertices',[X,Y,Z],'facecolor','r');
                        case 'gt'
                            patch('parent',hAx,'faces',1:4,'vertices',[X,Y,Z],'facecolor','k');
                    end
                end
            end
            
        end
        
        function  []  = plotBBox(obj,hF)
            
            hAx = get(hF,'Children'); hAx=hAx(end);
            for i = 1:obj.nfilters
                if ~strcmpi(obj.filters(i).type,'bbox')
                    continue;
                end
                nbbox = size(obj.filters(i).val,2);
                for j = 1:nbbox
                    x=obj.filters(i).val{j}(1,:);
                    y=obj.filters(i).val{j}(2,:);
                    z=obj.filters(i).val{j}(3,:);
                    
                    f1x = x(1)*ones(4,1);
                    f1y = [y(1);y(2);y(2);y(1)];
                    f1z = [z(1);z(1);z(2);z(2)];
                    
                    f2x = x(2)*ones(4,1);
                    f2y = [y(1);y(2);y(2);y(1)];
                    f2z = [z(1);z(1);z(2);z(2)];
                    
                    f3x = [x(1);x(2);x(2);x(1)];
                    f3y =  y(1)*ones(4,1);
                    f3z = [z(1);z(1);z(2);z(2)];
                    
                    f4x = [x(1);x(2);x(2);x(1)];
                    f4y =  y(2)*ones(4,1);
                    f4z = [z(1);z(1);z(2);z(2)];
                    
                    f5x = [x(1);x(2);x(2);x(1)];
                    f5y = [y(1);y(1);y(2);y(2)];
                    f5z = z(1)*ones(4,1);
                    
                    f6x = [x(1);x(2);x(2);x(1)];
                    f6y = [y(1);y(1);y(2);y(2)];
                    f6z = z(2)*ones(4,1);
                    
                    switch obj.filters(i).relate{j}
                        case 'lt'
                            patch('parent',hAx,'faces',1:4,'vertices',[f1x,f1y,f1z],'facecolor','r','FaceAlpha','flat','FaceVertexAlphaData',0.4);
                            patch('parent',hAx,'faces',1:4,'vertices',[f2x,f2y,f2z],'facecolor','r','FaceAlpha','flat','FaceVertexAlphaData',0.4);
                            patch('parent',hAx,'faces',1:4,'vertices',[f3x,f3y,f3z],'facecolor','r','FaceAlpha','flat','FaceVertexAlphaData',0.4);
                            patch('parent',hAx,'faces',1:4,'vertices',[f4x,f4y,f4z],'facecolor','r','FaceAlpha','flat','FaceVertexAlphaData',0.4);
                            patch('parent',hAx,'faces',1:4,'vertices',[f5x,f5y,f5z],'facecolor','r','FaceAlpha','flat','FaceVertexAlphaData',0.4);
                            patch('parent',hAx,'faces',1:4,'vertices',[f6x,f6y,f6z],'facecolor','r','FaceAlpha','flat','FaceVertexAlphaData',0.4);
                        case 'gt'
                            patch('parent',hAx,'faces',1:4,'vertices',[f1x,f1y,f1z],'facecolor','k','FaceAlpha','flat','FaceVertexAlphaData',0.4);
                            patch('parent',hAx,'faces',1:4,'vertices',[f2x,f2y,f2z],'facecolor','k','FaceAlpha','flat','FaceVertexAlphaData',0.4);
                            patch('parent',hAx,'faces',1:4,'vertices',[f3x,f3y,f3z],'facecolor','k','FaceAlpha','flat','FaceVertexAlphaData',0.4);
                            patch('parent',hAx,'faces',1:4,'vertices',[f4x,f4y,f4z],'facecolor','k','FaceAlpha','flat','FaceVertexAlphaData',0.4);
                            patch('parent',hAx,'faces',1:4,'vertices',[f5x,f5y,f5z],'facecolor','k','FaceAlpha','flat','FaceVertexAlphaData',0.4);
                            patch('parent',hAx,'faces',1:4,'vertices',[f6x,f6y,f6z],'facecolor','k','FaceAlpha','flat','FaceVertexAlphaData',0.4);
                    end
                end
            end
        end
        
        function  [elF] = filterWithDMesh(obj,p,el,num)
            d=obj.filters(num).func(p);
            val=d(el);
            elF=el(any(val<0,2),:);
        end
        
        function  [elF] = filterWithPlane(obj,p,el,num)
            %Plane: n'*x + d = 0
            
            %Set function handle defining acceptable regions
            for j = 1:length(obj.filters(num).val)
                switch obj.filters(num).relate{j}
                    case 'lt'
                        func{j} = @(x) x*obj.filters(num).n(:,j) - obj.filters(num).val(j) < 0;
                    case 'gt'
                        func{j} = @(x) x*obj.filters(num).n(:,j) - obj.filters(num).val(j) > 0;
                end
            end
            
            is_elF = false(size(el,1),1);
            for i = 1:size(el)
                for j = 1:length(obj.filters(num).val)
                    is_elF(i) = is_elF(i) || all(func{j}(p(el(i,:),:)));
                end
            end
            elF=el(is_elF,:);
        end
        
        function  [elF] = filterWithBBox(obj,p,el,num)
            %Set function handle defining acceptable regions
            val=obj.filters(num).val;
            for j = 1:length(val)
                switch obj.filters(num).relate{j}
                    case 'lt'
                        func{j} = @(x) all(bsxfun(@gt,x',val{j}(:,1)) & bsxfun(@lt,x',val{j}(:,2)),1);
                    case 'gt'
                        func{j} = @(x) all(bsxfun(@lt,x',val{j}(:,1)) & bsxfun(@gt,x',val{j}(:,2)),1);
                end
            end
            
            is_elF = false(size(el,1),1);
            for i = 1:size(el)
                for j = 1:length(val)
                    is_elF(i) = is_elF(i) || any(func{j}(p(el(i,:),:)));
                end
            end
            elF=el(is_elF,:);
        end
                
        %zFEM communication
        function  [fExtNodal0] = nodalLoad(obj,dim,val,msh,dbc)
            
            fFull=zeros(msh.nsd,msh.np);
            fFull(dim,unique(obj.elem(:)))=val;
            fExtNodal0 = fFull(dbc.loc);
            fExtNodal0=fExtNodal0(:);
            
        end
        
        function   [b,bdbc] = bodyLoad(obj,dim,val,msh,dbc)
            
            bFull=zeros(msh.nsd,msh.np);
            bFull(dim,unique(obj.elem(:)))=val;
            b    = bFull(dbc.loc);
            bdbc = bFull(~dbc.loc);
            
            b=b(:); bdbc=bdbc(:);
            
        end
        
        function  [dbc] = dispBCs(obj,dim,msh)
            
            nodesInElem = unique(obj.elem(:));
            
            dbc.loc=false(msh.nsd,msh.np);
            dbc.loc(dim,nodesInElem)=true;
            dbc.ndbc=int64(sum(sum(dbc.loc)));

%             if any(dbc.loc(dim,nodesInElem))
%                 warning('Setting DBCs that were previously set in dbc structure.  Overwritting with new values in output structure.');
%             end
%             
%             uFull = zeros(msh.nsd,msh.np);
%             uFull(dbc.loc) = dbc.udbc;
%             
%             dbc_new = dbc;
%             dbc_new.loc(dim,nodesInElem)=true;
%             dbc_new.ndbc=sum(sum(dbc_new.loc));
%             uFull(dim,nodesInElem) = val;
%             dbc_new.udbc = uFull(dbc_new.loc);
%             dbc_new.udbc=dbc_new.udbc(:);            
        end
        
        function  [elemlist] = MatProp(obj,msh)
            
            [~,elemlist] = intersect(msh.IEN'+1,obj.elem,'rows');
            
        end
    end
    
    methods (Static)
        %Set Operations
        %need to think about how to modify filters of surf2 so they
        %represent the complement.  need to use de morgan's law.
        function  [reg2] = complement(t,reg1)
            
            reg2=Region();
            reg2.name = [reg1.name, ' - complement'];
            reg2.elem = setdiff(t,reg1.elem,'rows','stable');
            reg2.numelem = size(reg2.elem,1);
        end
        
        function  [reg3] = setdiff(reg1,reg2)
            
            reg3=Region();
            reg3.name = [reg1.name, ' - ', reg2.name];
            reg3.elem = setdiff(reg1.elem,reg2.elem,'rows','stable');
            reg3.numelem = size(reg3.elem,1);
        end
        
        function [reg1] = union(varargin)
            
            reg1 = Region();
            reg1.name = varargin{1}.name;
            reg1.elem = varargin{1}.elem;
            for i = 2:length(varargin)
                reg1.name = [reg1.name,' + ',varargin{i}.name];
                reg1.elem = union(reg1.elem,varargin{i}.elem,'rows');
            end
            reg1.numelem = size(reg1.elem,1);
        end
        
        function [reg1] = intersect(varargin)
            
            reg1 = Region();
            reg1.name = varargin{1}.name;
            reg1.elem = varargin{1}.elem;
            for i = 2:length(varargin)
                reg1.name = [reg1.name,' & ',varargin{i}.name];
                reg1.elem = intersect(reg1.elem,varargin{i}.elem,'rows');
            end
            reg1.numelem = size(reg1.elem,1);
        end
        
        function [mat] =  physicsFromRegions(msh,varargin)
            mat=zeros(msh.nel,1,'uint32');
            for i = 1:length(varargin)
                for j = 1:length(varargin{i})
                    mat(varargin{i}(j).MatProp(msh)) = uint32(i-1);
                end
            end
        end
        
        function [mat,phys] = intersectMatRegions(mat1,phys1,mat2,phys2)
            %This function assumes that regions of lambda values and mu
            %values were created independently, with the possibilty that
            %the regions used to define them were not the same.  This
            %function defines new materials so every combination can be
            %represented with physics objects.
            
            nel=length(mat1);
            mat=zeros(nel,1,'uint32');
            mismatch = find(mat1~=mat2);
            match = setdiff(1:nel,mismatch);
            
            mat(match)=mat1(match);
            
            matchedmats = unique(mat1(match));
            for i = matchedmats(:)'
                phys(i+1).lam = phys1(i+1).lam;
                phys(i+1).mu  = phys2(i+1).mu;
            end
            
            nummat = length(matchedmats);
            while ~isempty(mismatch)
                newmatch = find((mat1 == mat1(mismatch(1))) & (mat2 == mat2(mismatch(1))));
                mat(newmatch)=nummat;
                phys(nummat+1).lam=phys1(mat1(mismatch(1))+1).lam;
                phys(nummat+1).mu =phys2(mat2(mismatch(1))+1).mu;
                
                nummat=nummat+1;
                mismatch = setdiff(mismatch,newmatch);
            end
        end
        
        function [geom] = intersectRegionsFromGeom(varargin)
            
            if length(varargin) == 1 || isempty(varargin{2})
                geom=varargin{1};
                return;
            end
            
            geom1=varargin{1};
            m = length(geom1);
            
            cnt = 0;
            for i = 1:m
                for j = 1:length(varargin{2})
                    intersect_regions = Region.intersect(geom1(i),varargin{2}(j));
                    if ~isempty(intersect_regions.elem)
                        cnt=cnt+1;
                        geom(cnt)=intersect_regions;
                    end
                end
            end
            
            if length(varargin) == 2
                return;
            end
            
            geom=Region.intersectRegionsFromGeom(geom,varargin{3:end});
%                 for j = 2:length(varargin)
%                     region1 = varargin{j};
%                     
%                     for k = 1:length(region1)
%                         intersect_regions = Region.intersect(geom1(i),region1 (k));
%                         if ~isempty(intersect_regions.elem)
%                             cnt=cnt+1;
%                             geom(cnt)=intersect_regions;
%                         end
%                     end
%                 end
%             end
        end
        
        function [mat,phys] = physFromGeomMatSnaps(geom,lam,mu,t)
            %Geometry (array of regions whose union is entire domain)
            %lam = nel array defining material distribution of lambda
            %mu  = nel array defining material distribution of mu
            
            mat=zeros(length(lam),1,'uint32');
            for i = 1:length(geom)
                [~,elnum] = intersect(t,geom(i).elem,'rows','stable');
                phys(i).lam = lam(elnum(1));
                phys(i).mu  =  mu(elnum(1));
                mat(elnum) = i-1;
            end
        end
        
        function [] = plotMatDist(geom,X,p,t,bnds,varargin)
            %additional arguments are regions to intersect with
            
            for i = 1:length(geom)
                geomNew(i) = geom(i).hardCopy();
                geomNew(i) = Region.intersect(geomNew(i),varargin{:});
            end
            
            %bnds=[min(X),max(X)];
            
            map = colormap('jet');
            %map = colormap('bone');
            ncol = size(map,1);
            
            [~,bin] = histc(X,linspace(bnds(1),bnds(2),ncol));
            
            for i = 1:length(geom)
                [~,elnum] = intersect(double(t),geom(i).elem,'rows','stable');
                if bin(elnum(1)) == 0
                    continue;
                end
                geomNew(i).plotRegion(gcf,p,geomNew(i).elem,map(bin(elnum(1)),:));
            end
            if bnds(2)>bnds(1)
                set(gca,'CLim',bnds);
                colorbar;
            end
        end
        
        function [Xlam,Xmu,Xrho]  = matphys2vec(mat,phys)
            Xlam=[phys(mat+1).lam];  Xlam=Xlam(:);
            Xmu =[phys(mat+1).mu];   Xmu =Xmu(:);
            Xrho=[phys(mat+1).rho0]; Xrho=Xrho(:);
        end
        
        function [mat,phys] = vec2matphys(Xlam,Xmu,Xrho)
            
            [matprop,~,mat] = unique([Xlam,Xmu,Xrho],'rows');
            mat=uint32(mat(:)-1);
            
            nmat=size(matprop,1);
            phys = repmat(struct('lam',[],'mu',[],'rho0',[]),1,nmat);
            for i = 1:nmat
                phys(i).lam = matprop(i,1);
                phys(i).mu  = matprop(i,2);
                phys(i).rho0= matprop(i,3);
            end
            
            return;
            
            nel=size(Xlam,1);
            mat=uint32(0:nel-1);mat=mat(:);
            
            tmp=num2cell(Xlam); clear Xlam;
            [phys(1:nel).lam]=deal(tmp{:});
            tmp=num2cell(Xmu); clear Xmu;
            [phys(1:nel).mu] =deal(tmp{:});
            tmp=num2cell(Xrho); clear Xrho;
            [phys(1:nel).rho0] =deal(tmp{:});
        end
    end
end
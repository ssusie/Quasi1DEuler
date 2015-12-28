classdef Surface < handle
    
    properties (SetAccess=private, GetAccess=public)
        name;
        
        tri;
        numtri;
        
        filters;
        nfilters;
    end
    
    properties (Hidden=true, SetAccess=private, GetAccess=public)
        bcol=[.8,.9,1];
        hcol=[.5,.6,1];
        
        hFig;
        hTri;
        
        fTri2aTri;

        xlim;
        ylim;
        zlim;
        xlimX;
        ylimX;
        zlimX;
    end
    
    methods
        function [obj] = Surface(surfname,p,varargin)
            
            obj.name = surfname;

            if nargin==1, return; end;
            
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
                    case 'normal'
                        nnormals = length(curr_filt)-1;
                        obj.filters(i).relate = [];
                        obj.filters(i).n = zeros(3,nnormals);
                        obj.filters(i).val=[];
                        
                        for j = 1:nnormals
                            obj.filters(i).n(:,j) = curr_filt{j+1};
                            obj.filters(i).n(:,j) = obj.filters(i).n(:,j)/norm(obj.filters(i).n(:,j));
                        end
                    case 'bbox'
                        nbbox = (length(curr_filt)-1)/2;
                        
                        for j = 1:nbbox
                            obj.filters(i).relate{j} = curr_filt{2+2*(j-1)};
                            obj.filters(i).val{j} = curr_filt{3+2*(j-1)};
                        end
                    case 'distmesh'
                        obj.filters(i).func = curr_filt{2};
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
        end
        
        function  [] = setProperty(obj,prop,val,varargin)
            if iscell(obj.(prop))
                obj.(prop){varargin{:}}=val;
            else
                obj.(prop)=val;
            end
        end
        
        function [act] = addremoveTriangle(obj,newtri)
            
            loc = ismember(obj.tri,newtri,'rows');
            if sum(loc) > 0
                %If the current element is already highlighted, 
                %unhighlight it and remove from triangle list
                act='remove';
                obj.numtri = obj.numtri-1;
                obj.tri(loc,:)=[];
            else
                act='add';
                obj.numtri = obj.numtri+1;
                obj.tri=[obj.tri;newtri];
            end
        end
        
        function  [] = generateAllPlots(obj,p,t)
            
            tri1=surftri(p,t);
            
            obj.hFig(1) = figure('pos',[316   677   560   420],...
                'DeleteFcn',{@(Esrc,Eobj,obj) obj.setProperty('hTri',[],1),obj});
            axes;
            obj.hTri{1}=obj.plotSurface(obj.hFig(1),p,t,false,tri1);
            
            obj.hFig(2) = figure('pos',[886   677   560   420],...
                'DeleteFcn',{@(Esrc,Eobj,obj) obj.setProperty('hTri',[],2),obj});
            axes;
            obj.hTri{2}=obj.plotSurface(obj.hFig(2),p,t,false,tri1);
            obj.plotPlanes(obj.hFig(2));
            
            obj.hFig(3) = figure('pos',[583   174   560   420],...
                'DeleteFcn',{@(Esrc,Eobj,obj) obj.setProperty('hTri',[],3),obj});
            axes;
            obj.hTri{3}=obj.plotSurface(obj.hFig(3),p,t,true,tri1);
            
            obj.associateFacesBetweenFigs();
        end
               
        function  [hT]   = plotSurface(obj,hF,p,t,filterflag,tri1)
            %Determine bounding box for view
            bbox(1,1) = min(min(p(:,1))); bbox(1,2) = max(max(p(:,1)));
            bbox(2,1) = min(min(p(:,2))); bbox(2,2) = max(max(p(:,2)));
            bbox(3,1) = min(min(p(:,3))); bbox(3,2) = max(max(p(:,3)));

            %Extract surface triangles
            if nargin < 6 || isempty(tri1), tri1=surftri(p,t); end;

            %Apply filters
            triF = obj.applyFilters(p,tri1);
            if filterflag, tri2=triF; else tri2=tri1; end;
            
            %Plot triangles
            hAxes = get(hF,'Children');
            for i = 1:size(tri2,1)
                tri3 = tri2(i,:);

                hT(i)=patch('faces',1:3,'vertices',p(tri3,:),...
                    'facevertexcdata',p(tri2,3),'facecolor',obj.bcol,...
                    'edgecolor','k','parent',hAxes,...
                    'ButtonDownFcn',{@obj.highlightTriangle,tri3});
            end
            view(hAxes,3); grid(hAxes,'on');
            axis equal;
            set(hAxes,'xlim',bbox(1,:),'ylim',bbox(2,:),'zlim',bbox(3,:));
            hold on;
        end

        function  [] = highlightTriangle(obj,Esrc,Edata,tri2)
            %Add or remove triangle to surface and highlight element appropriately
            act=obj.addremoveTriangle(tri2);
            
            switch act
                case 'add'
                    curr_color = obj.hcol;
                case 'remove'
                    curr_color = obj.bcol;
            end
            
            otherTriHan = getappdata(Esrc,'otherTriH');
            if size(otherTriHan,2) == 1 && strcmpi(act,'add')
                disp('Trying to add face that did not make it through filter to surface.  This is not allowed.  Skipping face.');
                return;
            end
            set([Esrc,otherTriHan(ishandle(otherTriHan))],'facecolor',curr_color);
        end
        
        function  [] = associateFacesBetweenFigs(obj)
            
            nfaces = length(obj.hTri{1});
            
            for j = 1:nfaces
                
                tmp=find(obj.fTri2aTri==j);
                filthan=[];
                if ~isempty(tmp), filthan=obj.hTri{3}(tmp); end;
                
                hanvec1 = [obj.hTri{2}(j),filthan];
                hanvec2 = [obj.hTri{1}(j),filthan];
                hanvec3 = [obj.hTri{1}(j),obj.hTri{2}(j)];
                
                setappdata(obj.hTri{1}(j),'otherTriH',hanvec1);
                setappdata(obj.hTri{2}(j),'otherTriH',hanvec2);
                if ~isempty(filthan), setappdata(filthan,'otherTriH',hanvec3); end;
            end
        end
                
        function  [] = selectAllTri(obj,p,t)
                        
            tri1=surftri(p,t);
            triF = obj.applyFilters(p,tri1);
            
            obj.tri = triF;
            obj.numtri = size(obj.tri,1);

            for i = 1:length(obj.hFig)
                if ~ishandle(obj.hFig(i)), continue; end;
                if i==3, allhTri=obj.hTri{i}; else allhTri = obj.hTri{i}(obj.fTri2aTri); end;
                set(allhTri,'facecolor',obj.hcol);
            end
        end
        
        function [] = unselectAllTri(obj)
            
            obj.tri=[];
            obj.numtri=0;
            
            for i = 1:length(obj.hFig)
                if ~ishandle(obj.hFig(i)), continue; end;
                set(obj.hTri{i},'facecolor',obj.bcol);
            end
        end
        
        function [] = surfUnion(obj,surf2)
            obj.name = [obj.name,' + ',surf2.name];
            obj.tri  = [obj.tri;surf2.tri];
            obj.numtri = obj.numtri + surf2.numtri;
        end
        
        %Filters
        function  [triF]  = applyFilters(obj,p,tri1)
            triF=tri1;
            for i = 1:obj.nfilters
                switch obj.filters(i).type
                    case 'plane'
                        triF = obj.filterWithPlane(p,triF,i);
                    case 'normal'
                        triF = obj.filterWithNormal(p,triF,i);
                    case 'bbox'
                        triF = obj.filterWithBBox(p,triF,i);
                    case 'distmesh'
                        triF = obj.filterWithDMesh(p,triF,i);
                end
            end
            [~,obj.fTri2aTri]=intersect(tri1,triF,'rows');
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
        
        function  [triF] = filterWithPlane(obj,p,tri1,num)
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
            
            is_triF = false(size(tri1,1),1);
            for i = 1:size(tri1)
                for j = 1:length(obj.filters(num).val)
                    is_triF(i) = is_triF(i) || all(func{j}(p(tri1(i,:),:)));
                end
            end
            triF=tri1(is_triF,:);
        end
        
        function [triF] = filterWithNormal(obj,p,tri1,num)
            
            %Set function handle defining acceptable regions
            for i = 1:size(obj.filters(num).n,2)
                func{i} = @(x) abs(1-obj.filters(num).n(:,i)'*(cross(x(2,:)'-x(1,:)',x(3,:)'-x(1,:)')/norm(cross(x(2,:)'-x(1,:)',x(3,:)'-x(1,:)')))) < 1e-1;
            end
            
            is_triF = false(size(tri1,1),1);
            for i = 1:size(tri1)
                for j = 1:size(obj.filters(num).n,2)
                    is_triF(i) = is_triF(i) || any(func{j}(p(tri1(i,:),:)));
                end
            end
            triF=tri1(is_triF,:);
        end
        
        function [triF] = filterWithBBox(obj,p,tri1,num)
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
            
            is_triF = false(size(tri1,1),1);
            for i = 1:size(tri1)
                for j = 1:length(val)
                    is_triF(i) = is_triF(i) || any(func{j}(p(tri1(i,:),:)));
                end
            end
            triF=tri1(is_triF,:);
        end
        
        function  [triF] = filterWithDMesh(obj,p,tri1,num)
            d=obj.filters(num).func(p);
            val=d(tri1);
            triF=tri1(any(val<0,2),:);
        end
        
        %zFEM communication
        function  [fExtNodal0] = distributedLoad(obj,dim,val,msh,dbc)
                        
            fFull=zeros(msh.nsd,msh.np);
            for i = 1:obj.numtri
                v1 = msh.X(:,obj.tri(i,2)) - msh.X(:,obj.tri(i,1));
                v2 = msh.X(:,obj.tri(i,3)) - msh.X(:,obj.tri(i,1));
                th = acos((v1'*v2)/(norm(v1)*norm(v2)));
                A = norm(v1)*norm(v2)*sin(th);
                
                fFull(dim,obj.tri(i,:)) = A*val/3;
            end
            fExtNodal0 = fFull(~dbc.loc);            
        end
        
        function  [fExtNodal0] = distributedLoad_func(obj,dim,b,msh,dbc)
            [W,Z] = gaussquad_tri(3,2,1);
            
            fFull = zeros(msh.nsd,msh.np);
            for i = 1:obj.numtri
                v1 = msh.X(:,obj.tri(i,1)); v1=v1(:);
                v2 = msh.X(:,obj.tri(i,2)); v2=v2(:);
                v3 = msh.X(:,obj.tri(i,3)); v3=v3(:);
                
                drdzeta_drdeta = norm(cross(v2-v1,v3-v1));
                f = 0;
                for j = 1:length(W)
                    f = f + 0.5*W(j)*b(Z(:,j))*drdzeta_drdeta;
                end
                fFull(dim,obj.tri(i,:)) = fFull(dim,obj.tri(i,:)) + f/3;
            end
            fExtNodal0 = fFull(~dbc.loc);
        end
        
        function  [fExtNodal0] = nodalLoad(obj,dim,val,msh,dbc)
            
            fFull=zeros(msh.nsd,msh.np);
            fFull(dim,unique(obj.tri(:)))=val;
            fExtNodal0 = fFull(dbc.loc);
            fExtNodal0=fExtNodal0(:);
            
        end
        
        function  [dbc] = dispBCs(obj,dim,msh)
            
            nodesOnSurf = unique(obj.tri(:));
            
            %uFull = zeros(msh.nsd,msh.np);
            %uFull(dbc.loc) = dbc.udbc(1); %udbc at time 1 (assuming static)
            
            %dbc_new = dbc;
            
            dbc.loc=false(msh.nsd,msh.np);
            dbc.loc(dim,nodesOnSurf)=true;
            dbc.ndbc=int64(sum(sum(dbc.loc)));
            
            %uFull(dim,nodesOnSurf) = val;
            %dbc_new.udbc = uFull(dbc_new.loc);
            %dbc_new.udbc=dbc_new.udbc(:);
            
        end
    end
    
    methods (Static)
        function  [surf3] = setdiff(surf1,surf2)
            Name = surf1.name;
            for i = 1:length(surf2)
                Name = [Name, ' - ',surf2(i).name];
            end
            surf3 = Surface(Name);
            
            surf3.tri = setdiff(surf1.tri,surf2(1).tri,'rows','stable');
            for i = 2:length(surf2)
                surf3.tri = setdiff(surf3.tri,surf2(i).tri,'rows','stable');
            end
            surf3.numtri = size(surf3.tri,1);
        end
        
        function  [surf3] = intersect(surf1,surf2)
            Name = surf1.name;
            for i = 1:length(surf2)
                Name = [Name, ' & ',surf2(i).name];
            end
            surf3 = Surface(Name);
            
            surf3.tri = intersect(surf1.tri,surf2(1).tri,'rows','stable');
            for i = 2:length(surf2)
                surf3.tri = intersect(surf3.tri,surf2(i).tri,'rows','stable');
            end
            surf3.numtri = size(surf3.tri,1);
        end
    end
end
classdef ppAXES < handle
   
    properties (SetAccess = private, GetAccess = public)
        %Input Properties
        han = [];
        objhan = [];
        leghan = [];
        
        figId = [];
        id   = [];
        setax = {};
        subplotNum = 1;
        
        xData = [];
        xKey = [];
        
        yData = [];
        yKey = [];
        
        zData = [];
        zKey = [];
        
        scaleData = [1,1,1];
        
        model = [];
        modelAux = [];
        modelFom = [];
        pType = {'plot'};
        numPlotObj = [];
        
        legend = [];
        
        connWithCurve = [];
        connCurveSpec = [];
        
        normType = 2;
        
        xLabel = [];
        yLabel = [];
        zLabel = [];
        Title = [];
        
        plotspec = {'k-','linewidth',2};
    end
    
    methods
        function  [obj] = ppAXES(fname,indR,indC,spShape,ppobj)
            %This is the constructor of the PostProcess class.  It reads
            %the appropriate entries from the input file, stores them in
            %the class instance, and generates the requested plots.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %fname - string containing the postprocessing file name
            %indR  - scalar indicating the index of the row of the
            %        desired AXES object that will come out of readInFile
            %indC  - scalar indicating the index of the column of the
            %        desired AXES object that will come out of readInFile
            %spShape - 1x2 array containing the subplot configuration of
            %          the current figure
            %ppobj   - PostProcess object
            %
            %Outputs:
            %--------
            %obj   - instance of the ppAXES object that was
            %        constructed
            %--------------------------------------------------------------
            
            %Extract parameters from PP file
            VARtext = readInFile('VAR',fname,1);
            
            [~,AXcell] = readInFile('AXES',fname,[]);
            obj.determineAX(AXcell{indR,indC},[ppobj.GVARtxt;VARtext],ppobj);
            %Set axes and establish the nextplot behavior
            if ~strcmpi(obj.pType{1},'animate')
                obj.han = subplot(spShape(1),spShape(2),obj.subplotNum);
                set(obj.han,'nextplot','add');
            end
            obj.generatePlots;
            obj.plotConnectionCurves;
            %Set the user specified figure properties as long as the cell
            %array is not empty
            if ~isempty(obj.setax)
                set(obj.han,obj.setax{:});
            end
        end
        
        function  [] = determineAX(obj,AXtext,VARtxt,ppobj)
            %This function converts the text from the AXES object into a
            %structure whose fields correspond to the ppAXES properties
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj    - ppAXES object
            %AXcell - a string containing the text from the AXES objects in
            %         the input file
            %ppobj  - PostProcess object
            %
            %Outputs:
            %--------
            %There are no outputs.  Properties are stored in the FIGURE
            %handle class.
            %--------------------------------------------------------------
            
            %Extract variables from VARtxt and evaluate them in the current
            %workspace.
            if ~isempty(VARtxt)                
                for i = 1:length(VARtxt)
                    eval([VARtxt{i},';']);
                end
            end
            
            %Determine id of current AXES object
            obj.id = checkPostInput(AXtext,1,'id',obj.id);
            
            %Determine subplotNum of current AXES object
            obj.subplotNum = checkPostInput(AXtext,1,'subplotNum',obj.subplotNum);
            
            %Determine table properies
            obj.setax = checkPostInput(AXtext,1,'setax',obj.setax); 
            
            %Determine the norm type to use in this AXES object
            obj.normType = checkPostInput(AXtext,1,'normType',obj.normType);
            
            %Determine number of objects to include in the axes
            obj.numPlotObj = checkPostInput(AXtext,1,'numPlotObj',obj.numPlotObj);
            N = obj.numPlotObj;
            
            %Determine data scaling
            obj.scaleData = checkPostInput(AXtext,N,'scaleData',obj.scaleData);
            if size(obj.scaleData,1) ~= N
                obj.scaleData = repmat(obj.scaleData,N,1);
            end
            
            %Determine number of connection curves in the axes
            obj.connWithCurve = checkPostInput(AXtext,[],'connWithCurve',obj.connWithCurve);
            M = length(obj.connWithCurve);
            
            %Determine the connnection curve properties
            connCurveSpecVec = checkPostInput(AXtext,M,'connCurveSpec',obj.connCurveSpec);
            obj.connCurveSpec= PostInputPlotFormat(connCurveSpecVec,M,'connCurveSpec');
            
            %Determine legend parameters
            obj.legend = checkPostInput(AXtext,[],'legend',obj.legend);
            
            %Determine model to extract information from.
            obj.model = checkPostInput(AXtext,N,'model',obj.model);
            
            %Determine auxilary model to extract information from.
            obj.modelAux = checkPostInput(AXtext,N,'modelAux',obj.modelAux);
            
            %Determine FOM model to extract information from.
            obj.modelFom = checkPostInput(AXtext,N,'modelFom',obj.modelFom);
            
            %Determine the type of plotting function to use
            pTypevec = checkPostInput(AXtext,N,'pType',obj.pType);
            if size(pTypevec,1) == 1
                obj.pType = cell(N,1);
                [obj.pType{:}] = deal(pTypevec{:});
            elseif size(pTypevec,1) == N
                obj.pType = pTypevec;
            else
                error('In the AXES field, pType must be either have 1 entry or numPlotObj entries');
            end
            
            %Determine the x, y, and z labels and the title
            obj.xLabel = checkPostInput(AXtext,1,'xlabel',obj.xLabel);
            obj.yLabel = checkPostInput(AXtext,1,'ylabel',obj.yLabel);
            obj.zLabel = checkPostInput(AXtext,1,'zlabel',obj.zLabel);
            obj.Title = checkPostInput(AXtext,1,'title',obj.Title);
            
            %Determine the plotspec
            plotspecvec = checkPostInput(AXtext,N,'plotspec',obj.plotspec);
            obj.plotspec= PostInputPlotFormat(plotspecvec,N,'plotspec');
            
            %Determine the keywords for plotting
            xkeyvec = checkPostInput(AXtext,N,'xData',obj.xData);
            if size(xkeyvec,1) == 1
                obj.xKey = cell(N,1);
                [obj.xKey{:}] = deal(xkeyvec{:});
            elseif size(xkeyvec,1) == N
                obj.xKey = xkeyvec;
            elseif isempty(xkeyvec)
                obj.xKey = cell(N,1);
            else
                error('In the AXES field, xData must be either have 1 entry or numPlotObj entries');
            end
            
            ykeyvec = checkPostInput(AXtext,N,'yData',obj.yData);
            if size(ykeyvec,1) == 1
                obj.yKey = cell(N,1);
                [obj.yKey{:}] = deal(ykeyvec{:});
            elseif size(ykeyvec,1) == N
                obj.yKey = ykeyvec;
            elseif isempty(ykeyvec)
                obj.yKey = cell(N,1);
            else
                error('In the AXES field,yData must be either have 1 entry or numPlotObj entries');
            end
            
            zkeyvec = checkPostInput(AXtext,N,'zData',obj.zData);
            if size(zkeyvec,1) == 1 
                obj.zKey = cell(N,1);
                [obj.zKey{:}] = deal(zkeyvec{:});
            elseif size(zkeyvec,1) == N
                obj.zKey = zkeyvec;
            elseif isempty(zkeyvec)
                obj.zKey = cell(N,1);
            else
                error('In the AXES field, zData must be either have 1 entry or numPlotObj entries');
            end
            
            %Generate data for plotting
            obj.xData = cell(N,1);
            obj.yData = cell(N,1);
            obj.zData = cell(N,1);
            
            for k = 1:obj.numPlotObj
                modelobj = determineModelFromText(ppobj,obj.model,k,1,false);
                modelAuxobj = determineModelFromText(ppobj,obj.modelAux,k,1,true);
                modelFomobj = determineModelFromText(ppobj,obj.modelFom,k,1,false);
                
                obj.xData{k,1} = extractDataFromKey(obj,obj.xKey{k,1},modelobj,modelAuxobj,modelFomobj);
                obj.xData{k,1} = obj.scaleData(k,1)*obj.xData{k,1};
                
                obj.yData{k,1} = extractDataFromKey(obj,obj.yKey{k,1},modelobj,modelAuxobj,modelFomobj);
                obj.yData{k,1} = obj.scaleData(k,2)*obj.yData{k,1};
                
                obj.zData{k,1} = extractDataFromKey(obj,obj.zKey{k,1},modelobj,modelAuxobj,modelFomobj);
                obj.zData{k,1} = obj.scaleData(k,3)*obj.zData{k,1};
            end
        end
        
        function  [] = generatePlots(obj)
            %This functions generates all of the plots specified in the
            %PostProcess input file fname.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj   - ppAXES object
            %
            %Outputs:
            %--------
            %This function generates no outputs.  It generates the
            %requested plots.
            %--------------------------------------------------------------
            
            if ~strcmpi(class(obj.xData{1,1}),'function_handle')
                %If the user has requested an animation plot, the xData
                %will be a function handle to the problem animation plot
                %function.  In this case, the bounds of the plot will be
                %set internally.
                xmin = min(min(obj.xData{1,1})); xmax = max(max(obj.xData{1,1}));
                ymin = min(min(obj.yData{1,1})); ymax = max(max(obj.yData{1,1}));
                zmin = min(min(obj.zData{1,1})); zmax = max(max(obj.zData{1,1}));
                
                for k = 2:obj.numPlotObj
                    xmin = min(xmin,min(min(obj.xData{k,1})));
                    xmax = max(xmax,max(max(obj.xData{k,1})));
                    ymin = min(ymin,min(min(obj.yData{k,1})));
                    ymax = max(ymax,max(max(obj.yData{k,1})));
                    zmin = min(zmin,min(min(obj.zData{k,1})));
                    zmax = max(zmax,max(max(obj.zData{k,1})));
                end
            else
                %This is the case where the animation was requested in this
                %axes.  We assume that the user requested ONLY this plot in
                %the axes.
                f = obj.xData{1,1};
                f(obj.plotspec(1,:),2);
                return;
            end

            for k = 1:obj.numPlotObj
                x = obj.xData{k,1};
                y = obj.yData{k,1};
                z = obj.zData{k,1};
                
                if length(y) == 1 && length(x) > 1
                    y = repmat(y,length(x),1);
                elseif length(x) == 1 && length(y) > 1
                    x = repmat(x,length(y),1);
                end
                
                switch obj.pType{k} 
                    case 'plot'
                        obj.objhan(k) = plot(obj.han,x,y,obj.plotspec{k}{:});
                        if xmax > xmin && ymax > ymin
                            set(obj.han,'xlim',[xmin,xmax],'ylim',[ymin,ymax]);
                        end
                    case 'semilogx'
                        obj.objhan(k) = plot(obj.han,x,y,obj.plotspec{k}{:});
                        if xmax > xmin && ymax > ymin
                            set(obj.han,'xscale','log','xlim',[xmin,xmax],'ylim',[ymin,ymax]);
                        end
                    case 'semilogy'
                        obj.objhan(k) = plot(obj.han,x,y,obj.plotspec{k}{:});
                        if xmax > xmin && ymax > ymin
                            set(obj.han,'yscale','log','xlim',[xmin,xmax],'ylim',[ymin,ymax]);
                        end
                    case 'loglog'
                        obj.objhan(k) = plot(obj.han,x,y,obj.plotspec{k}{:});
                        if xmax > xmin && ymax > ymin
                            set(obj.han,'xscale','log','yscale','log','xlim',[xmin,xmax],'ylim',[ymin,ymax]);
                        end
                    case 'plot3'
                        obj.objhan(k) = plot3(obj.han,x,y,z,obj.plotspec{k}{:},'xlim',[xmin,xmax],'ylim',[ymin,ymax],'zlim',[zmin,zmax]);
                    case 'surf'
                        obj.objhan(k) = surf(obj.han,x,y,z,obj.plotspec{k}{:});
                        if xmax > xmin && ymax > ymin && zmax > zmin
                            set(obj.han,'xlim',[xmin,xmax],'ylim',[ymin,ymax],'zlim',[zmin,zmax]);
                        end
                    case 'surfc'
                        obj.objhan(k) = surfc(obj.han,x,y,z,obj.plotspec{k}{:});
                        if xmax > xmin && ymax > ymin && zmax > zmin
                            set(obj.han,'xlim',[xmin,xmax],'ylim',[ymin,ymax],'zlim',[zmin,zmax]);
                        end
                    case 'mesh'
                        obj.objhan(k) = mesh(obj.han,x,y,z,obj.plotspec{k}{:});
                    case 'contour'
                        obj.objhan(k) = contour(obj.han,x,y,z,obj.plotspec{k}{:});
                        if xmax > xmin && ymax > ymin && zmax > zmin
                            set(obj.han,'xlim',[xmin,xmax],'ylim',[ymin,ymax],'zlim',[zmin,zmax]);
                        end
                    case 'contourf'
                        obj.objhan(k) = contourf(obj.han,x,y,z,obj.plotspec{k}{:});
                        if xmax > xmin && ymax > ymin && zmax > zmin
                            set(obj.han,'xlim',[xmin,xmax],'ylim',[ymin,ymax],'zlim',[zmin,zmax]);
                        end
                end
            end
            if ~isempty(obj.xLabel)
                xlabel(obj.xLabel{:});
            end
            if ~isempty(obj.yLabel)
                ylabel(obj.yLabel{:});
            end
            if ~isempty(obj.zLabel)
                zlabel(obj.zLabel{:});
            end
            if ~isempty(obj.Title)
                title(obj.Title{:});
            end
            if ~isempty(obj.legend)
                legprops = obj.legend{1};
                for j = 1:length(obj.legend)-1
                    %Set up vector of object handles to be included in the
                    %legend
                    temp = obj.legend{j+1};
                    temphan(j) = obj.objhan(temp{1});
                    str(j) = temp(2);
                end
                obj.leghan = legend(temphan,str{:});
                set(obj.leghan,legprops{:});
            end
        end
        
        function  [] = plotConnectionCurves(obj)
            %This function plots the curves that connect the points
            %indicated in the input file.  This is only meant to be used
            %with the plot is 2D or plot3 (i.e. not mesh or surf, etc.
            %because it doesn't make sense in that context).
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj - ppAXES object
            %
            %Outputs:
            %--------
            %There are no outputs.
            %--------------------------------------------------------------
            
            if isempty(obj.connWithCurve)
                %There is nothing to connect so just return;
                return;
            end
            
            for j = 1:length(obj.connWithCurve)
                connectInfo = obj.connWithCurve{j};
                M = length(connectInfo);
                x = []; y = []; z = [];
                for k = 1:M
                    plotInd = connectInfo(k);
                    x = [x; obj.xData{plotInd}];
                    y = [y; obj.yData{plotInd}];
                    z = [z; obj.zData{plotInd}];
                end
                
                plotFormat = obj.connCurveSpec{j};
                switch obj.pType{plotInd}
                    case {'plot','loglog','semilogx','semilogy'}
                        plot(x,y,plotFormat{:});
                    case 'plot3'
                        plot(x,y,z,plotFormat{:});
                    otherwise
                        error('If your plot is not 2D or plot3 (i.e. mesh, surf, contour, etc), connWithCurve must be empty: [] or {}');
                end
            end
        end
    end
end
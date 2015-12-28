classdef ppTABLE < handle
    
    properties (SetAccess = private, GetAccess = public)
        %Input Properties
        han = [];
        panelHan = [];
        
        id   = [];
        settab = {};
        subtableNum = 1;
        normType = 2;
        
        numRows;
        numCols;
        numElem;
        
        rowKey;
        colKey;
        elemKey;
        rowData;
        colData;
        elemData;
        
        model;
        modelAux;
        modelFom;
        
        title;
        rowNames;
        colNames;
    end
    
    methods
        function  [obj] = ppTABLE(fname,indR,indC,spShape,ppobj)
            %This is the constructor of the ppTABLE class.  It reads
            %the appropriate entries from the input file, stores them in
            %the class instance, and generates the requested table.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %fname - string containing the postprocessing file name
            %indR  - scalar indicating the index of the row of the
            %        desired TABLE object that will come out of readInFile
            %indC  - scalar indicating the index of the column of the
            %        desired TABLE object that will come out of readInFile
            %spShape - 1x2 array containing the subtable configuration of
            %          the current figure
            %ppobj - PostProcess object
            %
            %Outputs:
            %--------
            %obj   - instance of the ppTABLE object that was
            %        constructed 
            %--------------------------------------------------------------
            
            panelGrid = spShape;
            
            %Extract parameters from PP file
            VARtext = readInFile('VAR',fname,1);
            
            [~,TABcell] = readInFile('TABLE',fname,[]);
            obj.panelHan = uipanel;
            data = obj.determineTAB(TABcell{indR,indC},[ppobj.GVARtxt;VARtext],ppobj);
            
            boundTol    = 0.02;
            lengthPanel = (1 - boundTol*(panelGrid(2)+1))/panelGrid(2);
            heightPanel = (1 - boundTol*(panelGrid(1)+1))/panelGrid(1);
            indMat = reshape(1:prod(panelGrid),panelGrid(2),panelGrid(1))';
            [numV,numH] = find(indMat == obj.subtableNum);
            LeftEdge = lengthPanel*(numH-1) + boundTol*(numH);
            BottomEdge = heightPanel*(panelGrid(1) - numV) + boundTol*(numV);
            set(obj.panelHan,'position',[LeftEdge,BottomEdge,lengthPanel,heightPanel]);
            
            
            obj.han = uitable(obj.panelHan,'Data',data);
            set(obj.han,'ColumnName',obj.colNames,'RowName',obj.rowNames,'Units','normalized');
            newPos = [0.05, 0.05, 0.9, 0.9];
            
            set(obj.han,'position',newPos);
            %Set the user specified figure properties as long as the cell
            %array is not empty
            if ~isempty(obj.settab)
                set(obj.han,obj.settab{:});
            end
        end
            
        function  [data] = determineTAB(obj,TABtext,VARtxt,ppobj)
            %This function converts the text from the TABLE object into a
            %structure whose fields correspond to the ppTABLE properties
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - ppTABLE object
            %TABtext - a string containing the text from the ppTABLE objects
            %          in the input file 
            %ppobj   - PostProcess object
            %
            %Outputs:
            %--------
            %There are no outputs.  Properties are stored in the ppTABLE
            %handle class.
            %--------------------------------------------------------------
            
            %Extract variables from VARtxt and evaluate them in the current
            %workspace.
            if ~isempty(VARtxt)                
                for i = 1:length(VARtxt)
                    eval([VARtxt{i},';']);
                end
            end
            
            %Determine id of current FIGURE object
            obj.id = checkPostInput(TABtext,1,'id',obj.id);

            %Determine subtable configuration in this FIGURE.
            obj.subtableNum = checkPostInput(TABtext,1,'subtableNum',obj.subtableNum);

            %Determine table properies
            obj.settab = checkPostInput(TABtext,1,'settab',obj.settab); 
            
            %Determine norm type to use
            obj.normType = checkPostInput(TABtext,1,'normType',obj.normType);

            %Determine number of rows and columns
            obj.numRows = checkPostInput(TABtext,1,'numRows',obj.numRows);
            obj.numCols = checkPostInput(TABtext,1,'numCols',obj.numCols);
            obj.numElem = obj.numRows*obj.numCols;
            
            %Determine data
            obj.rowKey = checkPostInput(TABtext,0,'rowData',obj.rowKey);
            obj.colKey = checkPostInput(TABtext,0,'colData',obj.colKey);
            obj.elemKey = checkPostInput(TABtext,0,'elemData',obj.elemKey);
            
            checkRowColElemData(obj.rowKey,obj.numRows,[],'rowData');
            checkRowColElemData(obj.colKey,[],obj.numCols,'colData');
            checkRowColElemData(obj.elemKey,obj.numRows,obj.numCols,'elemData');
            
            %Determine models to use in plots
            obj.model = checkPostInput(TABtext,obj.numRows,'model',obj.model);
            obj.modelAux = checkPostInput(TABtext,obj.numRows,'modelAux',obj.modelAux);
            obj.modelFom = checkPostInput(TABtext,obj.numRows,'modelFom',obj.modelFom);
            
            %Error checks!
            A = isempty(obj.rowKey); B = isempty(obj.colKey); C = isempty(obj.elemKey);
            if A && B && C
                error('In PP file, must specify either colData, rowData, or elemData');
            elseif (~A && ~B && C) || (A && ~B && ~C) || (~A && B && ~C) || (~A && ~B && ~C)
                error('In PP file, may only specify one of the following fields: rowData, colData, elemData.  Others must be empty');
            elseif ~A
                checkRowColElemData(obj.model,obj.numRows,[],'model');
                checkRowColElemData(obj.modelAux,obj.numRows,[],'modelAux');
                checkRowColElemData(obj.modelFom,obj.numRows,[],'modelFom');
            elseif ~B
                checkRowColElemData(obj.model,[],obj.numCols,'model');
                checkRowColElemData(obj.modelAux,[],obj.numCols,'modelAux');
                checkRowColElemData(obj.modelFom,[],obj.numCols,'modelFom');
            elseif ~C
                checkRowColElemData(obj.model,obj.numRows,obj.numCols,'model');
                checkRowColElemData(obj.modelAux,obj.numRows,obj.numCols,'modelAux');
                checkRowColElemData(obj.modelFom,obj.numRows,obj.numCols,'modelFom');
            end
            
            %Determine title, row and column names
            obj.title = checkPostInput(TABtext,1,'title',obj.title);
            obj.rowNames = checkPostInput(TABtext,1,'rowNames',obj.rowNames);
            obj.colNames = checkPostInput(TABtext,1,'colNames',obj.colNames);
                  
            %Extract the data for filling the entries of the table
            if ~isempty(obj.rowKey)
                m = obj.numRows;
                keys = obj.rowKey;
                
                for i = 1:m
                    modelobj = determineModelFromText(ppobj,obj.model,i,1,false);
                    modelAuxobj = determineModelFromText(ppobj,obj.modelAux,i,1,true);
                    modelFomobj = determineModelFromText(ppobj,obj.modelFom,i,1,false);
                    
                    temp = extractDataFromKey(obj,keys{i,1},modelobj,modelAuxobj,modelFomobj);
                    obj.rowData{i,1} = temp(:)';
                end
                data = cell2mat(obj.rowData);
            elseif ~isempty(obj.colKey)
                n = obj.numCols;
                keys = obj.colKey;
                
                for j = 1:n
                    modelobj = determineModelFromText(ppobj,obj.model,1,j,false);
                    modelAuxobj = determineModelFromText(ppobj,obj.modelAux,1,j,true);
                    modelFomobj = determineModelFromText(ppobj,obj.modelFom,1,j,false);
                    
                    temp = extractDataFromKey(obj,keys{1,j},modelobj,modelAuxobj,modelFomobj);
                    obj.colData{1,j} = temp(:);
                end
                data = cell2mat(obj.colData);
            elseif ~isempty(obj.elemKey)
                m = obj.numRows;
                n = obj.numCols;
                
                keys = convertColRow2ElemData(obj.elemKey,m,n);
                obj.model = convertColRow2ElemData(obj.model,m,n);
                obj.modelAux = convertColRow2ElemData(obj.modelAux,m,n);
                obj.modelFom = convertColRow2ElemData(obj.modelFom,m,n);

                for i = 1:m
                    for j = 1:n
                        modelobj = determineModelFromText(ppobj,obj.model,i,j,false);
                        modelAuxobj = determineModelFromText(ppobj,obj.modelAux,i,j,true);
                        modelFomobj = determineModelFromText(ppobj,obj.modelFom,i,j,false);

                        obj.elemData{i,j} = extractDataFromKey(obj,keys{i,j},modelobj,modelAuxobj,modelFomobj);
                    end
                end
                data = cell2mat(obj.elemData);
            end
        end
    end
end
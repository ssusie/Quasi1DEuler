classdef ppFIGURE < handle
    
    properties (SetAccess = private, GetAccess = public)
        %Input Properties
        han = [];
        ppaxHan = [];
        pptabHan = [];
        
        id   = 1;
        axesLayout = [1,1];
        setfig = {'Position',[322   541   831   565]};
        fname = [];
    end
    
    methods
        function  [obj] = ppFIGURE(fname,ind,ppobj)
            %This is the constructor of the ppFIGURE class.  It reads
            %the appropriate entries from the input file, stores them in
            %the class instance, and generates the requested plots.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %fname - string containing the postprocessing file name 
            %ind   - scalar indicating the index of the desired FIGURE
            %        object that will come out of readInFile
            %ppobj - PostProcess object
            %
            %Outputs:
            %--------
            %obj   - instance of the PostProcess object that was
            %        constructed 
            %--------------------------------------------------------------
            
            %Extract parameters from CFG file
            VARtext = readInFile('VAR',fname,1);
            [FIGcell,AXcell] = readInFile('AXES',fname,[]);
            [~,TABcell] = readInFile('TABLE',fname,[]);
            
            obj.determineFIG(FIGcell{ind,1},[ppobj.GVARtxt;VARtext]);
            
            obj.han = figure('Position',[322   541   831   565]);
            set(obj.han,'name',[fname,': id #',num2str(obj.id)]);
            %Set the user specified figure properties as long as the cell
            %array is not empty
            if ~isempty(obj.setfig)
                set(obj.han,obj.setfig{:});
            end
            
            n = size(AXcell,2);
            m = size(TABcell,2);
            
            obj.ppaxHan = [];
            obj.pptabHan = [];
            
            for i = 1:n
                if ~isempty(AXcell{ind,i}) 
                    obj.ppaxHan  = [obj.ppaxHan; ppAXES(fname,ind,i,obj.axesLayout,ppobj)];
                end
            end
            
            for j = 1:m
                if ~isempty(TABcell{ind,j}) 
                    obj.pptabHan = [obj.pptabHan; ppTABLE(fname,ind,j,obj.axesLayout,ppobj)];
                end
            end
            
            if ~isempty(obj.fname)
                for i = 1:size(obj.fname,1)
                    tmp = regexp(obj.fname{i},'\.','end');
                    ext = obj.fname{i}(tmp(end)+1:end);
                    
                    if strcmp(ext,'fig')
                        saveas(obj.han,obj.fname{i});
                    else
                        print(obj.han,['-d',ext],obj.fname{i});
                    end
                end
            end
        end
            
        function  [] = determineFIG(obj,FIGtext,VARtxt)
            %This function converts the text from the FIGURE object into a
            %structure whose fields correspond to the ppFIGURE properties
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj     - ppFIGURE object
            %FIGtext - a string containing the text from the FIGURE objects
            %          in the input file 
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

            %Determine id of current FIGURE object
            obj.id = checkPostInput(FIGtext,1,'id',obj.id); 
            
            %Determine subplot configuration in this FIGURE.
            obj.axesLayout = checkPostInput(FIGtext,2,'axesLayout',obj.axesLayout); 
            
            %Determine figure properies
            obj.setfig = checkPostInput(FIGtext,1,'setfig',obj.setfig); 
            
            %Determine filename to save figure (can't use checkPostInput
            %because it can have an arbitrary number of rows)
            tmp = determinePropText(FIGtext,'fname');
            if ~isempty(tmp)
                fnameVal = eval(tmp);
                if ~iscell(fnameVal)
                    obj.fname = {fnameVal};
                else
                    obj.fname = fnameVal;
                end
            end
        end
    end
end
classdef PostProcess < handle
    
    properties (SetAccess = private, GetAccess = public)
        %Input Properties
        ppfigHan = [];
        
        specialPlot = [];
        numPlotObj = [];
        
        fom = [];
        rom = [];
        gnat = [];
        tpwl = [];
        prob = [];
    end
    
    properties (Hidden = true, SetAccess = private, GetAccess = public)
        GVARtxt;
    end
    
    methods
        function  [obj] = PostProcess(fname,fomobj,romobj,gnatobj,tpwlobj)
            %This is the constructor of the PostProcess class.  It reads
            %the appropriate entries from the input file, stores them in
            %the class instance, and generates the requested plots.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %fname   - string containing the postprocessing file name 
            %fomobj  - FOM object
            %romobj  - vector of ROM objects
            %gnatobj - vector of GNAT or locGNAT objects
            %tpwlobj - vector of TPWL objects
            %
            %Outputs:
            %--------
            %obj     - instance of the PostProcess object that was
            %          constructed 
            %--------------------------------------------------------------
            
            obj.fom  = fomobj;
            obj.rom  = romobj;
            obj.gnat = gnatobj;
            obj.tpwl = tpwlobj;
            
            %%%This section handles the global variable issue%%%
            obj.GVARtxt = [];
            %Add global variable from all fom's (if these come from
            %different CFG files, they will not necessarily be the same)
            %and all rom's, gnat's, tpwl's (in case their foms (they have
            %the same global variables as the fom) come from different CFG
            %files)
            for i = 1:length(fomobj)
                obj.GVARtxt = [obj.GVARtxt; fomobj{i}.GVARtxt];
            end
            
            for i = 1:length(romobj)
               
                obj.GVARtxt = [obj.GVARtxt; romobj{i}.GVARtxt];
            end
            
            for i = 1:length(gnatobj)
                obj.GVARtxt = [obj.GVARtxt; gnatobj{i}.GVARtxt];
            end
            
            for i = 1:length(tpwlobj)
                obj.GVARtxt = [obj.GVARtxt; tpwlobj{i}.GVARtxt];
            end
%             obj.GVARtxt = unique(obj.GVARtxt);
            %%%End of global variable section%%%
            
            [FIGcell,~] = readInFile('AXES',fname,[]);  
            
            obj.ppfigHan = [];
            for i = 1:length(FIGcell)
                obj.ppfigHan = [obj.ppfigHan; ppFIGURE(fname,i,obj)];
            end
        end
    end
end
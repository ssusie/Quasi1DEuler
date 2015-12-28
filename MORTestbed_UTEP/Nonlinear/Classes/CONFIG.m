classdef CONFIG < handle
    
    properties (SetAccess = private, GetAccess = public)
        %Input Properties
        problem;
        id;
        type; %Configuration type: 'train' or 'online'
        desc;
        time;
        form = 'non_descriptor';
        altFile = [];
        DLim = [];
        nNodes; 
        ndof; %Number of unknown dofs
        param; %Parameters
        
        staticFlag=false;
        BFunc = @() 0; %Function name or function handle describing input matrix
        CFunc = @() 0; %Function name or function handle describing output matrix
        DFunc = @() 0; %Function name or function handle describing feedthrough matrix
        GFunc = @(x,y,z) 0; %Function name or function handle describing source term (fcn of x,y,z)
        icFunc = @(x,y,z) 0; %Function name or function handle describing initial condition (fcn of x,y,z)
        inFunc = @(t) 0; %Function name or function handle describing input function (fcn of t)
    end
    
    properties (Hidden = true, SetAccess = private, GetAccess = public)
        %Global variable text
        GVARtxt = [];
        defaultFile = [];
    end
    
    methods
        function  [obj] = CONFIG(CFGfile,cfgid)
            %This is the constructor of the CONFIG class.  It reads the
            %appropriate entries from the input file and stores them in the
            %class instance.  
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %CFGfile - string indicating the filename of the cfg file to
            %          use
            %cfgid   - integer indicating the CONFIG object in CFGfile to
            %          use (corresponds to the id number of the CONFIG
            %          block)
            %
            %Outputs:
            %--------
            %obj     - instance of the CONFIG object that was constructed
            %--------------------------------------------------------------
            
            %Extract parameters from CFG file
            obj.GVARtxt = readInFile('GVAR',CFGfile,1);
            VARtext = readInFile('VAR',CFGfile,1);
            %Extract the FOM text from the cfg file
            FOMtext = readInFile('FOM',CFGfile,cfgid);
            %Extract the CONFIG text with id number = cfgid from cfg file
            CONFIGtext = readInFile('CONFIG',CFGfile,cfgid);
            
            txt = [FOMtext; CONFIGtext];
            %Determine the problem (apriori) so we can set all of the
            %defaults properly.
            probNum = extractModelPropMultInput(txt,1,1,'problem',[]);
            %Make sure user specifed problem number
            if isempty(probNum)
               error(['In ',CFGfile,': The problem number was not specified.  This field must be specified for program to have meaning.']); 
            end
            
            %Initialize time structure
            obj.time = struct('T',[],'dt',[],'nstep',[],'quiet',[]);
            
%             obj.setAllValues2Default(probNum);
            
            %Determine the CONFIG properties based on the text in the CONFIG
            %section of the input file
            determineCONFIG(obj,txt,[obj.GVARtxt;VARtext]);
            
        end
        
        function  [] = determineCONFIG(obj,CONFIGtext,VARtxt)
            %This function computes and stores properties of CONFIG class
            %from the char array CONFIGtext
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj       - instance of CONFIG class
            %CONFIGtxt - a string containing the text from the FOM and
            %            CONFIG objects in the input file
            %VARtxt    - a string containing the text from the VAR block in
            %            the CFG file  
            %
            %Outputs:
            %--------
            %No outputs.  Properties are stored in CONFIG handle class.
            %--------------------------------------------------------------
            
            %Extract variables from VARtxt and evaluate them in the current
            %workspace.
            if ~isempty(VARtxt)                
                for i = 1:length(VARtxt)
                    eval([VARtxt{i},';']);
                end
            end
            
            %%%%%%%%% Determine problem number %%%%%%%%%
            %Determine the problem number specified in the input file
            %At this point in the program, I know that the below line has
            %meaning since the default file has problem defined and so does
            %the input file (I already checked outside this function).
            obj.problem = extractInputRobust(CONFIGtext,'problem',[]);
             
            %%%%%%%%% Determine id %%%%%%%%%
            obj.id = extractModelPropMultInput(CONFIGtext,1,1,'id',[]);
            if isempty(obj.id)
                error('In CFG file, the ID must be specified.  This controls the connection between different aspects of the program.');
            end
            
            %%%%%%%%% Determine form %%%%%%%%%
            tmp = determinePropText(CONFIGtext,'form');
            if ~isempty(tmp)
                obj.form = tmp;
            end
            
            %%%%%%%%% Determine staticFlag %%%%%%%%%
            obj.staticFlag = determinePropText(CONFIGtext,'staticFlag');

            %%%%%%%%% Determine desc %%%%%%%%%
            obj.desc = determinePropText(CONFIGtext,'desc');
            
            %%%%%%%%% Determine type %%%%%%%%%
            obj.type = determinePropText(CONFIGtext,'type');
            
            %%%%%%%%% Determine number of dofs in FOM %%%%%%%%%
            obj.ndof = extractInputRobust(CONFIGtext,'ndof',obj.ndof);
            if length(obj.ndof) > 1
                error('ndof must be a scalar indicating the number of unknowns in the problem');
            end

            %%%%%%%%% Determine time structure %%%%%%%%%
            obj.time = struct('T',[],'dt',[],'nstep',[],'quiet',[]);
            obj.time.T     = extractInputRobust(CONFIGtext,'T',[]);
            obj.time.dt    = extractInputRobust(CONFIGtext,'dt',[]);
            obj.time.nstep = extractInputRobust(CONFIGtext,'nstep',1000);
            obj.time.quiet = extractInputRobust(CONFIGtext,'timeQuiet',false);
            obj.time.steadyconverge = extractInputRobust(CONFIGtext,'steadyconverge',false);
            obj.time.cfl = extractInputRobust(CONFIGtext,'cfl',[]);
            
%             %Compute missing time quantity
%             if isempty(obj.time.T) && ~isempty(obj.time.dt) && ~isempty(obj.time.nstep) %Calculate T from dt, nstep
%                 obj.time.T = [0,obj.time.nstep*obj.time.dt]; %Assume start time is zero
%             elseif ~isempty(obj.time.T) && isempty(obj.time.dt) && ~isempty(obj.time.nstep) %Calculate dt from T, nstep
%                 obj.time.dt = (obj.time.T(2) - obj.time.T(1))/obj.time.nstep;
%             elseif ~isempty(obj.time.T) && ~isempty(obj.time.dt) && isempty(obj.time.nstep) %Calculate nstep from T, dt
%                 obj.time.nstep = ceil((obj.time.T(2) - obj.time.T(1))/obj.time.dt);
%             else
%                 %User needs to exactly 2 of the 3 time quantities
%                 error('Must specify exactly 2 of the 3 following fields: T, dt, nstep');
%             end
            
            %%%%%%%%% Determine number of nodes in the FOM %%%%%%%%%
            obj.nNodes = extractInputRobust(CONFIGtext,'nNodes',[]); 
            
            %%%%%%%%% Determine the limits of the domain %%%%%%%%%
            obj.DLim = extractInputRobust(CONFIGtext,'DLim',[]);

            %%%%%%%%% Determine alternate file if it exists %%%%%%%%% 
            obj.altFile = determinePropText(CONFIGtext,'altFile');
            if ~isempty(obj.altFile)
                if ~strcmpi(obj.altFile(1:2),'[]') && ~strcmpi(obj.altFile(end-1:end),'.m')
                    error('All file names must end with a .m extension');
                end
                if strcmpi(obj.altFile(end-1:end),'.m')
                    obj.altFile = obj.altFile(1:end-2);
                end
            end
            
            %%%%%%%%% Determine and store the input matrix function %%%%%%%%% 
            temp = determinePropText(CONFIGtext,'BFunc');
            if ~isempty(temp) && isempty(strfind(temp,'.m'))
                %This is the case where a valid matlab statement is
                %specified inline
                obj.BFunc = eval(['@() ',temp]);
            elseif isempty(temp)
                %If BFunc is empty, we want to use the default
                %which has already been set, so we do nothing here
                obj.BFunc = @() 0;
            else
                %This is the case where an input file is specified
                obj.BFunc = temp(1,1:end-2);
            end
            
            %%%%%%%%% Determine and store the output matrix function %%%%%%%%% 
            temp = determinePropText(CONFIGtext,'CFunc');
            if ~isempty(temp) && isempty(strfind(temp,'.m'))
                %This is the case where a valid matlab statement is
                %specified inline
                obj.CFunc = eval(['@() ',temp]);
            elseif isempty(temp)
                %If CFunc is empty, we want to use the default
                %which has already been set, so we do nothing here
                obj.CFunc = @() 0;
            else
                %This is the case where an input file is specified
                obj.CFunc = temp(1,1:end-2);
            end
            
            %%%%%%%%% Determine and store the feedthrough matrix function %%%%%%%%% 
            temp = determinePropText(CONFIGtext,'DFunc');
            if ~isempty(temp) && isempty(strfind(temp,'.m'))
                %This is the case where a valid matlab statement is
                %specified inline
                obj.DFunc = eval(['@() ',temp]);
            elseif isempty(temp)
                %If DFunc is empty, we want to use the default
                %which has already been set, so we do nothing here
                obj.DFunc = @(x,y,z) 0;
            else
                %This is the case where an input file is specified
                obj.DFunc = temp(1,1:end-2);
            end
            
            %%%%%%%%% Determine the source term function of x, y, z %%%%%%%%% 
            temp = determinePropText(CONFIGtext,'GFunc');
            if ~isempty(temp) && isempty(strfind(temp,'.m'))
                %This is the case where a function of x, y, z is specified
                %inline
                obj.GFunc = eval(['@(x,y,z) ',temp]);
            elseif isempty(temp)
                %If GFunc is empty, we want to use the default
                %which has already been set, so we do nothing here
                obj.GFunc = @(x,y,z) 0;
            else
                %This is the case where an input file is specified
                obj.GFunc = temp(1,1:end-2);
            end
            
            %%%%%%%%% Determine the initial condition function of x, y, z %%%%%%%%% 
            temp = determinePropText(CONFIGtext,'icFunc');
            if ~isempty(temp) && isempty(strfind(temp,'.m'))
                %This is the case where a function of x, y, z is specified
                %inline
                obj.icFunc = eval(['@(x,y,z) ',temp]);
            elseif isempty(temp)
                %If icFunc is empty, we want to use the default
                %which has already been set, so we do nothing here
                obj.icFunc = @(x,y,z) 0;
            else
                %This is the case where an input file is specified
                obj.icFunc = temp(1,1:end-2);
            end

            %%%%%%%%% Determine the parameters for the problem at hand %%%%%%%%%
            obj.param = extractInputRobust(CONFIGtext,'param',obj.param);

            %%%%%%%%% Determine the input file name %%%%%%%%%
            temp = determinePropText(CONFIGtext,'inFunc');
            if ~isempty(temp) && isempty(strfind(temp,'.m'))
                %This is the case where a function of t is specified
                %inline
                obj.inFunc = eval(['@(t) ',temp]);
            elseif isempty(temp)
                %If icFunc is empty, we want to use the default
                %which has already been set, so we do nothing here
                obj.inFunc = @(t) 0;
            else
                %This is the case where an input file is specified
                obj.inFunc = temp(1,1:end-2);
            end
            
            %if ~isempty(temp)
            %  obj.inFunc = temp;
            %end
            
            %Remove .m extension if it exists from input function
            if strcmpi(obj.inFunc(end-1:end),'.m')
                obj.inFunc(end-1:end) = [];
            end
        end
        
        function  [] = setAllValues2Default(obj,probNum)
            %This function sets all values to default and then goes through
            %and overwrites the values the user specifies.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %
            %Outputs:
            %--------
            %
            %--------------------------------------------------------------

            %Get path and filename of the default input file
            folder = 'usr';
            tmp = pwd;
            ind = regexp(tmp,folder);
            dir2add = tmp(1:ind+length(folder)-1);
            fname = mapProbNum2DefaultFile(probNum);
            obj.defaultFile = [dir2add,SlashBS,'defaults',SlashBS,fname];
            fname = [obj.defaultFile,'.cfg'];
            
            %Extract the relevant properties and set the properties
            VARtext = readInFile('VAR',fname,1);
            %Extract the FOM text from the cfg file
            FIMtext = readInFile('FOM',fname,1);
            %Extract the CONFIG text with id number = cfgid from cfg file
            CONFIGtext = readInFile('CONFIG',fname,1);
            
            %Determine the CONFIG properties based on the text in the CONFIG
            %section of the input file
            determineCONFIG(obj,[FOMtext; CONFIGtext],VARtext);
            
            %All of the default values are now set.
        end
        
        function [] = modifyTimeProperties(obj,configobj)
            
           obj.time = configobj.time;    
        end
        
        function  [] = setProperties(obj,prop,val)
            obj.(prop) = val;
        end
        
        function  [] = setTimeObj(obj,prop,val)
            obj.time.(prop) = val;
        end
    end  
end

%             %Determine the problem number specified in the input file
%             %At this point in the program, I know that the below line has
%             %meaning since the default file has problem defined and so does
%             %the input file (I already checked outside this function).
%             obj.problem = extractModelPropMultInput(CONFIGtext,1,1,'problem',[]);
%             %Determine and store configuration id
%             tmp = determinePropText(CONFIGtext,'id');
%             if isempty(tmp)
%                 error('In CFG file, the ID must be specified.  This controls the connection between different aspects of the program.');
%             end
%             obj.id = str2double(tmp);
%             %Determine and store configuration type
%             temp = determinePropText(CONFIGtext,'type');
%             if ~isempty(temp)
%                 obj.type = temp;
%             end
%             %Determine and store configuration description
%             temp = determinePropText(CONFIGtext,'desc');
%             if ~isempty(temp)
%                 obj.desc = temp;
%             end
%             
%             %Determine the upper and lower bounds on the time interval
%             obj.time.T = extractModelPropMultInput(CONFIGtext,1,1,'T',obj.time.T);
%             %Determine the time step size
%             obj.time.dt = extractModelPropMultInput(CONFIGtext,1,1,'dt',obj.time.dt);
%             %Determine the number of time steps
%             obj.time.nstep = extractModelPropMultInput(CONFIGtext,1,1,'nstep',obj.time.nstep);
%             %Determine whether or not to time progress
%             obj.time.quiet = extractModelPropMultInput(CONFIGtext,1,1,'timeQuiet',obj.time.quiet); 
%             
%             %Compute missing time quantity
%             if isempty(obj.time.dt) && isempty(obj.time.nstep)
%                 %User needs to specify either nstep or dt
%                 error('Must specify either time increment or number of time steps');
%             elseif isempty(obj.time.dt) && ~isempty(obj.time.nstep) %Calculate dt from nstep
%                 obj.time.dt = (obj.time.T(2) - obj.time.T(1))/obj.time.nstep;
%             elseif isempty(obj.time.nstep) && ~isempty(obj.time.dt) %Calculate nstep from dt
%                 obj.time.nstep = ceil((obj.time.T(2) - obj.time.T(1))/obj.time.dt);
%             else %if both are specified, return an error
%                 if obj.time.dt ~= (obj.time.T(2) - obj.time.T(1))/obj.time.nstep
%                     error('time increment (dt) and number of time steps (nstep) both specified and are not consistent (dt must equal (T(2) - T(1))/nstep)');
%                 end
%             end
% 
%             %Determine alternate file if it exists
%             temp = determinePropText(CONFIGtext,'altFile');
%             if ~isempty(temp)
%                 obj.altFile = temp;
%                 if ~isempty(obj.altFile) && ~strcmpi(obj.altFile(1:2),'[]') && ~strcmpi(obj.altFile(end-1:end),'.m')
%                     error('All file names must end with a .m extension');
%                 end
%                 if strcmpi(obj.altFile(end-1:end),'.m')
%                     obj.altFile = obj.altFile(1:end-2);
%                 end
%             end
%             %Determine and store the number of nodes in the FOM
%             obj.ndof = extractModelPropMultInput(CONFIGtext,1,1,'ndof',obj.ndof); 
%             if length(obj.ndof) > 1
%                 error('ndof must be a scalar indicating the number of unknowns in the problem');
%             end
%             %Determine and store the number of nodes in the FOM
%             obj.nNodes = extractModelPropMultInput(CONFIGtext,1,1,'nNodes',obj.nNodes); 
%             %Determine the limits of the domain
%             temp = determinePropText(CONFIGtext,'DLim');
%             if ~isempty(temp)
%                 obj.DLim = eval(temp);
%             end
%             %Determine and store the input matrix function
%             temp = determinePropText(CONFIGtext,'BFunc');
%             if ~isempty(temp) && isempty(strfind(temp,'.m'))
%                 %This is the case where a valid matlab statement is
%                 %specified inline
%                 obj.BFunc = eval(['@() ',temp]);
%             elseif isempty(temp)
%                 %If BFunc is empty, we want to use the default
%                 %which has already been set, so we do nothing here
%             else
%                 %This is the case where an input file is specified
%                 obj.BFunc = temp(1,1:end-2);
%             end
%             %Determine and store the output matrix function
%             temp = determinePropText(CONFIGtext,'CFunc');
%             if ~isempty(temp) && isempty(strfind(temp,'.m'))
%                 %This is the case where a valid matlab statement is
%                 %specified inline
%                 obj.CFunc = eval(['@() ',temp]);
%             elseif isempty(temp)
%                 %If CFunc is empty, we want to use the default
%                 %which has already been set, so we do nothing here
%             else
%                 %This is the case where an input file is specified
%                 obj.CFunc = temp(1,1:end-2);
%             end
%             %Determine and store the feedthrough matrix function
%             temp = determinePropText(CONFIGtext,'DFunc');
%             if ~isempty(temp) && isempty(strfind(temp,'.m'))
%                 %This is the case where a valid matlab statement is
%                 %specified inline
%                 obj.DFunc = eval(['@() ',temp]);
%             elseif isempty(temp)
%                 %If DFunc is empty, we want to use the default
%                 %which has already been set, so we do nothing here
%             else
%                 %This is the case where an input file is specified
%                 obj.DFunc = temp(1,1:end-2);
%             end
%             %Determine the source term function of x, y, z
%             temp = determinePropText(CONFIGtext,'GFunc');
%             if ~isempty(temp) && isempty(strfind(temp,'.m'))
%                 %This is the case where a function of x, y, z is specified
%                 %inline
%                 obj.GFunc = eval(['@(x,y,z) ',temp]);
%             elseif isempty(temp)
%                 %If GFunc is empty, we want to use the default
%                 %which has already been set, so we do nothing here
%             else
%                 %This is the case where an input file is specified
%                 obj.GFunc = temp(1,1:end-2);
%             end
%             %Determine the initial condition function of x, y, z
%             temp = determinePropText(CONFIGtext,'icFunc');
%             if ~isempty(temp) && isempty(strfind(temp,'.m'))
%                 %This is the case where a function of x, y, z is specified
%                 %inline
%                 obj.icFunc = eval(['@(x,y,z) ',temp]);
%             elseif isempty(temp)
%                 %If icFunc is empty, we want to use the default
%                 %which has already been set, so we do nothing here
%             else
%                 %This is the case where an input file is specified
%                 obj.icFunc = temp(1,1:end-2);
%             end
% 
%             %Determine the parameters for the problem at hand
%             obj.param = extractInputRobust(CONFIGtext,'param',obj.param);
%             %obj.param = extractModelPropMultInput(CONFIGtext,1,1,'param',obj.param); 
%             %Determine the input file name
%             temp = determinePropText(CONFIGtext,'inFunc');
%             if ~isempty(temp) && isempty(strfind(temp,'.m'))
%                 %This is the case where a function of t is specified
%                 %inline
%                 obj.inFunc = eval(['@(t) ',temp]);
%             elseif isempty(temp)
%                 %If icFunc is empty, we want to use the default
%                 %which has already been set, so we do nothing here
%             else
%                 %This is the case where an input file is specified
%                 obj.inFunc = temp(1,1:end-2);
%             end
%             
%             %if ~isempty(temp)
%             %  obj.inFunc = temp;
%             %end
%             
%             %Remove .m extension if it exists from input function
%             if strcmpi(obj.inFunc(end-1:end),'.m')
%                 obj.inFunc(end-1:end) = [];
%             end
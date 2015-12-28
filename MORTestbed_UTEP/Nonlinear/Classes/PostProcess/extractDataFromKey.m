function  [data] = extractDataFromKey(obj,key,modelobj,modelAuxobj,modelFomobj)
%This function extracts the appropriate data from the model
%indicated by the string model.
%--------------------------------------------------------------------------
%Inputs:
%-------
%obj         - ppAXES or ppTABLE object
%key         - string containing the keyword used in
%              determining the data to extract
%modelobj    - model object to pull data from.
%modelAuxobj - model object to use as the auxillary model (only needed when
%              modelobj is of class GNAT or locGNAT) 
%modelFomobj - FOM object to use (i.e. which FOM are the ROMs supposed to
%              match)  
%
%Outputs:
%--------
%data        - vector or scalar containing the information that was
%              requested in the input file. 
%--------------------------------------------------------------------------

if isempty(key)
    data = [];
    return;
end

%If the user requests a animation, parse the input
if length(key) >= 7 && strcmpi(key(1:7),'animate')
    key = key(1:7);
    if length(key) >= 8
        animNum = str2double(key(8:end));
    else
        animNum = 1;
    end
end


%If the user requests a model property, appropriately parse the input
%statement
if length(key) > 3 && strcmpi(key(1:4),'prop')
    key = key(6:end);
    
    %Handle the case where the user wants to extract from a structure
    propEndInd = regexp(key,'\.','once');
    if isempty(propEndInd)
        structField = [];
    else
        structField = key(propEndInd:end);
        key = key(1:propEndInd-1);
    end
    
    %Handle the case where the user wants to extract from a cell
    propEndInd = regexp(key,'[0-9]');
    if isempty(propEndInd)
        index = [];
    else
        if sum(abs(diff([propEndInd,length(key)+1]) - 1)) ~= 0
            %To be valid, all numbers must be specified at the END of the
            %property name.  This statement determines if there are any
            %numbers in the property name not consecutive and not at the
            %end.
            error('There was an error in your postprocessing request in the xData or yData field.  There cannot be a property name with a number anywhere but the end.');
        end
        index = str2double(key(propEndInd:end));
        key = key(1:propEndInd-1);
    end
    
    
    prop = properties(modelobj);
    flag = false;
    for i = 1:length(prop)
        if strcmpi(key,prop{i})
            flag = true;
            break;
        end
    end
    
    if flag
        data = eval(['modelobj.',key,structField]);
    else
        data = NaN;
    end
    
    if ~isempty(index)
       %If the user requested a specific entry from the data
       if iscell(data)
           data = data{index};
       else
           data = data(index);
       end
    end
end

%If the user requests a state vector, appropriately parse
%the input statement
if length(key) > 1 && strcmpi(key(1:2),'sv')
    %Determine which unknown to use (if multiple are
    %stacked on each other to form the state vector)
    ind = regexp(key,'[0-9][@]');
    if isempty(ind)
        SVid = 1;
    else
        SVid = str2double(key(3:ind));
    end
    
    %Determine the position or time to extract the state and
    %alter key into a usable form
    indKey = regexp(key,'[^0-9@.=]');
    indNum = regexp(key,'[=]');
    num = eval(key(indNum+1:end));
    key = ['sv_',key(indKey(3))];
end

%If the user requests a pod vector, appropriately parse
%the input statement
if length(key) > 2 && strcmpi(key(1:3),'pod')
    %Need to check if modelobj is ROM, GNAT, or TPWL.
    
    %Determine which unknown to use (if multiple are
    %stacked on each other to form the state vector)
    ind = regexp(key,'[0-9][@]');
    if isempty(ind)
        SVid = 1;
    else
        SVid = str2double(key(4:ind));
    end
    
    %Determine the position or time to extract the state and
    %alter key into a usable form
    b_indNum = regexpi(key,'b=(\d+)','tokenExtents');
    n_indNum = regexpi(key,'n=(\d+)','tokenExtents');
    
    bnum = eval(key(unique(b_indNum{1})));
    nnum = eval(key(unique(n_indNum{1})));
    
    if round(bnum) ~= bnum
        %Need to issue warning of the change that is taking place!
        bnum = round(bnum);
        if bnum == 0
            bnum = 1;
        elseif bnum == modelobj.nBases
            bnum = modelobj.nBases;
        end
    end
    
    if round(nnum) ~= nnum
        %Need to issue warning of the change that is taking place!
        nnum = round(nnum);
        if nnum == 0
            nnum = 1;
        elseif nnum == modelobj.nY
            nnum = modelobj.nY;
        end
    end
    key = 'pod';
end

%If the user requests a cluster center vector, appropriately parse
%the input statement
if length(key) > 9 && strcmpi(key(1:10),'clusCenter')
    %Need to check if modelobj is ROM (nBases > 1) or locGNAT.
    
    %Determine which unknown to use (if multiple are
    %stacked on each other to form the state vector)
    ind = regexp(key,'[0-9][@]');
    if isempty(ind)
        SVid = 1;
    else
        SVid = str2double(key(4:ind));
    end
    
    %Determine the position or time to extract the state and
    %alter key into a usable form
    indNum = regexp(key,'[=]');
    num = eval(key(indNum+1:end));
    if round(num) ~= num
        %Need to issue warning of the change that is taking place!
        num = round(num);
        if num == 0
            num = 1;
        elseif num == modelobj.nBases
            num = modelobj.nBases;
        end
    end
    key = 'clusCenter';
end

%Prepare key and SVid for xdomain or ydomain data request
if ~isempty(strfind(key,'xdomain')) || ~isempty(strfind(key,'ydomain'))
    SVid = str2double(key(8:end));
    key = key(1:7);
end

%If percent was requested for error, multiply by 100
percentFlag = false;
if ~isempty(findstr(key,'err')) && ~isempty(findstr(key,'%'))
    ind = findstr(key,'%');
    percentFlag = true;
    key(ind) = [];
end

%Extract the appropriate data depending on the value of key
if ~ischar(key)
    %If user requests to plot literals
    data = key;
    return;
end

switch key
    case 'animate'
        data = @(pstr,freq) modelobj.prob.problemAnimate(modelobj,modelAuxobj,pstr,animNum,freq);
    case 'sv_t'
        %This won't work for 2D problems!
        [~,~,temp] = modelobj.prob.returnPosSV(modelobj,modelAuxobj,SVid,'sv');
        
        t = modelobj.time.T(1):modelobj.time.dt:modelobj.time.T(2);
        data = interp1(t',temp',num);
        data = data';
    case 'sv_x'
        %This won't work for 2D problems!
        [x,~,svTemp] = modelobj.prob.returnPosSV(modelobj,modelAuxobj,SVid,'sv');
        data = interp1(x,svTemp,num);
    case 'pod'
        [~,~,temp] = modelobj.prob.returnPosSV(modelobj,modelAuxobj,SVid,'pod',bnum);
        data = temp(:,nnum);
    case 'clusCenter'
        [~,~,temp] = modelobj.prob.returnPosSV(modelobj,modelAuxobj,SVid,'clusCenter');
        data = temp(:,num);
    case 'time'
        data = modelobj.time.T(1):modelobj.time.dt:modelobj.time.T(2);
    case 'timestep'
        data = 1:modelobj.time.nstep;
    case 'xdomain'
        [data,~,~] = modelobj.prob.returnPosSV(modelobj,modelAuxobj,SVid,'sv');
    case 'ydomain'
        [~,data,~] = modelobj.prob.returnPosSV(modelobj,modelAuxobj,SVid,'sv');
    case 'output'
        svF = modelobj.reconstructFullState(modelAuxobj);
        data = modelobj.prob.C*svF;
    case 'ontime'
        data = modelobj.ontime;
    case 'speedup'
        data = modelFomobj.ontime/modelobj.ontime;
    case 'Aerr'
        svROM = modelobj.reconstructFullState(modelAuxobj);
        svFOM = modelFomobj.sv;
        data = computeError(svROM,svFOM,obj.normType,'abs');
    case 'MAerr'
        svROM = modelobj.reconstructFullState(modelAuxobj);
        svFOM = modelFomobj.sv;
        tmp = computeError(svROM,svFOM,obj.normType,'abs');
        data = max(tmp(~isnan(tmp)));
    case 'AvgAerr'
        svROM = modelobj.reconstructFullState(modelAuxobj);
        svFOM = modelFomobj.sv;
        tmp = computeError(svROM,svFOM,obj.normType,'abs');
        data = mean(tmp(~isnan(tmp)));
    case 'Rerr'
        svROM = modelobj.reconstructFullState(modelAuxobj);
        svFOM = modelFomobj.sv;
        data = computeError(svROM,svFOM,obj.normType,'rel');
    case 'MRerr'
        svROM = modelobj.reconstructFullState(modelAuxobj);
        svFOM = modelFomobj.sv;
        tmp = computeError(svROM,svFOM,obj.normType,'rel');
        data = max(tmp(~isnan(tmp)));
    case 'AvgRerr'
        svROM = modelobj.reconstructFullState(modelAuxobj);
        svFOM = modelFomobj.sv;
        tmp = computeError(svROM,svFOM,obj.normType,'rel');
        data = mean(tmp(~isnan(tmp)));    
    case 'linpt'
        if strcmpi(class(modelobj),'TPWL')
            data = modelobj.LinPt;
        else
            error('For data to be linpt, model must be TPWL');
        end
end
%If percent was requested for error, multiply by 100
if ~isempty(findstr(key,'err')) && percentFlag
    data = data*100;
end
    
end
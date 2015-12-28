function  [prop] = returnCSdata(txt)
%This function converts the property values specified in the char array txt
%into their appropriate MATLAB types
%--------------------------------------------------------------------------
%Inputs:
%-------
%txt  - a string containing the text from the current property value
%
%Outputs:
%--------
%prop - cell array of MATLAB types containing the property values in txt
%--------------------------------------------------------------------------

%If txt is empty, there is nothing to do, so just return empty
if isempty(txt)
    prop = []; return;
end
%Find the indices of all commas in the property value of interest
comind = findstr(txt,',');
%Store the number of commas located (n commas -> n+1 different values)
ncomma = length(comind);
%Find the indices of all left and right brackets in the property value of
%interest
lBrac = findstr(txt,'[');
rBrac = findstr(txt,']');

if ncomma == 0
    %If there are no commas, the property value is a "scalar" (in the sense
    %that only one value is specified), so just evaluate the "command" in
    %txt (command means any valid MATLAB statement)
    prop{1,1} = eval(txt);
elseif ncomma == 1
    %If there is 1 comma, the property value contains two "scalars", so
    %parse them property values appropriately, convert them from strings
    %to the appropriate MATLAB type, and store them in prop
    prop{1,1} = eval(txt(lBrac+1:comind-1));
    prop{1,2} = eval(txt(comind+1:rBrac-1));
else
    %If there are n commas, we have n "scalars", so parse them
    %apprpriately, convert them from strings to the appropriate MATLAB
    %type, and store them in prop
    
    %Initialize the output to be a cell array of appropriate size
    prop = cell(1,ncomma+1);
    
    %Convert 1st property "scalar" value into MATLAB type and store
    prop{1,1} = eval(txt(lBrac+1:comind(1,1)-1));
    for i = 2:ncomma
        %Convert 2nd through n-1st property "scalar" value into MATLAB type
        %and store 
        prop{1,i} = eval(txt(comind(1,i-1)+1:comind(1,i)-1));
    end
    %Convert last property "scalar" value into MATLAB type and store
    prop{1,end} = eval(txt(comind(1,end)+1:rBrac-1));
end
end
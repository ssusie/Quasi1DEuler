function  [out] = checkPostInput(txt,N,flag,curr_val)
%This code extracts the model properties when using the multiple input
%notation (i.e. in the ROM file).
%--------------------------------------------------------------------------
%Inputs:
%-------
%txt      - string containing the text from the input file
%N        - scalar that indicates the size of the array to be put in out
%          (i.e. this requires 'out' to be of size 1 or N)
%ind      - index of model array to extract 
%flag     - string declaring which property to extract
%curr_val - the current value of the property
%
%Outputs:
%--------
%out - value of the field at the ind_th row
%--------------------------------------------------------------------------

%Extract the line of code from txt that is going to be evaluted.
%Note that it should already be in MATLAB syntax.
tmp = determinePropText(txt,flag);
if isempty(tmp)
    %The property was not specified in the input file.  We need to use the
    %default!
    out = curr_val;
    return;
end

%Evaluate the MATLAB expression in tmp to get a (possibly) vector or cell
%array of entries.  We need to determine which row to use.
vec  = evalin('caller',tmp);
if isempty(N)
    %This indicates the user does not care how what the size of vec is.
    out = vec; return;
end

%Handle cases separately
if (size(vec,1) == N && size(vec,1) > 1) || (size(vec,1) == 1)
    out = vec;
elseif isempty(vec)
    if iscell(vec)
        out = {};
    else
        out = [];
    end
else
    error('Descriptive error message here!!!!');
end

end

% function  [out] = checkPostInput(txt,flag,N,default)
% tmp = determinePropText(txt,flag);
% if ~isempty(tmp)
%     out = evalin('caller',tmp);
%     if ~(size(out,1) == 1 || size(out,1) == N) && ~isempty(out)
%         error(['In the AXES field, ',flag, ' must be either have 1 entry or numPlotObj entries']);
%     end
% else
%     switch flag
%         case 'id'
%             error(['In the AXES field, ',flag, ' must be either have 1 entry or numPlotObj entries']);
%         otherwise
%             out = default;
%     end
% end
% 
% end
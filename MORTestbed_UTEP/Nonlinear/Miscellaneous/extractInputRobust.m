function  [out] = extractInputRobust(txt,flag,curr_val)
%This code extracts the model properties when using the multiple input
%notation (i.e. in the ROM file).
%--------------------------------------------------------------------------
%Inputs:
%-------
%txt      - string containing the text from the input file
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
    out = curr_val;
    return;
end

%Evaluate the MATLAB expression in tmp to get a (possibly) vector or cell
%array of entries.  We need to determine which row to use.
out  = evalin('caller',tmp);
end
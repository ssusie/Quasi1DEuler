function  [] = clearAllExcept(varargin)
%This function clears all variables in the workspace except for those
%specified in varargin
%--------------------------------------------------------------------------
%Inputs:
%-------
%varargin - a cell array of strings of the variables to KEEP in the
%           workspace.  All other variables will remain.
%
%Outputs:
%--------
%There are no outputs.
%--------------------------------------------------------------------------

s = evalin('base','who');
for i = 1:length(s)
    flag = false;
    for j = 1:length(varargin)
        if strcmpi(s{i},varargin{j})
            flag = true;
        end
    end
    if ~flag
        evalin('base',['clearvars(''',s{i},''')']);
    end
end
end
function [fname] = mapProbNum2DefaultFile(num)
%This function maps the problem number into the name of its default input
%file.
%--------------------------------------------------------------------------
%Inputs:
%-------
%num   - integer indicating the problem number 
%Outputs:
%--------
%fname - string containing the root of the default filename for this
%        problem
%--------------------------------------------------------------------------

switch num
    case 1
        fname = 'dBurger';
    case 2
        fname = 'dNLTrans';
    case 3
        fname = 'dFHN';
    case 4
        fname = 'dHNL2dSS';
    case 5
        fname = 'dConvDiffReact';
    case 6
        fname = 'dMEMS';
    case 7
        fname = 'dSteadyNozzle';
    case 8
        fname = 'dquasiEuler1D';
end
end
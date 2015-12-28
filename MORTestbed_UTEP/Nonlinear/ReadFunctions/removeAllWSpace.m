function  [txtNWS] = removeAllWSpace(txt)
%This function removes ALL white space (including \n, \t, etc) from txt.
%--------------------------------------------------------------------------
%Inputs:
%-------
%txt    - char array from which to removes all white spaces.
%
%Outputs:
%--------
%txtNWS - char array containing the txt char array without white spaces.
%--------------------------------------------------------------------------

ind = isstrprop(txt,'wspace');
txtNWS = txt(~ind);
end
function  [out] = round_nearest(num,ind)
%This function rounds num to the nearest decimal place indicated by ind
%--------------------------------------------------------------------------
%Inputs:
%-------
%num - any scalar to be rounded to a particular decimal place
%ind - number indicating decimal place you want to round to
%
%Outputs:
%--------
%out - num rounded to the decimal place indicated by indicator
%--------------------------------------------------------------------------

%Set up multiplication factor
factor = 1/ind;
%Multiply number by the factor, round to the closest whole number, then
%divide by factor
out = round(num*factor)/factor;

end
function [out] = ColumnwiseNorm(X,p)
%This function computes the nType-norm of each column of the the matrix X
%--------------------------------------------------------------------------
%Inputs:
%-------
%X   - n x k matrix
%p   - positive double (or inf) specifying the type of norm to use in
%      calculating the maximum relative error
%
%Outputs:
%--------
%out - 1 x k vector where each entry is the nType-norm of the corresponding
%      column of X 
%--------------------------------------------------------------------------

if nargin < 2, p=2; end;

switch p %Handle the two special p-norms separately
    case inf %Infinity norm (a.k.a. maximum norm, uniform norm, supremum norm)
        out = max(X,[],1); %Compute infinity norm of each column by max-ing
                           %over all rows in each column
    case 1 %1-norm (a.k.a. Taxicab norm or Manhattan norm)
        out = sum(abs(X),1); %Compute 1-norm by adding  the abs of all of the 
                             %entries in each column 
    otherwise %General p-norm
        out = (sum(abs(X).^p,1)).^(1/p);
end
end
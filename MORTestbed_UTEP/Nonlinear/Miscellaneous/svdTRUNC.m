function  [V] = svdTRUNC(X, tol)
%This function computes the svd of X and returns the left singluar vectors
%corresponding to singular values with magnitude greater than tol
%--------------------------------------------------------------------------
%Inputs:
%-------
%X   - matrix that we want to the compute the left singular values
%      containing the energy greater than tol
%tol - scalar specifying the minimum energy corresponding to left singular
%      values we will choose 
%Outputs:
%--------
%V   - matrix of left singular vectors with singular values greater than
%      tol
%--------------------------------------------------------------------------

%Compute SVD of X
[U,S,~] = svd(X);

%Find the index of the greatest singular value with magnitude less than tol
ind = find(diag(S) < tol,1);

%Truncate the left singular vectors with energy less than tol
V = U(:,1:ind-1);

end
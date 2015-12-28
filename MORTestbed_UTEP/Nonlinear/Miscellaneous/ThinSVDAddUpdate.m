function  [Un,Sn,Vn] = ThinSVDAddUpdate(U,S,V,A,B)
%This code will take the SVD of X = USV* and return the SVD of 
%Y = X + A^T*B = UnSnVn*.  U is N x r, S is r x r, and V is n x r, where r
%is the rank of X.
%--------------------------------------------------------------------------
%Inputs:
%-------
%U  - (N x ny) is the truncated left singular vectors
%S  - (ny x ny) are the singular values
%Vt1 - (ny x 1) is the column sum of the right singular vectors
%
%Outputs:
%--------
%Un  - updated (approx) left singular vectors (N x (ny+1))
%Sn  - updated (approx) singular values (ny+1 x ny+1)
%Vt1n - updated column sums of right singular vectors (nsnap x 1)
%--------------------------------------------------------------------------

%Rank of original SVD
r = size(S,1);

%Compute the quantities for the left update
m = U'*A;
M = A - U*m;
[P,Ra] = qr(M,0);
% Ra = P'*M;

n = V'*B;
N = B - V*n;
[Q,Rb] = qr(N,0);
% Rb = Q'*N;

%Form intermediate matrix
K = [m;Ra]*[n;Rb]';
K(1:r,1:r) = K(1:r,1:r) + diag(S,0);

%SVD of smaller matrix
[C,Sn,D] = svd(K,0);

%Update output quantities
Un = [U,P]*C;
Vn = [V,Q]*D;
end
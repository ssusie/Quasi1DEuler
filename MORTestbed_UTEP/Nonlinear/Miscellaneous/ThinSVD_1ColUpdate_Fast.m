function  [Un,Sn,Vn] = ThinSVD_1ColUpdate_Fast(U,S,V,c)
%This code will take the SVD of A = USV* and return the SVD of 
%B = [A,c] = UnSnVn*.
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

%Compute the quantities for the left update
m = U'*c;
p = c - U*m;
Ra = norm(p,2);
P = p/Ra;

%Rank of initial SVD
r = size(S,1);

%Form intermediate matrix
K = [diag(S,0),m;zeros(1,r),Ra];

%SVD of smaller matrix
[C,Sn,D] = svd(K,0);

%Update output quantities
Un = [U,P]*C;
% Vn = [V]*D;
Vn = [V*D(1:end-1,:);D(end,:)];
end
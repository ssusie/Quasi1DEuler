function  [Un] = ThinSVD_AddVec2Cols_Fast(U,S,Vt1,normq,a)
%This code will take the SVD of A = USV* and return the SVD of 
%B = A + a(1^T) = UnSnVn* (approximately).  The intended use is in the
%context of model order reduction.  Namely, we assume U (N x ny) is a
%reduced basis (truncation of left sv of snapshot matrix X), S (ny x ny) is
%are the singular values, and Vt1 (nsnap x 1) are the row sums of the right
%singular vectors (we don't need the whole V because in our SVD update to X
%+ ab^T, we know b = 1).
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
m = U'*a;
p = a - U*m;
Ra = norm(p,2);
P = p/Ra;

%Rank of initial SVD
ny = size(S,1);

%Form intermediate matrix
K = [m;Ra]*[Vt1',normq];

K(1:ny,1:ny) = K(1:ny,1:ny) + diag(S,0);

%SVD of smaller matrix
[C,~,~] = svd(K,0);

%Update output quantities
Un = [U,P]*C;

% sumQ = sum(Q);
% Vt1n = D'*[Vt1;sumQ];
% OneMinusVVt1n = OneMinusVVt1 + sum(Q)*Q;
end
function  [U,S,V] = LowRankSVDApprox(X,options)

%Extract the rank of the approximation
n = size(X,2);
k = min(options.approxRank,n);

%Generate Gaussian random test matrix and linear combination of the columns
%of X
Omega = randn(n,k);
Y = X*Omega;

%Apply power iteration for better convergence (if requested)
for q = 1:options.powerIters
    Y = X*(X'*Y);
end

%Extract orthonormal basis for Y 
[Q,~] = qr(Y,0);

%Compute low-rank SVD of X
B = Q'*X;
[Uhat,S,V] = svd(B,0);
U = Q*Uhat;

end
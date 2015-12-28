function  [Un,Sn] = ThinSVDUpdate_AppendVecs(U0,S0,X)

r = size(U0,2);

U0tX = U0'*X;

Pbar = X-U0*U0tX;

[P,R] = qr(Pbar,0);
ind = find(abs(diag(R)) < 1e-10,1);
if ~isempty(ind)
    P = P(:,1:ind-1);
end
RA = P'*Pbar;

c = size(RA,1);
K = [S0,U0tX;zeros(c,r),RA];
[Up,Sn] = svd(K,0);

Un = [U0,P]*Up;

end
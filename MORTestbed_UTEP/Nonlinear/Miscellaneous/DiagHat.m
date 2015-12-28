function Jtemp = DiagHat(x,M,N)
% Jtemp = diag(reshape(repmat(x',M,1),N*M,1));
Jtemp = spdiags(reshape(repmat(x',M,1),N*M,1),0,N*M,N*M);
end
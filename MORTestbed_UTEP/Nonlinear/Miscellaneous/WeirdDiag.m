function Jtemp2 = WeirdDiag(A,M,N)
Jtemp2 = reshape(repmat(A(:)',M,1),M*N,N);
end
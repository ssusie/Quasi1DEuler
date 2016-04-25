function residual = mach_eq(M,Mref,A,Aref,gamma) % for constant mass flow rate
%    residual = M^-2*((1+(gamma-1)/2*M^2)/(1+(gamma-1)/2))^((gamma+1)/(gamma-1)) - (A/Aref)^2;
    residual = (Mref/M)^2*((1+(gamma-1)/2*M^2)/(1+(gamma-1)/2*Mref^2))^((gamma+1)/(gamma-1)) - (A/Aref)^2;
%    residual = Mref^2*Aref^2*(1+(gamma-1)/2*M^2)^((gamma+1)/(gamma-1)) - M^2*A^2*(1+(gamma-1)/2*Mref^2)^((gamma+1)/(gamma-1));
%    residual = 0.5*(gamma+1)/(gamma-1)*(log(1+(gamma-1)/2*M^2)-log(1+(gamma-1)/2*Mref^2))+log(Mref)-log(M)-log(A)+log(Aref);
end

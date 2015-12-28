function residual = mach_eq(M,A,Athroat,gamma)
    residual = M^-2*(2/(gamma+1)*(1+(gamma-1)/2*M^2))^((gamma+1)/(gamma-1))-(A/Athroat)^2;
end
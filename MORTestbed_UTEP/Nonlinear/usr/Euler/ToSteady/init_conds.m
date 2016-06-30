function [M P rho u T] = init_conds(A, A_throat, Mach_ref, gamma, R, Tt, Pt, nVol, ShockLoc, PexitIncrease)
% function [M P rho u T] = init_conds(A_inlet, A_exit, A_throat, gamma, R, Tt, Pt, nVol)

A_inlet = A(1);
A_exit = A(end);

M_inlet = fzero(@(Mach) mach_eq(Mach,A_inlet,A_throat,gamma), .1);
M_exit = fzero(@(Mach) mach_eq(Mach,A_exit,A_throat,gamma), 2);

M = linspace(M_inlet,M_exit,nVol);

P   = Pt*(1+((gamma-1)/2)*M.^2).^(-gamma/(gamma-1));
T   = Tt./(1+(gamma-1)/2*M.^2);
rho = P./(R*T);
c   = sqrt(gamma*P./rho);
u   = M.*c;

deltaP = PexitIncrease*P(end);
P(end) = P(end) + deltaP;
rho(end) = rho(end) + deltaP/c(end)^2;
u(end) = rho(1)*A(1)*u(1)/(rho(end)*A(end));


end
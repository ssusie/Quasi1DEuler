function [M P rho u T] = init_conds(A_inlet, A_exit, A_throat, gamma, R, Tt, Pt, nVol)

M_inlet = fzero(@(Mach) mach_eq(Mach,A_inlet,A_throat,gamma), .1);
M_exit = fzero(@(Mach) mach_eq(Mach,A_exit,A_throat,gamma), 2);

M = linspace(M_inlet,M_exit,nVol);

P   = Pt*(1+((gamma-1)/2)*M.^2).^(-gamma/(gamma-1));
T   = Tt./(1+(gamma-1)/2*M.^2);
rho = P./(R*T);
c   = sqrt(gamma*P./rho);
u   = M.*c;

end
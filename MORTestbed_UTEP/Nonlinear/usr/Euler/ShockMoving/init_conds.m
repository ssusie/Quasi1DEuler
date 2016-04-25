function [M P rho u T] = init_conds(A, A_ref, Mach_ref, gamma, R, Tt, Pt, nVol, ShockLoc, PexitIncrease)

M = zeros(1,nVol);

for i=1:ShockLoc-1
    M(i) = fzero(@(Mach) mach_eq(Mach,Mach_ref,A(i),A_ref,gamma), 2.6);
    M0 = 0.1; N = 1;
    while ( (M(i) < 1 || isnan(M(i))) && N < 10)
        M(i) = fzero(@(Mach) mach_eq(Mach,Mach_ref,A(i),A_ref,gamma), 1+M0);
        M0 = M0 / 10;    N = N + 1;
    end
end

P   = Pt*(1+((gamma-1)/2)*M.^2).^(-gamma/(gamma-1));
T   = Tt./(1+(gamma-1)/2*M.^2);
rho = P./(R*T);
c   = sqrt(gamma*P./rho);
u   = M.*c;

e   = 0.5*rho(ShockLoc-1)*u(ShockLoc-1)^2 + P(ShockLoc-1)/(gamma-1);
m = rho(ShockLoc-1)*u(ShockLoc-1);
n = rho(ShockLoc-1)*u(ShockLoc-1)^2 + P(ShockLoc-1);


h = (e+P(ShockLoc-1))/rho(ShockLoc-1);

a = 0.5-gamma/(gamma-1);
b = gamma/(gamma-1)*n/m;
c = -h;

u(ShockLoc) = min((-b+sqrt(b^2-4*a*c))/(2*a),(-b-sqrt(b^2-4*a*c))/(2*a));
rho(ShockLoc) = m/u(ShockLoc);
P(ShockLoc) = n-rho(ShockLoc)*u(ShockLoc)^2;
M(ShockLoc) = u(ShockLoc)/sqrt(gamma*P(ShockLoc)/rho(ShockLoc));

for i=ShockLoc+1:nVol
    M(i) = fzero(@(Mach) mach_eq(Mach,M(ShockLoc),A(i),A(ShockLoc),gamma), 0.1);
%     M0 = 0.1; N = 1;
%     while ( (M(i) < 1 || isnan(M(i))) && N < 10)
%         M(i) = fzero(@(Mach) mach_eq(Mach,M(ShockLoc),A(i),A(ShockLoc),gamma), 0.01);
%         M0 = M0 / 10;    N = N + 1;
%     end
end

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

function [alpha_star,success] = linesrchBackNewton(phi)
% linesearch algorithm with Lagrange interpolation
% input: phi is a function that returns both the value of the 1D
% directional function phi(alpha)=f(x+alpha*p) and its first derivative
% LS_param contains parameters for the Line search (optional)
% output: alpha_star is a step length that satisfies sufficient
% decrease condition

%Extract parameters
rho = 0.5;
c = 0.15;
iterMax=50;
alpha_max=1;

%Linesearch
[phi0,dphi0] = phi(0);
% phi0=phi(0);
% dphi0=0;

if (dphi0 > 0)
    disp('Warning: this is not a descent direction');
    sign=-1;
else
    sign = 1;
end

%Start with Newton step if allowed
alpha_new =  sign*min(1,alpha_max);

%Backtracking linesearch
iter = 0;
while (iter<iterMax)
    phi_new = phi(alpha_new);
    
    if ( phi_new <= phi0 + c*alpha_new*dphi0 )
        alpha_star = alpha_new;
        success = 1;
        return;
    else
        alpha_new = rho*alpha_new;
        iter = iter + 1;
    end
end

disp('Sufficient Descent condition is not satisfied !');

if phi_new > phi0
    disp('--- Phi did not decrease...Setting set length to 0 ---');
    alpha_star=0;
    success=0;
    return;
end

success = 0;
alpha_star = alpha_new;
end
function [alpha_star,success,nFunCalls,alphaCalls] = linesearchBacktracking(phi_initial,LSparam)
% linesearch algorithm with Lagrange interpolation
% input: phi is a function that returns both the value of the 1D
% directional function phi(alpha)=f(x+alpha*p) and its first derivative
% LS_param contains parameters for the Line search (optional)
% output: alpha_star is a step length that satisfies sufficient
% decrease condition

%Extract parameters
c1 = LSparam(1);
c2 = LSparam(2);
alpha_max = LSparam(3);
iterMax = LSparam(4);
maxZoom = LSparam(5);
tolInterp = LSparam(6);
eps = LSparam(7);
extrapolBeta = LSparam(8);

%Linesearch
[phi0,dphi0] = phi_initial(0);
nFunCalls = 1;
alphaCalls = 0;
if (dphi0 > 0)
    disp('Warning: this is not a descent direction');
    sign=-1;
    alpha_star = 0;
    success = 0;
%     return;
else
    sign = 1;
end
phi = @(x) phi_initial(sign*x);

%Start with Newton step if allowed
alpha_new =  min(1,alpha_max);

%Backtracking linesearch
iter = 0;
success = 0;
while (iter<iterMax)
    [phi_new,~] = phi(alpha_new);
    nFunCalls = nFunCalls + 1;
    alphaCalls = [alphaCalls,alpha_new];
    if ( phi_new <= phi0 + c1*alpha_new*dphi0 )
        alpha_star = alpha_new;
        success = 1;
        return;
    else
        alpha_new = alpha_new/2;
        iter = iter + 1;
    end
end
% disp('Sufficient Descent condition is not satisfied !');
success = 0;
alpha_star = sign*alpha_new;
end
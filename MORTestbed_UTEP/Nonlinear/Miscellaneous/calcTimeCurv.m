function  [rho] = calcTimeCurv(x,dt)
%This function computes the pointwise curvature of the curve defined by x.
%--------------------------------------------------------------------------
%Inputs:
%-------
%x   - N x linptMAX matrix where the ith column is a vector in R^N at time
%      step i.  All columns define a curve in R^N parameterized by time t.
%dt  - scalar specifying time step size
%
%Outputs:
%--------
%rho - linptMAX x 1 vector defining the pointwise curvature (curvature of
%      curve defined by x at each time step).  Curvature not defined at end
%      points so they are set to zero.
%--------------------------------------------------------------------------

%Select a norm to use
nNum = 2;

%Determine the number of time steps
n = size(x,2);
%Initialize curvature vector
rho = zeros(n,1);

%Loop over each non-boundary point of the curve
for i = 2:n-1
    %Compute the first derivative of curve wrt t via finite differences
    gamma_p  = (x(:,i+1) - x(:,i-1))/(2*dt);
    %Compute the second derivative of curve wrt t via finite differences
    gamma_pp = (x(:,i+1) - 2*x(:,i) + x(:,i-1))/(dt^2);
    %Compute the norm of the first derivative of the curve wrt t
    n_gamma_p = norm(gamma_p,nNum);
    %Compute the curvature at each interior point of the curve
    rho(i,1) = n_gamma_p^3/sqrt(n_gamma_p^2*norm(gamma_pp,nNum)^2 - (gamma_p'*gamma_pp)^2);
end
end
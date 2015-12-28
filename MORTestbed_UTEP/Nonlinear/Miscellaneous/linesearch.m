function [alpha_star,success,nFunCalls,alphaCalls] = linesearch(phi_init,LSparam)
%This function performs a linesearch with Lagrange interpolation.  The
%algorithm returns a step length tha satisfies the strong Wolfe conditions.
%Numerical Optimization - Nocedal and Wright, pg. 60-61.
%--------------------------------------------------------------------------
%Input:
%------
%obj        - model object
%phi        - function handle that returns the 1D directional funtion
%             phi(alpha) = f(x+alpha*p)
%dphi       - function handle that returns the first derivative of the 1D
%             directional function phi(alpha) = f(x+alpha*p)
%Output:
%-------
%alpha_star - step length that satisfies the strong Wolfe conditions.  Upon
%             exit, the model has been evaluated at alpha_star.
%--------------------------------------------------------------------------

%Extract parameters
c1 = LSparam(1);
c2 = LSparam(2);
alpha_max = LSparam(3);
iterMax = LSparam(4);
maxZoom = LSparam(5);
tolInterp = LSparam(6);
eps = LSparam(7);
extrapolBeta = LSparam(8);
      
%Determine local function properties at current iterate
alpha_old    = 0;
[phi0,dphi0] = phi_init(0);
% nFunCalls = 1;
% alphaCalls = 0;
nFunCalls  = 0;
alphaCalls = [];
%Check for decent direction.  Adjust if necessary.
if (dphi0>0)
    disp('Warning: this is not a descent direction: switching direction');
    dphi0 = -dphi0;
    phi = @(x) phi_init(-x);
    sign = -1;
else
    phi = phi_init;
    sign = 1;
end
      
% zoom_loc = @(alpha_lo_loc,alpha_hi_loc) ...
%     zoom(alpha_lo_loc,alpha_hi_loc,phi0,dphi0,c1,c2,maxZoom,tolInterp,eps,phi);
zoom_loc = @(alpha_lo_loc,alpha_hi_loc,phi_lo_loc) ...
                 zoom(alpha_lo_loc,alpha_hi_loc,phi_lo_loc,phi0,dphi0,c1,c2,maxZoom,tolInterp,eps,phi);
phi_old = phi0;
%dphi_old = dphi0;

%Start with Newton step, if allowed
% alpha_new = min(1,0.5*alpha_max);
alpha_max=min(alpha_max,1);
alpha_new = 0.8*alpha_max;

%Iterate to find appropriate step length
iter = 0;      
success = 0;            
       
% %Check if maximum step satisfies the strong Wolfe conditions
% [phimax,dphimax]=phi(alpha_max);
% if (phimax <= phi0 + c1*alpha_max*dphi0) && (abs(dphimax)<=c2*abs(dphi0))
%     alpha_star=sign*alpha_max;
%     success=1;
%     nFunCalls=1;
%     alphaCalls = sign*alpha_max;
%     return;
% end


while (iter<iterMax)
    [phi_new,dphi_new] = phi(alpha_new);
    nFunCalls = nFunCalls + 1;
    alphaCalls = [alphaCalls,alpha_new];
    if (phi_new > phi0 + c1*alpha_new*dphi0 || (phi_new >= phi_old && iter>0))
        [alpha,success,nCalls,aCalls]  = zoom_loc(alpha_old,alpha_new,phi_old);
%         [alpha,success]  = zoom(alpha_old,alpha_new,phi0,dphi0,c1,c2,maxZoom,tolInterp,eps,phi);
        nFunCalls = nFunCalls + nCalls;
        alphaCalls = [alphaCalls,aCalls]; %#ok<*AGROW>
        if (success)
            alpha_star = sign*alpha;
            return;
        end
    end
    
    if (abs(dphi_new) <= -c2*dphi0)
        alpha_star = sign*alpha_new;
        success = 1;
        return;
    end
    
    if (dphi_new >= 0)
        [alpha,success,nCalls,aCalls]  = zoom_loc(alpha_new,alpha_old,phi_new);
%         [alpha,success]  = zoom(alpha_new,alpha_old,phi0,dphi0,c1,c2,maxZoom,tolInterp,eps,phi);
        nFunCalls = nFunCalls + nCalls;
        alphaCalls = [alphaCalls,aCalls];
        if (success)
            alpha_star = sign*alpha;
            return;
        end
    end
    alpha_old = alpha_new;
    phi_old = phi_new;
    %   dphi_old = dphi_new;
    % alpha_new = extrapol(alpha_old,alpha_max,extrapolBeta);
    alpha_new = extrapol2(alpha_old,alpha_max);
    iter = iter+1;
end

if alpha_max < 1
    disp('Using alpha_max!');
else
    disp('Line search Wolfe conditions not satisfied !');
end
success = 0;
% alpha_star=sign*alpha_old;
alpha_star = sign*alpha_max;
phi(alpha_star);
end

function [alpha_star,success,nFunCalls,aCalls] = zoom(alpha_lo,alpha_hi,phi_lo,phi0,dphi0,c1,c2,maxZoom,tollInterp,eps,phi)

%zoom in interval [alpha_lo,alpha_hi]

% interpol_loc = @(alpha_lo_loc,alpha_hi_loc)...
%            interpol(alpha_lo_loc,alpha_hi_loc,phi,tollInterp,eps);
interpol_loc = @(alpha_lo_loc,alpha_hi_loc)...
    interpol2(alpha_lo_loc,alpha_hi_loc);
iterMax = maxZoom; %max of zoom steps
iter = 0;
success = 1;
nFunCalls = 0;
aCalls = [];
%  disp('zoom');

while (iter<iterMax)
    
    alpha_new = interpol_loc(alpha_lo,alpha_hi);
    [phi_new,dphi_new] = phi(alpha_new);
%     [phi_lo,~] = phi(alpha_lo);
    
    nFunCalls = nFunCalls + 1;
    aCalls = [aCalls,alpha_new];
    
    if (phi_new > phi0 + c1*alpha_new*dphi0 || phi_new >= phi_lo)
        alpha_hi = alpha_new;
    else
        if (abs(dphi_new) <= -c2*dphi0)
            alpha_star = alpha_new;
            return;
        end
        if (dphi_new*(alpha_hi-alpha_lo) >= 0 )
            alpha_hi = alpha_lo;
        end
        alpha_lo = alpha_new;
        phi_lo = phi_new;
    end
    iter = iter+1;
end

disp('Leaving zoom : max iterations were exceeded');
success = 0;
alpha_star = alpha_new;
end

function alpha_new = interpol2(alpha_lo,alpha_hi)

alpha_new = (alpha_hi+alpha_lo)/2;

end

function alpha_new = extrapol2(alpha_lo,alpha_hi)
if (alpha_lo>alpha_hi)
    alpha_new = alpha_hi;
else
    alpha_new = interpol2(alpha_lo,alpha_hi);
end
end
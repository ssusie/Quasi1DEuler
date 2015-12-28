function [B] = updateHessian(Bk,y,s,QNflag,tol)
%This function computes the inverse hessian update using the
%specified QuasiNewton update
%--------------------------------------------------------------------------
%Inputs:
%-------
%Hp     - previous Hessian inverse approximation
%y      - step taken in parameter space, i.e. alpha*dp
%s      - gradient difference
%QNflag - string indicating which QN inverse hessian update to use
%
%Outputs:
%--------
%H      - updated Hessian inverse approximation
%--------------------------------------------------------------------------

switch QNflag
    case 'BFGS'
        %BFGS - Broyden
        %  - rank 2 update
        %  - symmetric
        %  - guaranteed SPD if H0 SPD
        %  - may eventually suffer from ill-conditioning
        
        B = Bk;
        rho = y'*s;
        if (rho>tol)
            v = Bk*s;
            B = Bk - v*v'/(s'*v) + (1/rho)*(y*y');
        end
    case 'DFP'
        %DFP - Davidson, Fletcher, Powell update
        %  - rank 2 update
        %  - symmetric
        %  - guaranteed SPD if H0 SPD
        %  - may eventually suffer from ill-conditioning
        
    case 'SR1'
        %Symmetric Rank 1 update
        %  - rank 1 update
        %  - symmetric
        %  - not guaranteed SPD (even if H0 SPD)
        %  - may eventually suffer from ill-conditioning
        
end
end
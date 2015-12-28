function [H] = updateInverseHessian(Hp,y,s,QNflag,tol)
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
        
        H = Hp;
        rho = y'*s;
%         if(rho>1e-3)
%         if (rho>1e-12)
        if (rho>tol)
            M = speye(size(Hp)) - 1/rho*s*y';
            H = M*Hp*M' + 1/rho*(s*s');
        end
        
%         rho = 1/(y'*s);
%         if (abs(rho)>1e-3)
%             M = speye(size(Hp)) - rho*s*y';
%             H = M*Hp*M' + rho*(s*s');
%         end
    case 'DFP'
        %DFP - Davidson, Fletcher, Powell update
        %  - rank 2 update
        %  - symmetric
        %  - guaranteed SPD if H0 SPD
        %  - may eventually suffer from ill-conditioning
        
        Hy = Hp*y;
        rho = 1/(y'*s);
        gamma = y'*Hy;
        if (abs(rho)>1e-3) && (abs(gamma)>1e-3)
            H = Hp - (Hy/gamma)*Hy' + (s/rho)*s';
        end
    case 'SR1'
        %Symmetric Rank 1 update
        %  - rank 1 update
        %  - symmetric
        %  - not guaranteed SPD (even if H0 SPD)
        %  - may eventually suffer from ill-conditioning
        
        Hy = Hp*y;
        v = (s - Hy);
        rho = v'*y;
        if (abs(rho)>1e-16)
            H = Hp + (v/rho)*v';
        end
end
end

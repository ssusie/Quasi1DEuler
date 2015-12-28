function  [V,W,ndofr] = krylov(A,b,c,k,flag)
%This function computes bases V and W of the krylov subspaces of the SISO
%linear, dynamical system defined by A, b, and c.  
%--------------------------------------------------------------------------
%Inputs:
%-------
%A     - nxn dynamics matrix of the system
%b     - b is the nx1 input matrix of the system
%c     - c is the 1xn output matrix of the system
%k     - integer specifying the order of the reduced system
%flag  - integer specifying the algorithm to use for constructing the
%        Krylov subspaces (1 = two-sided Lanczos, 2 = Arnoldi)
%
%Outputs:
%--------
%V     - right Krylov subspace of the dynamical system defined by A, b, c
%W     - left Krylov subspace of the dynamical system defined by A, b, c
%ndofr - integer indicating the size of the reduced order model (may be
%        less than k due to the rank of V and W).
%--------------------------------------------------------------------------

n = size(A,1);
V = zeros(n,k);
switch flag
    case 1 %Two-sided Lanczos Algorithm (algorithm from Antoulas pg. 334)
        W = zeros(n,k);
        beta = sqrt(abs(b'*c'));
        gamma = sign(b'*c')*beta;
        V(:,1) = b/beta;
        W(:,1) = c'/gamma;
        
        for j = 1:k
            AV = A*V(:,j);
            alpha = W(:,j)'*AV;
            
            if j > 1
                r = AV - alpha*V(:,j) - gamma*V(:,j-1);
                q = A'*W(:,j) - alpha*W(:,j) - beta*W(:,j-1);
            else
                r = AV - alpha*V(:,j);
                q = A'*W(:,j) - alpha*W(:,j);
            end
            
            beta = sqrt(abs(r'*q));
            gamma = sign(r'*q)*beta;
            
            if j < k
                V(:,j+1) = r/beta;
                W(:,j+1) = q/gamma;
            end
        end
        
        R = min(rank(V),rank(W));
            
        V = V(:,1:R);
        W = W(:,1:R);
    case 2 %Arnoldi Algorithm (algorithm from Antoulas pg. 335 with modified Gram-Schmidt)
        if norm(b) == 0
            V(1,1) = 1;
        else
            V(:,1) = b/norm(b);
        end
        w = A*V(:,1);
        alpha = V(:,1)'*w;
        f = w-V(:,1)*alpha;
        
        for j = 1:k-1
            beta = norm(f);
            V(:,j+1) = f/beta;
            f = A*V(:,j+1);
            for l = 1:j+1   % modified Gram-Schmidt
                h = V(:,l)'*f;
                f = f - V(:,l)*h;
            end
        end
        R = rank(V);
        V  = V(:,1:R);
        W = V;
    otherwise
        error('Only Two-Sided Lanczos (1) and Arnoldi (2) Methods are supported');
end
ndofr = min([R,k]);
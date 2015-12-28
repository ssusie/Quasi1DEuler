function U=iterative_svd2(A,k)
    A1=A(:,1:k);
    %% Step 1: URV decomposition (using SVD here)
    [U,R,V]=svd(A1,0);

    for i=k:size(A,2)-1
        %% Step 2: Project b
        b=A(:,i+1);
        alpha=U'*b;
        u=(eye(size(U,1))-U*U')*b;
        alpha1=norm(u);
        u1=u./alpha1;
        
        %% Step 3: Update U, R and V
        U1=[U u1];
        
        R1=[R alpha];
        R1(length(R1),length(R1))=alpha1;
        
        % Must to have the following re-orthogonalize step
        % Think of an efficient way of doing this (do once in a while)
        [U1,temp]=qr(U1,0); 
        R1=temp*R1;
        
        V1=V'; % Transposed to follow convention in paper A=URV
        V1(length(R1),length(R1))=1;
        
        %% Step 4: Reduce R1 to Upper tridiagonal
        [Q1g,R1g,Q2g]=givens_rotations1(R1);
        [Q3g,R2g,Q4g]=givens_rotations2(R1g);
        
        % get the householder matrix to update URV to USV (svd)   
        [R1u,R1r,R1v]=svd(R2g,0);
        x=R1u(:,length(R1));
        x=Q1g'*Q3g'*x;
                
        [v,beta,s] = gallery('house',x);
        Gu = eye(length(R1)) - beta*v*v';
        Gu=flipud(Gu);
        Gu=Gu./s;
        S1=Gu*R1;
        %% Step 4: Compute RQ factorization of S1
        [R2,H]=rq(S1);
        %% Step 5: Update the matrices
        U=U1*Gu';
        V=H*V1;
        R=R2(1:k,1:k);
        U=U(:,1:k);
        V=V(1:k,1:k);
    end
    [Ru,Rr,Rv]=svd(R);
    U=U*Ru;
      
    %% RQ factorization of a matrix
    function [R,Q]=rq(B)
       [Q R]=qr(flipud(B).');
       R=flipud(R.');
       R=fliplr(R);
       Q=flipud(Q.');  
    end
end
function U=iterative_svd(A,k)
tic
    A1=A(:,1:k);
    %% Step 1: URV decomposition (using SVD here)
    [U,R,V]=svd(A1,0);

    for i=k:size(A,2)-1
        %% Step 2: Update via modified Gram-Schmidt
        b=A(:,i+1);
        r=U'*b;
        b1=b-U*r;
        row=norm(b1);
        u1=b1/row;
        
        %% Step 3: get last left singular vector of R1
        R1=[R r];
        R1(length(R1),length(R1))=row;
        
        U1=[U u1];
        
        % Must to have the following re-orthogonalize step
        % Think of an efficient way of doing this (do once in a while)
        [U1,temp]=qr(U1,0); 
        R1=temp*R1;

        % get the householder matrix to update URV to USV (svd)   
        [R1u,R1r,R1v]=svd(R1,0);
        x=R1u(:,length(R1));
        
        [v,beta,s] = gallery('house',x);
        Gu = eye(length(R1)) - beta*v*v';
        Gu=flipud(Gu);
        Gu=Gu./s;
        S1=Gu*R1;

        %% Step 4: Compute RQ factorization of S1
        [R2]=rq1(S1);

        %% Step 5: Update the matrices
        U=U1*Gu';
        R=R2(1:k,1:k);
        U=U(:,1:k);
    end
    [Ru,Rr,Rv]=svd(R);
    U=U*Ru;
toc    

    %% RQ factorization of a matrix just to get R
    function [R]=rq1(B)
       [Q R]=qr(flipud(B).');
       R=flipud(R.');
       R=fliplr(R); 
    end

    %% RQ factorization of a matrix
    function [R,Q]=rq(B)
       [Q R]=qr(flipud(B).');
       R=flipud(R.');
       R=fliplr(R);
       Q=flipud(Q.');  
    end
end
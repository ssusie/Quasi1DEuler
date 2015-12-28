function  [nnzeros,irow_crs,irow_coor,jcol_coor] = createSparsityStructure(msh,flag)
%This assumes that the matrix (jacobian) does not have an all zero row.

%Assumes LM is zero based.  Hence the adding of ones everywhere.

if nargin == 1 || isempty(flag)
    flag=1;
end

switch flag
    case 1
        %This algorithm is far more memory efficient than algorithm 2 and
        %should be used on large problems.  May not be as fast as algorithm
        %2 for smaller problems.
        irow_coor=zeros(double(msh.nel)*double(msh.nsd)*double(msh.nen),1);
        jcol_coor=zeros(double(msh.nel)*double(msh.nsd)*double(msh.nen),1);
        cnt=1;
        for e = 1:double(msh.nel)
            LMpos = msh.LM(:,e)>-1;
            ncol = sum(LMpos);
            
            irow_coor(cnt:cnt+ncol^2-1) = reshape(repmat(msh.LM(LMpos,e)',ncol,1),ncol^2,1);
            jcol_coor(cnt:cnt+ncol^2-1) = repmat(msh.LM(LMpos,e),ncol,1);
            cnt=cnt+ncol*ncol;
        end
        C=unique([irow_coor,jcol_coor],'rows');
        irow_coor=C(:,1)+1; jcol_coor=C(:,2)+1; clear C; %1-based indexing
        nnzeros=length(irow_coor);
        
        %Compressed row format from coordinate format
        irow_crs=[0;find(diff(irow_coor));nnzeros+1];
        
        %zero-based indexing (because this will be used by C++)
        %irow_coor, jcol_coor are 1-based indexing (because they are used most by
        %MATLAB)
        irow_crs = int32(irow_crs);
        irow_coor=int32(irow_coor);
        jcol_coor=int32(jcol_coor);
    case 2
        %This algorithm is very memory inefficient but relatively fast for
        %small problems.  Not feasible for large meshes.
        df=false(msh.ndof,msh.ndof);
        for e = 1:msh.nel
            LMpos = msh.LM(:,e)>-1;
            
            df(msh.LM(LMpos,e)+1,msh.LM(LMpos,e)+1)=true;
        end
        [jcol_coor,irow_coor] = find(df');
        nnzeros=length(irow_coor);
        
        irow_crs=zeros(msh.ndof+1,1);
        irow_crs(1) = 0;
        for i = 2:msh.ndof
            irow_crs(i) = irow_crs(i-1)+nnz(df(i-1,:));
        end
        irow_crs(end) = nnzeros+1;
        
        %zero-based indexing (because this will be used by C++)
        %irow_coor, jcol_coor are 1-based indexing (because they are used most by
        %MATLAB)
        irow_crs = int32(irow_crs);
        irow_coor=int32(irow_coor);
        jcol_coor=int32(jcol_coor);
end
end
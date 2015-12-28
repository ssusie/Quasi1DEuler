function  [msh] = buildMatrices(msh,dbc,t)
%Build matrices for mapping between local and global quantities

%ID array: map from dof# and Gnode# to Gdof#
msh.ID=zeros(msh.nsd,msh.np);
dofcnt=0;
dbccnt=0;
for j=1:double(msh.np)
    for i=1:double(msh.nsd)
        if dbc.loc(i,j)
            dbccnt=dbccnt-1;
            msh.ID(i,j)=int32(dbccnt);
        else
            dofcnt=dofcnt+1;
            msh.ID(i,j)=int32(dofcnt-1);
        end
        
    end
end
msh.ID = int32(msh.ID);
%IEN array: map from Lnode# and elem# to Gnode#
msh.IEN=uint32(t)';

%LM array: map from Lnode# and elem# to dof#
msh.LM=zeros(double(msh.nsd)*double(msh.nen),double(msh.nel));
for i=1:double(msh.nsd)
    for j=1:double(msh.nen)
        for e=1:double(msh.nel)
            msh.LM(i+double(msh.nsd)*(j-1),e)=msh.ID(i,msh.IEN(j,e));
        end
    end
end
msh.IEN=uint32(msh.IEN-1);
msh.LM=int32(msh.LM);
end
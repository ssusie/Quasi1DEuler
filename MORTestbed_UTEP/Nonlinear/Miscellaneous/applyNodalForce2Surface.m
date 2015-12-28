function  f=applyNodalForce2Surface(ts,dir,val,ID,ndof)

surfnodes=unique(ts);

f=zeros(ndof,1);
for i = 1:length(dir)
    id = 1+ID(dir(i),surfnodes)';
    id(id<0,:)=[];
    
    f(id)=f(id)+val;
end

end
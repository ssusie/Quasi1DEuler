function   [U,S] = incrementalUpdateSVD(X,options)

m = size(X,2);
numPartitions = ceil(m/options.maxvecs);


[U,S]=svd(X(:,1:min(options.maxvecs,m)),0);
truncSize=min(options.truncIntermediate,size(U,2));
Ut=U(:,1:truncSize);
St=S(1:truncSize,1:truncSize);

for i = 2:numPartitions
    
    [U,S] = ThinSVDUpdate_AppendVecs(Ut,St,X(:,(i-1)*options.maxvecs+1:min(m,i*options.maxvecs)));
    truncSize=min(size(Ut,2)+options.truncIntermediate,size(U,2));
    Ut=U(:,1:truncSize);
    St=S(1:truncSize,1:truncSize);
    
end

end
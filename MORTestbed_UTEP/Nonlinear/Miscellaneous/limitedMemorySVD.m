function   [U,S] = limitedMemorySVD(X,options)
%This function computes the limited memory of X according to the options in
%the options structure.
%--------------------------------------------------------------------------
%Inputs:
%-------
%X       - an m x n matrix
%options - structure containing the options for the algorithm.  Fields are:
%          algorithm, approxRank, maxvecs, trucIntermediate,
%          includeSingVals, and truncFinal
%
%Outputs:
%--------
%--------------------------------------------------------------------------

%Define number of partitions of the matrix X based on the maximum number of
%vectors we are allowed to take in an SVD
n = size(X,2);
numPartitions = ceil(n/options.maxvecs);

Uconcat=[];
for i = 1:numPartitions
    %For each partition, compute the SVD and truncate to the intermediate
    %truncation value
    [Uf,Sf] = svd(X(:,(i-1)*options.maxvecs+1:min(i*options.maxvecs,n)),0);
    truncSize = min(options.truncIntermediate,size(Uf,2));
    
    %Concatenate remaining singular vectors in a matrix
    if options.includeSingVals
        Uconcat = [Uconcat,Uf(:,1:truncSize)*Sf(1:truncSize,1:truncSize)];
    else
        Uconcat = [Uconcat,Uf(:,1:truncSize)];
    end
end
%Compute svd of matrix of truncated SVDs to get an approximation of the
%true SVD
[U,S] = svd(Uconcat,0);

end
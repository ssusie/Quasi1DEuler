function  [snap] = snapshotCollection(X,coll,ic,Y,z)
%This function determines a snapshot collection from a set of state
%vectors.
%--------------------------------------------------------------------------
%Inputs:
%-------
%X    - state vector matrix with nstep+1 columns (first column is ic)
%coll - string specifying which state vector snapshot collection to use
%ic   - vector defining the initial condition that is to be used if coll =
%       'ic'
%Y    - state vector matrix of "previous states" in the event that the
%       columns of X do not correspond to consecutive snapshots
%z    - some state vector to reference from which to reference the snapshots.
%Outputs:
%--------
%snap - snapshot matrix
%--------------------------------------------------------------------------

switch coll
    case 'zero'
        %If we are using unreference state vectors as snapshots, use all
        %state vectors except initial condition
        snap = X;
    case 'ic'
        %If we are using state vectors that reference the initial condition
        %as snapshots, use all state vectors except initial condition and
        %subtract the initial condition from each one
        snap = bsxfun(@minus,X,ic); 
        %This is just a faster command to accomplish:
        %X - repmat(ic,size(X,2))
    case 'prev'
        %If we are using state vectors that reference the previous state
        %as snapshots, use all state vectors except initial condition and
        %subtract the previous one
        snap = X-Y;
    case {'somevec','clustcenters'}
        %We use state vectors that reference the cluster center
        snap = bsxfun(@minus,X,z);
end

end
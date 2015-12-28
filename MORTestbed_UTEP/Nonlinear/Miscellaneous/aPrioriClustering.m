function [closClusIndex,clustCenters,additionalPoints] = aPrioriClustering(nClusters,pointsCoord,addElemThresh)
%This function computes the a-priori clustering for  a set of points.
%--------------------------------------------------------------------------
%Inputs:
%-------
%nClusters          - initial estimate for the number of clusters
%pointsCoord        - coordinates of the points
%additionalElements - boolean indicating which whether or not to add
%                     additional elements to each cluster
%addElemThresh      - scalar indicating the percentage of elements to add
%                     from neighboring clusters
%
%Outputs:
%--------
%closClusIndex      - index of the closest cluster center to a given point
%clustCenters       - coordinates of the clusters centers
%additionalPoints   - vector containing indices of additional snapshots to
%                     add to each local snapshot matrix
%--------------------------------------------------------------------------
%Copyright David Amsallem 2013.

%Determine whether to allow sharing

addElem = (addElemThresh>0);

nPoints = size(pointsCoord,2);

closClusIndex = zeros(nPoints,1);
clustCenters = zeros(size(pointsCoord,1),nClusters);

nPointsPerCluster = floor(nPoints/nClusters);
nPointsInLastCluster = nPoints - (nClusters-1)*floor(nPoints/nClusters);

for iCluster = 1:nClusters-1
    block = (iCluster-1)*nPointsPerCluster+1:iCluster*nPointsPerCluster;
    closClusIndex(block,1) = iCluster*ones(nPointsPerCluster,1); 
    clustCenters(:,iCluster) = 1/nPointsPerCluster*sum(pointsCoord(:,block),2);
end
block = (nClusters-1)*nPointsPerCluster+1:nPoints;
closClusIndex(block,1) = nClusters*ones(nPointsInLastCluster,1); 
clustCenters(:,nClusters) = 1/nPointsInLastCluster*sum(pointsCoord(:,block),2);


% add closest elements to each cluster if required
if addElem   
    % find how many elements should  be added to each cluster
    additionalPoints(nClusters).Indices = 0;
    for iCluster=1:nClusters
        addPointsInCluster = ceil(addElemThresh*nPointsPerCluster);
        
        switch iCluster
            case 1
                additionalPoints(iCluster).Indices = (nPointsPerCluster+1:nPointsPerCluster+addPointsInCluster);
            case nClusters
                additionalPoints(iCluster).Indices = (nPoints-nPointsInLastCluster-addPointsInCluster+1:nPoints-nPointsInLastCluster);
            otherwise
                addPointsInClusterBefore = floor(addPointsInCluster/2);
                addPointsInClusterAfter =  addPointsInCluster - addPointsInClusterBefore;
                additionalPoints(iCluster).Indices = [(iCluster-1)*nPointsPerCluster-addPointsInClusterBefore+1:(iCluster-1)*nPointsPerCluster,iCluster*nPointsPerCluster+1:iCluster*nPointsPerCluster+addPointsInClusterAfter];
        end
    end
    
else
    additionalPoints = 0;
end
end

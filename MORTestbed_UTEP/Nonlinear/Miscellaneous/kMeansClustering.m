function [closClusIndex,clustCenters,additionalPoints] = kMeansClustering(nClusters,pointsCoord,addElemThresh)
%This function computes the k-means clustering for  a set of points.
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
%Copyright David Amsallem 2010. Modifications by Matthew Zahr 2011 for
%integration with MORTestbed.

%Determine whether to allow sharing
addElem = (addElemThresh>0);

stream = RandStream('mrg32k3a','seed',71489);
nPoints = size(pointsCoord,2);

closClusIndex = zeros(nPoints,1);
if addElemThresh>0
    secondClosClusIndex = zeros(nPoints,1);
    secondDist = zeros(nPoints,1);
else
    secondClosClusIndex = 0;
    secondDist = 0;
end
X = pointsCoord;

% random initial centers
[~,Irand] = sort(rand(stream,1,nPoints));
Irand = Irand(1:nClusters);
clustCenters = X(:,Irand);
clustCentersOld = clustCenters;

%k-means algorithm
iter = 0;
tol = 1e-6; % tolerance for convergence
while(norm(clustCenters-clustCentersOld,'inf')> tol || iter==0)
    iter = iter+1;
    for iPoint=1:nPoints
        [closClusIndex(iPoint,1),~,~] = closestCenter(X(:,iPoint),clustCenters,addElem);
    end
    clustCentersOld = clustCenters;
    clustCenters = zeros(size(clustCenters));
    pointsInCluster = zeros(nClusters,1);
    
    
    for iPoint=1:nPoints
        clustCenters(:,closClusIndex(iPoint,1)) = clustCenters(:,closClusIndex(iPoint,1)) + X(:,iPoint);
        pointsInCluster(closClusIndex(iPoint,1)) = pointsInCluster(closClusIndex(iPoint,1)) + 1;
    end
    
    for iCluster=1:nClusters
        clustCenters(:,iCluster) = clustCenters(:,iCluster) / pointsInCluster(iCluster,1);
    end
end

%Reorder the indices
[~,I] = unique(closClusIndex,'first');
[~,A] = sort(I);
tmp = clustCenters;
tmp2 = closClusIndex;
tmp3 = pointsInCluster;
for i = 1:size(tmp,2)
    clustCenters(:,i) = tmp(:,A(i));
    closClusIndex(tmp2 == A(i)) = i;
    pointsInCluster(i) = tmp3(A(i));
end

% add closest elements to each cluster if required
if addElem
    for iPoint=1:nPoints
        [~,secondClosClusIndex(iPoint,1),secondDist(iPoint,1)] = closestCenter(X(:,iPoint),clustCenters,addElem);
    end
    
    % find how many elements should  be added to each cluster
    additionalPoints(nClusters).Indices = 0;
    for iCluster=1:nClusters
        addPointsInCluster = ceil(addElemThresh*pointsInCluster(iCluster));
        
        % for each cluster find the elements second closest
        IElem = find(secondClosClusIndex==iCluster);
        [~,J] = sort(secondDist(IElem),'ascend');
        % keep the closest
        Iclosest = IElem(J);
        additionalPoints(iCluster).Indices = Iclosest(1:min(addPointsInCluster,length(Iclosest)));
    end
    
else
    additionalPoints = 0;
end
end

function [jmin,jmin2,distmin2] = closestCenter(X,clustCenters,addElem)
%This function determines the closest cluster center to the element X.
%--------------------------------------------------------------------------
%Inputs:
%-------
%X                  - vector containing the state we wish to find the
%                     closest cluster center to.
%clustCenters       - matrix whose columns contain the state vector
%                     corresponding to the cluster centers
%additionalElements - boolean indicating which whether or not to add
%                     additional elements to each cluster 
%
%Outputs:
%--------
%jmin1              - index of the closest cluster center to X
%jmin2              - index of the second closest cluster center to X
%distmin2           - distance between X and the second closest cluster
%                     center. 
%--------------------------------------------------------------------------
%Copyright David Amsallem 2010. Modifications by Matthew Zahr 2011 for
%integration with MORTestbed.

nClusters = size(clustCenters,2);
distCenters = zeros(nClusters,1);
for iCluster=1:nClusters
    distCenters(iCluster,1) = dist(clustCenters(:,iCluster),X);
end
[d,I] = sort(distCenters,'ascend');
jmin = I(1);
if addElem
    jmin2 = I(2);
    distmin2 = d(2);
else
    jmin2 = 0;
    distmin2 = 0;
end
end

function res = dist(X,Y)
%This function computes the 2-norm distance between the vectors or matrices
%X and Y.
%--------------------------------------------------------------------------
%Inputs:
%-------
%X   - vector or matrix
%Y   - vector or matrix
%
%Outputs:
%--------
%res - distance between X and Y
%--------------------------------------------------------------------------
%Copyright David Amsallem 2010. Modifications by Matthew Zahr 2011 for
%integration with MORTestbed.

res = norm(X-Y,2);

% disp('Warning...Change kMeansClusting distance back to the original');
% res = norm((X-mean(X))./norm(X-mean(X),2)-(Y-mean(Y))./norm(Y-mean(Y),2),2);

% normalize = norm(X,2)*norm(Y,2);
% res = acos(abs(X'*Y)/normalize);

end
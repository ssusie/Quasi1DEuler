function  [PHI,S,V,options] = pod(X,k,options)
%This function calculates the Proper Orthogonal Decomposition for a problem
%in the time domain.  If the k is large enough such that nearly zero (<
%1e-15) energy vectors are included, they will be removed and k updated
%accordingly. 
%--------------------------------------------------------------------------
%Inputs:
%-------
%X          - snapshot matrix
%k          - scalar specifying the number of left singular vectors to retain
%approxRank - scalar (or empty); if specified, a rank approxRank SVD of X
%             will be computed instead of the full SVD.  This is useful if
%             X is large and the SVD is too time consuming or memory
%             intensive to use.
%
%Outputs:
%--------
%PHI  - the POD compressed basis
%newk - the true number of left singular vectors retained based on the
%       energy tolerance of the left singular vectors of PHI 
%--------------------------------------------------------------------------

if isempty(X)
    PHI = [];
    S = [];
    V = [];
    return;
end

%Lines 35 - 55 are included purely to handle the case where X may be a cell
%array containing matrices we want to perform the POD on (in this case, k
%may be a vector to indicated the number of left sigular vectors to retain
%for each POD).

if iscell(X)
    %Handle the case where the user passes in a cell array of snapshot
    %matrices, where POD needs to be performed on each.
    PHI = cell(size(X));
    S   = cell(size(X));
    V   = cell(size(X));
    options = cell(size(X));
    if nargin > 1 && length(k) == 1
        khat = k*ones(size(X));
    elseif nargin > 1 && length(k) == length(X)
        khat = k;
    end
    for i = 1:length(X)
        if nargin == 1 || isempty(k)
            %If an order is not indicated, keep all left singular vectors
            khat(i) = min(size(X{i},1),size(X{i},2));
        end
        [PHI{i},S{i},V{i},options{i}] = pod(X{i},khat(i));
    end
    return;
end

%Default values for options
if nargin < 3 || isempty(options)
    options.algorithm='exactSVD';
    options.approxRank=NaN;
    options.powerIters=NaN;
    options.maxvecs=NaN;
    options.truncIntermediate=NaN;
    options.includeSingVals=NaN;
    options.normalize=false;
    options.sqrtMass=[];
%     fprintf('No SVD options specified, using: \n');
%     disp(options);
end
    
%Make sure options has all fields, with appropriate defaults
OptionFields = {'algorithm','approxRank','powerIters','maxvecs','truncIntermediate','includeSingVals','truncFinal','normalize','sqrtMass'};
for i = 1:length(OptionFields)
    if ~isfield(options,OptionFields{i})
        switch i
            case 1
                val = 'exactSVD';
            case 2
                %If algorithm is the low rank SVD approximation, set val to
                %the largest possible rank of the matrix.  Otherwise, set
                %it to NaN
                if strcmpi(options.algorithm,'approxRankSVD')
                    val = min(size(X));
                else
                    val = NaN;
                end
            case 3
                %If algorithm is the low rank SVD approximations, set
                %default number of power iterations to zero. Otherwise, NaN
                if strcmpi(options.algorithm,'approxRankSVD')
                    val = 0;
                else
                    val = NaN;
                end
            case 4
                %If algorithm is limited memory SVD or SVD updates, set
                %maxvecs to the number of columns in X.  Otherwise, NaN.
                if strcmpi(options.algorithm,'incrementalUpdateSVD') || strcmpi(options.algorithm,'limitedMemorySVD')
                    val = size(X,2);
                else
                    val = NaN;
                end
            case 5
                %If algorithm is limited memory SVD or SVD updates, set
                %truncIntermediate to the maximum possible rank of X.
                %Otherwise, NaN. 
                if strcmpi(options.algorithm,'incrementalUpdateSVD') || strcmpi(options.algorithm,'limitedMemorySVD')
                    val = min(size(X));
                else
                    val = NaN;
                end
            case 6
                %If algorithm is limited memory SVD, include singular
                %values.
                if strcmpi(options.algorithm,'limitedMemorySVD')
                    val = true;
                else
                    val = NaN;
                end
            case 7
                val = inf;
            case 8
                val = false;
            case 9
                val = [];
        end
        options=setfield(options,OptionFields{i},val);
    end
end

if options.normalize
    X = bsxfun(@rdivide,X,ColumnwiseNorm(X,2));
end

if ~isempty(options.sqrtMass)
    disp('Here 1');
    X = options.sqrtMass*X;
end

%Compute the SVD of the snapshot matrix
if strcmpi(options.algorithm,'approxRankSVD')
    [U,S,V] = LowRankSVDApprox(X,options);
elseif strcmpi(options.algorithm,'limitedMemorySVD')
    [U,S] = limitedMemorySVD(X,options); V = [];
elseif strcmpi(options.algorithm,'incrementalUpdateSVD')
    [U,S] = incrementalUpdateSVD(X,options); V=[];
elseif strcmpi(options.algorithm,'exactSVD')
    [U,S,V] = svd(X,0);
end

if nargin == 1 || isempty(k)
    %If an order is not indicated, keep all left singular vectors
    k = size(U,2);
end

if k > size(U,2);
    warning('Truncation of Left Singular Vectors is larger than the number of Lefts Singular Values... Returning all Left Singular Vectors');
    k = size(U,2);
end
options.truncFinal=k;

%Store the first k left singular values as our reduced basis (i.e.
%truncated U)
if ~isempty(options.sqrtMass) 
    disp('Here 2');
    PHI = options.sqrtMass\U(:,1:k);
else
    PHI = U(:,1:k);
end

end
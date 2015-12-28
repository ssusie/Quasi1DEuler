function  [] = checkRowColElemData(data,M,N,flag)
%This function issues an error is the size of data is not consistent with
%M, N.
%--------------------------------------------------------------------------
%Input:
%------
%data - matrix, vector, or scalar whose size will be checked for
%       consistency with M and N.
%M    - Number of rows (ignored if empty)
%N    - Number of columns (ignored if empty)
%flag - string used to issue error message if necessary.
%
%Outputs:
%--------
%There are no outputs.  Errors will be generated if necessary.
%--------------------------------------------------------------------------

if isempty(M) && ~isempty(N)
    if ~(size(data,2) == N || size(data,2) == 1) && ~isempty(data)
        error(['Dimensions in ',flag,' field not consistent with numCols']);
    end
elseif ~isempty(M) && isempty(N)
    if ~(size(data,1) == M || size(data,1) == 1) && ~isempty(data)
        error(['Dimensions in ',flag,' field not consistent with numRows']);
    end
else
    if ~(((size(data,1) == M  && size(data,2) == N)) || ((size(data,1) == 1  && size(data,2) == 1))) && ~isempty(data)
        error(['Dimensions in ',flag,' field not consistent with numRows/numCols']);
    end
end

end
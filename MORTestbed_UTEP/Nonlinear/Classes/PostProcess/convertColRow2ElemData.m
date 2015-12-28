function  [ElemData] = convertColRow2ElemData(ColRowData,m,n)
%This function converts the column/row data in ColRowData (cell array) into
%element data (cell array), i.e. if ColRowData is a 1xn cell (i.e. a row),
%it will copy this row m times to make it an mxn cell array.
%--------------------------------------------------------------------------
%Inputs:
%-------
%
%Outputs:
%--------
%ElemData - m x n cell array containing the element data
%--------------------------------------------------------------------------

if length(ColRowData) == 1
    ElemData = repmat(ColRowData,m,n);
    return;
end

if size(ColRowData,1) == m && size(ColRowData,2) == n
    ElemData = ColRowData;
elseif size(ColRowData,1) == m && size(ColRowData,2) == 1
    ElemData = repmat(ColRowData,1,n);
elseif size(ColRowData,1) == 1 && size(ColRowData,2) == n
    ElemData = repmat(ColRowData,m,1);
else
    error('The size of ColRowData is not consistent with m and n');
end

end
function  [out] = appendNum2StrInCellArray(str,arr)
%This function creates a cell array of strings where the (i,j) is str
%appended with the corresponding element from arr (converted to a string)
%--------------------------------------------------------------------------
%Inputs:
%-------
%str - string that will be the root of each entry
%arr - array that defines the suffix of each entry
%
%Outputs:
%--------
%out - cell array containing the information described above
%--------------------------------------------------------------------------

m = size(arr,1);
n = size(arr,2);

out = cell(m,1);
for i = 1:m
    for j = 1:n
        out{i,j} = [str,num2str(arr(i,j))];
    end
end

end
function  [out] = loadObjsIntoCell(fnameRoot,varName,suffArray)
%For each entry in the double array suffArray, this function loads the file
%named [fnameRoot,num2str(suffArray(i)),'.mat'] and saves the variable
%named varName into a cell array entry.
%--------------------------------------------------------------------------
%Inputs:
%-------
%fnameRoot - string containing the root of the list of files to be loaded
%            and stored.  It is assumed that the files have the filename:
%            [fnameRoot,num2str(suffArray(i)),'.mat'].
%varName   - string containing the variable name to be extracted from the
%            file.
%suffArray - double array containing the suffix of the files to load
%
% if you want to load ROM1.mat, ROM2.mat, and ROM5.mat, and save the
% variable rom from each of them in a cell array, then fnameRoot = 'ROM',
% varName = 'rom', and suffArray = [1,2,5]
%
%Outputs:
%--------
%out        - cell array containing the variable varName from each file.
%--------------------------------------------------------------------------

N = length(suffArray);

for i = 1:N
    load([fnameRoot,num2str(suffArray(i)),'.mat']);
    eval(['out{i} = ',varName,';']);
end


end
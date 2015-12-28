function  [cellArray] = convertObjArray2CellArray(objArray)
%This function converts an array of objects into a cell Array of object.
%--------------------------------------------------------------------------
%Inputs:
%-------
%objArray - array of objects to be converted into a cell array of objects
%
%Outputs:
%--------
%cellArray - cell array containing the objects in objArray
%--------------------------------------------------------------------------

cellArray = cell(size(objArray));
for i = 1:length(objArray)
    cellArray{i} = objArray(i);
end

end

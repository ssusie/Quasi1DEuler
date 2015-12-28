function  [] = copyProperties(obj,oldobj)
%This function copies the properties from oldobj to obj, where both oldobj
%and obj are instances of a Problem class.
%--------------------------------------------------------------------------
%Inputs:
%-------
%obj    - Problem class
%oldobj - Problem class whose properties are to be copied into obj
%
%Outputs:
%--------
%There are no outptus.
%--------------------------------------------------------------------------

props = properties(oldobj);
for i = 1:length(props)
    % Use Dynamic Expressions to copy the required property.
    % For more info on usage of Dynamic Expressions, refer to
    % the section "Creating Field Names Dynamically" in:
    % web([docroot '/techdoc/matlab_prog/br04bw6-38.html#br1v5a9-1'])
    obj.(props{i}) = oldobj.(props{i});
end


end
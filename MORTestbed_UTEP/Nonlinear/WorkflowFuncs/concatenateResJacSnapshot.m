function   [res,jac] = concatenateResJacSnapshot(romobjArray)
%This function takes in an array of ROM objects and concatenates the residual
%snapshots into a single matrix (or cell array of matrices if the Local ROB
%method is used).  Same process is also done for jacobian snapshots.
%--------------------------------------------------------------------------
%Inputs:
%-------
%romobjArray - array of ROM objects
%
%Outputs:
%--------
%res - matrix or cell array of matrices containing the concatenated
%      residual snapshots from each rom object
%jac - matrix or cell array of matrices containing the concatenated
%      jacobian snapshots from each rom object
%--------------------------------------------------------------------------

n = 1:length(romobjArray);

if ~iscell(romobjArray(1).res)
    %Handle case where the res/jac snapshots are not in cell arrays (i.e.
    %this is not a local ROB simulation)
    res = [];
    jac = [];
    
    for j = 1:n
        res = [res, romobjArray(j).res];
        jac = [jac, romobjArray(j).jac];
    end
else
    %Handle case where the res/jac snapshots are in cell arrays (i.e.
    %this is a local ROB simulation)
    
    m = romobjArray(1).nBases;
    
    res = cell(m,1);
    jac = cell(m,1);
    
    for i = 1:m
        for j = n
            res{i} = [res{i}, romobjArray(j).res{i}];
            jac{i} = [jac{i}, romobjArray(j).jac{i}];
        end
    end
end
end
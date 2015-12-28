function  [modelobj] = determineModelFromText(ppobj,model,i,j,flag)
%This function determines the model and number of the kth model specified
%in the cell array model.
%--------------------------------------------------------------------------
%Inputs:
%-------
%ppobj    - PostProcess object
%model    - cell array of model type indicators
%i        - row index of model to extract text from
%j        - column index of model to extract text from
%flag     - boolean specifying whether model corresponds to auxilary models
%
%Outputs:
%--------
%modelobj - FOM, ROM, GNAT or locGNAT, or TPWL object
%--------------------------------------------------------------------------

%Extract the entry from the modelobj cell array.
if length(model) == 1
    temp = model{1};
else
    temp = model{i,j};
end

indWord = regexp(temp,'\D');
indDigit = regexp(temp,'\d');
Model = temp(indWord);
num = str2double(temp(indDigit));

if flag
    if strcmpi(Model,'rom')
        modelobj = ppobj.rom{num};
    else
        modelobj = [];
    end
else
    switch lower(Model)
        case 'fom'
            modelobj = ppobj.fom{num};
        case 'rom'
            modelobj = ppobj.rom{num};
        case 'tpwl'
            modelobj = ppobj.tpwl{num};
        case 'gnat'
            modelobj = ppobj.gnat{num};
        case ''
            modelobj = [];
        otherwise
            error('In PP file in AXES field, model must be fom, rom#, gnat#, or tpwl#');
    end
end
end
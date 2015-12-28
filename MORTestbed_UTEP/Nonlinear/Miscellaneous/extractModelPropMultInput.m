function  [varargout] = extractModelPropMultInput(txt,N,ind,flag,def_val)
%This code extracts the model properties when using the multiple input
%notation (i.e. in the ROM file).
%--------------------------------------------------------------------------
%Inputs:
%-------
%txt      - string containing the text from the input file
%N        - number of id's for this model (i.e. number of instances of this
%           model)
%ind      - index of model array to extract 
%flag     - string declaring which property to extract
%def_val  - the default value of the property
%
%Outputs:
%--------
%out - value of the field at the ind_th row
%--------------------------------------------------------------------------

%Create the output
varargout = cell(1,nargout);

%Extract the line of code from txt that is going to be evaluted.
%Note that it should already be in MATLAB syntax.
tmp = determinePropText(txt,flag);
if isempty(tmp)
    %The property was not specified in the input file.  We need to use the
    %default!
    if nargout > 1 && iscell(def_val)
        [varargout{:}] = deal(def_val{:});
    else
        varargout{1} = def_val;
    end
    return;
end

%Evaluate the MATLAB expression in tmp to get a (possibly) vector or cell
%array of entries.  We need to determine which row to use.
% try
vec  = evalin('caller',tmp);
% catch
%     %If evaluation fails, evaluate as text
%     varargout{1} = evalin('caller',['''', tmp,'''']);
%     return;
% end
    
if iscell(vec)
    if size(vec,1) == N && size(vec,1) > 1
        %If the vec has N rows, it is consistent with the size of id, so
        %we just take the appropriate row.
        [varargout{:}] = deal(vec{ind,:});
    elseif size(vec,1) == 1
        %If vec only has one row, use it
        [varargout{:}] = deal(vec{1,:});
    elseif isempty(vec)
        %If vec is empty, set the output to empty
        [varargout{:}] = deal([]);
    else
        %The above cases are the only valid cases, anything else is
        %ambiguous and we must generate and error.
        error('Descriptive error message here!!!!');
    end
else
    if size(vec,1) == N && size(vec,1) > 1
        varargout{1} = vec(ind,:);
    elseif size(vec,1) == 1
        varargout{1} = vec;
    elseif isempty(vec)
        varargout{1} = [];
    else
        error('Descriptive error message here!!!!');
    end
    
    if nargout > 1
        [varargout{:}] = deal(varargout{1});
    end 
end
end
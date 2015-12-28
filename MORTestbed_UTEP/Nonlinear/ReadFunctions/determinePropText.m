function  [propval]= determinePropText(lines,str)
%This function extracts a specific property from the specified object.  If
%str is empty, it returns a cell array where each cell contains the text
%for a certain property (i.e. this gets you all of the properties in case
%the property names are unknown...this was designed to handle the VAR
%section of the input files).
%--------------------------------------------------------------------------
%Inputs:
%-------
%txt      - a string containing the text from the current object (i.e. FOM,
%           GROM, PGROM, GNAT, etc. objects)
%str      - a string indicating the property to be extracted from the
%           current object
%DELIMind - array of indices specifying the locations of the line breaks in
%           txt
%
%Outputs:
%--------
%propval  - char array containing the requested property value of the
%           current object
%--------------------------------------------------------------------------

%Initialize output in case the loop is not entered (i.e. n = 0).
propval = [];

if isempty(str)
    return;
end

n = length(lines);

for i = 1:n
    %Determine if the line contains the appropriate left-hand side.
    txtind = regexp(lines{i},[str,'(\s*)=(\s*)'],'end');
    if length(txtind) > 1
        %Handle the case where there may be more than one equals sign in
        %the line (i.e. title of plot)
        txtind = txtind(1);
    end
    if ~isempty(txtind)
        propval = lines{i}(txtind+1:end);
        return;
    end
end

% %Handle the case where str is empty first.
% if isempty(str)
%     delim  = [0,regexp(txt,'//','end')];
%     
%     propval = cell(length(delim)-1,1);
%     for i = 1:length(propval)
%         propval{i} = txt(delim(i)+1:delim(i+1)-2);
%     end
%     
%     return;
% end

% %Find the location of the property of interest plus its assignment operators
% % txtind = findstr(txt,[str,'=']);
% txtind = regexp(txt,[str,'(\s*)='],'end');
% %If there are no assignment operators, then there are no valid properties
% %and/or no valid property values, so just return an empty matrix
% if isempty(txtind)
%     propval = []; return;
% end
% %Find the line break that corresponds to the line in which the property of
% %interest is contained
% propDELIMind = min(DELIMind(DELIMind>txtind));
% 
% %The property value is all the text in txt after the property's assignment 
% %operator and before its line break 
% % propval = txt(txtind + (length(str)+1): propDELIMind-1);
% propval = txt(txtind+1:propDELIMind-1);
end
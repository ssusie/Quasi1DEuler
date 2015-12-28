function  [] = convertOldInput2New(fname,overwriteFlag)
%This function overwrite the file in fname with the new format.
%Implmentation is very basic and may need to check file to ensure that the
%original right curly brace } is removed.
%--------------------------------------------------------------------------
%Inputs:
%-------
%fname - string containing the file to overwrite with the new input file
%        format.
%overwriteFlag - boolean specifying whether or not to overwrite the current
%                file (if false, a new file named [fname,'.new'] is
%                created.
%
%Outputs:
%--------
%There are no outputs.  The file in fname is overwritten.
%--------------------------------------------------------------------------

keys = {'FOM';'CONFIG';'GROM';'PGROM';'GNAT';...
    'TPWL';'FIGURE';'AXES';'TABLE';'VAR';'GVAR'};

%Store all of file text in MATLAB variable
txt = fileread(fname);

oldEndPP = findstr(txt,'//}}');
ind = union(union(union(oldEndPP,oldEndPP+1),oldEndPP+2),oldEndPP+3);
txt(ind) = [];

oldEnd = findstr(txt,'//}');
ind = union(union(oldEnd,oldEnd+1),oldEnd+2);
txt(ind) = [];

oldDelim = findstr(txt,'//');
ind = union(oldDelim,oldDelim+1);
txt(ind) = [];

for i = 1:length(keys)
    keyindE = regexp(txt,['\<',keys{i}],'end')+1;
    txt(keyindE) = ':';
end

if overwriteFlag
    fid = fopen(fname,'w');
else
    fid = fopen([fname,'.new'],'w');
end
fprintf(fid,txt);
fclose(fid);

% lines = regexp(txt,'\n','split')';
% %Remove blank lines
% for j = length(lines):-1:1
%     if strcmpi(lines{j},'')
%         lines(j) = [];
%     end
% end
% lines = strtrim(lines);
% 
% for i = 1:length(keys)
%     secText = determineSectionText(lines,keys{i},1,keys);
%     
% end


end
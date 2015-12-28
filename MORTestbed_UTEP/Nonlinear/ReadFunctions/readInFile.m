function  [SecText,SecTextChild] = readInFile(sec,fname,ind)
%This function read the text file with filename fname and returns the
%information in this file in MATLAB form
%--------------------------------------------------------------------------
%Inputs:
%-------
%sec     - string specifying the section to extract input text
%fname   - filename of the input file of interest (needs to include .txt
%          extension)
%ind     - scalar indicating the instance of the specific object to extract
%          from the input file and store in obj.
%
%Outputs:
%--------
%SecText - string containing the text in the section specified by the
%          inputs from the specified input file.
%--------------------------------------------------------------------------

%This is a list of the key words in the input files that define blocks
keys = {'FOM:';'CONFIG:';'GROM:';'PGROM:';'GNAT:';...
    'TPWL:';'FIGURE:';'AXES:';'TABLE:';'VAR:';'GVAR:';'OPT:'};

%Store all of file text in MATLAB variable
txt = fileread(fname);

%%%%%%%%%%%%%%%%%Remove Comments%%%%%%%%%%%%%%%%%%%%%
nComm = length(findstr(txt,'/*'));
%Find indices of all begin and end comment specifiers
cStart = findstr(txt,'/*');
cEnd = findstr(txt,'*/');
%If the number of begin comment specifiers (/*) is not the same as the
%number of end comment specifiers (*/), issue an error because it means
%comments were not enclosed in properly
if length(cStart) ~= length(cEnd)
    error('All comments must start with /* and end with */');
end
%Remove comments from readable text
for i = nComm:-1:1
    Begin = cStart(i);
    temp = find(cEnd > Begin);
    End = cEnd(temp(1));
    
    txt(Begin:End+1) = [];
    cStart = findstr(txt,'/*');
    cEnd = findstr(txt,'*/');
end
% %Remove comments from readable text
% for i = nComm:-1:1
%     txt(cStart(i):cEnd(i)+1) = [];
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Combine all multi-line inputs to a single line%%%
mlineS = regexp(txt,'(\.\.\.)\s*\n*\s*','start');
mlineE = regexp(txt,'(\.\.\.)\s*\n*\s*','end');
if length(mlineS) ~= length(mlineE)
    error(['In ',fname,' : Something went wrong with multi-line statements.']);
else
    nMLine = length(mlineS);
end
%Remove ellipses and extra spaces.
for j = nMLine:-1:1
    txt(mlineS(j):mlineE(j)) = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lines = regexp(txt,'\n','split')';
%Remove blank lines
for j = length(lines):-1:1
    if strcmpi(lines{j},'')
        lines(j) = [];
    end
end
lines = strtrim(lines);

switch sec
    case 'FOM'
        %Determine all of the text in the FOM section
        [SecText,SecTextChild] = determineSectionText(lines,'FOM:',1,keys);
    case 'OPT'
        %Determine all of the text in the FOM section
        [SecText,SecTextChild] = determineSectionText(lines,'OPT:',ind,keys);
    case 'CONFIG'
        %Determine all of the text in the CONFIG section
        [SecText,SecTextChild] = determineSectionText(lines,'CONFIG:',ind,keys);
    case 'GROM'
        %Determine all of the text in the GROM section
        [SecText,SecTextChild] = determineSectionText(lines,'GROM:',ind,keys);
    case 'PGROM'
        %Determine all of the text in the PGROM section
        [SecText,SecTextChild] = determineSectionText(lines,'PGROM:',ind,keys);
    case 'GNAT'
        %Determine all of the text in the GNAT section
        [SecText,SecTextChild] = determineSectionText(lines,'GNAT:',ind,keys);
    case 'TPWL'
        %Determine all of the text in the TPWL section
        [SecText,SecTextChild] = determineSectionText(lines,'TPWL:',ind,keys);
    case 'AXES'
        %Determine all of the text in the FIGURE section
        [SecText,SecTextChild] = determineSectionText(lines,'AXES:',ind,keys);
    case 'TABLE'
        %Determine all of the text in the FIGURE section
        [SecText,SecTextChild] = determineSectionText(lines,'TABLE:',ind,keys);
    case 'VAR'
        %Determine all of the text in the VAR section
        [SecText,SecTextChild] = determineSectionText(lines,'VAR:',ind,keys);
    case 'GVAR'
        %Determine all of the text in the VAR section
        [SecText,SecTextChild] = determineSectionText(lines,'GVAR:',ind,keys);
end
end
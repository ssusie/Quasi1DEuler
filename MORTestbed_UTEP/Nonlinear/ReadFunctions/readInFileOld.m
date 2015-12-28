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

%Store all of file text in MATLAB variable
txt = fileread(fname);

if ~(strcmpi(sec,'AXES') ||  strcmpi(sec,'TABLE') || strcmpi(sec,'VAR')) || strcmpi(sec,'GVAR')
    %Remove all white space %
    txt(findstr(txt,' ')) = [];
end


%Remove all comments
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
    txt(cStart(i):cEnd(i)+1) = [];
end

switch sec
    case 'FOM'
        %Determine all of the text in the FOM section
        [SecText,SecTextChild] = determineSectionText(txt,'FOM',1);
    case 'CONFIG'
        %Determine all of the text in the CONFIG section
        [SecText,SecTextChild] = determineSectionText(txt,'CONFIG',ind);
    case 'GROM'
        %Determine all of the text in the GROM section
        [SecText,SecTextChild] = determineSectionText(txt,'GROM',ind);
    case 'PGROM'
        %Determine all of the text in the PGROM section
        [SecText,SecTextChild] = determineSectionText(txt,'PGROM',ind);
    case 'GNAT'
        %Determine all of the text in the GNAT section
        [SecText,SecTextChild] = determineSectionText(txt,'GNAT',ind);
    case 'TPWL'
        %Determine all of the text in the TPWL section
        [SecText,SecTextChild] = determineSectionText(txt,'TPWL',ind);
    case 'AXES'
        %Determine all of the text in the FIGURE section
        [SecText,SecTextChild] = determineSectionText(txt,'AXES',ind);
    case 'TABLE'
        %Determine all of the text in the FIGURE section
        [SecText,SecTextChild] = determineSectionText(txt,'TABLE',ind);
    case 'VAR'
        %Determine all of the text in the VAR section
        [SecText,SecTextChild] = determineSectionText(txt,'VAR',ind);
    case 'GVAR'
        %Determine all of the text in the VAR section
        [SecText,SecTextChild] = determineSectionText(txt,'GVAR',ind);
end
end
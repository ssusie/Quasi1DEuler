function  [SECtext,SECtextChild] = determineSectionText(lines,sec,idnum,keys)
%This function determines all of the text in the section sec.
%--------------------------------------------------------------------------
%Inputs:
%-------
%txt          - string containing the text from the input file
%sec          - string specifying which section to read the objects from
%keys         - cell array of strings defining the keywords that define
%               blocks in the input files
%
%Outputs:
%--------
%SECtext      - string or cell array of strings that contains the text from
%               the objects of the section specified in sec
%SECtextChild - string or cell array of strings that contains the text from
%               the children objects of the section specified in sec
%--------------------------------------------------------------------------

%Loop over all lines and find all instances of the block definition of
%interest (i.e. defined by sec)
% temp = strfind(keys,sec);
% thisBlkInd = [];
% for i = 1:length(keys)
%     if ~isempty(temp{i})
%         thisBlkInd = [thisBlkInd; i];
%     end
% end
%Loop over all lines and find all instances of the block definition of
%interest (i.e. defined by sec)
SECtext = [];
SECtextChild = [];

m = length(keys);
n = length(lines);

blkInd = []; figInd = []; %figInd will only be used if sec is AXES or TABLE
thisBlkInd = [];
temp = cell(n,1);
for i = 1:m
    %Determine the lines with block declarations.
    for k = 1:n
        temp{k} = regexp(lines{k},['\<',keys{i}]);
    end
    for j = 1:n
        if ~isempty(temp{j})
            %If we encouter a block declaration, store its line number
            blkInd = [blkInd; j];
            if strcmpi(keys{i},sec)
                %If we encounter the block declaration sec, store its line
                %number
                thisBlkInd = [thisBlkInd; j];
            elseif strcmpi(keys{i},'FIGURE:')
                %If we encounter the block declaration FIGURE
                figInd = [figInd; j];
            end
        end
    end
end
%Include the last line + 1 in the blkInd (will be helpful when determining
%where to "stop" grabbing lines for the last block defined)
figInd = [figInd; n+1];
blkInd = [sort(blkInd);n+1];
thisBlkInd = sort(thisBlkInd);

if ~strcmpi(sec,'CONFIG:') && ~strcmpi(sec,'AXES:') && ~strcmpi(sec,'TABLE:') && ~strcmpi(sec,'OPT:')
    %CONFIG is the only module where mutliple blocks may be specified.
    %Issue error if the user violates this.
    if length(thisBlkInd) > 1
        error('In an input file, CONFIG & OPT is the only modules where multiple blocks are allowed.');
    end
    %Pull out the lines between this block and the next.  These lines
    %constitute the inside the block.
    if isempty(thisBlkInd)
        return;
    end
    ind = find(blkInd == thisBlkInd);
    SECtext = {lines{thisBlkInd+1:blkInd(ind+1)-1}}';
elseif strcmpi(sec,'CONFIG:') || strcmpi(sec,'OPT:')
    for i = 1:length(thisBlkInd)
        %Find the first section header past the current one (will allow use
        %to define the limits of the block.
        ind = find(blkInd > thisBlkInd(i),1,'first');
        %Store the line indices corresponding to the current block
        lineIndBlock = [thisBlkInd(i)+1, blkInd(ind)-1];
        %Extract the appropriate lines and determine what the id number is
        %for this CFG block
        BlkLines = {lines{lineIndBlock(1):lineIndBlock(2)}}';
        tmp =  determinePropText(BlkLines,'id');
        if isempty(tmp)
            error('In the CFG file, all CONFIG blocks MUST have an id number!');
        end
        ID = eval(tmp);
        %If we have found the configuration with the correct id number, set
        %the output.
        if ID == idnum
            SECtext = BlkLines;
            break;
        end
    end
    %If we did not find a configuration with the right id number, issue a
    %warning
    if isempty(SECtext)
        error(['In the initialize command in your workflow, ',...
              'you did not specify a valid id number (needs to ',...
              'correspond to an id number of one of your configurations ',...
              'in the CFG file)']);
    end
elseif strcmpi(sec,'AXES:') || strcmpi(sec,'TABLE:')
    %In this section, the notation PANEL: will be used to denote either a
    %AXES: or TABLE: block
    
    %Determine the number of FIGURE sections
    p = length(figInd) - 1;
    SECtext = cell(p,1);
    SECtextChild = cell(p,length(thisBlkInd));
    for i = 1:p
        %The following lines will: determine the lines each FIGURE: block
        %and the children of each
        
        %First, get the text from the FIGURE: block
        %Find the first section header past the current FIGURE (will allow use
        %to define the limits of the block.
        ind = find(blkInd > figInd(i),1,'first');
        %Store the line indices corresponding to the current block
        lineIndBlock = [figInd(i)+1, blkInd(ind)-1];
        %Extract the appropriate lines and determine what the id number is
        %for this CFG block
        SECtext{i} = {lines{lineIndBlock(1):lineIndBlock(2)}}';
        
        %Second, determine all PANEL in this FIGURE.  Include the location
        %of the next FIGURE object to determine the bounds of the contents
        %of the last PANEL
        tmp = (thisBlkInd > figInd(i) & thisBlkInd < figInd(i+1));
        lineIndBlock = [thisBlkInd(tmp);figInd(i+1)];
        q = length(lineIndBlock)-1;
        for j = 1:q
            SECtextChild{i,j} = {lines{lineIndBlock(j)+1:lineIndBlock(j+1)-1}}';
        end
    end
end









% %Determine the location of all right curly brackets in the input file
% rCurly = regexpi(txt,'//(\s*)}','end');
% 
% %Determine the location of the beginning of the objects defined in sec (for
% %example, if we are in the GROM section and we have specified a GROM1{} and
% %GROM2{}, then SEClCurly will be a vector with two entries that specify the
% %location in txt of the beginning of these objects.
% switch upper(sec)
%     case 'FOM'
%         SEClCurly = regexpi(txt,'FOM{','end');
%     case 'CONFIG'
%         SEClCurly = regexpi(txt,'CONFIG{','end');
%     case 'GROM'
%         a = regexpi(txt,'GROM{','end');
%         b = regexpi(txt,'PGROM{','end');
%         SEClCurly = setdiff(a,b);
%     case 'PGROM'
%         SEClCurly = regexpi(txt,'PGROM{','end');
%     case 'GNAT'
%         SEClCurly = regexpi(txt,'GNAT{','end');
%     case 'TPWL'
%         SEClCurly = regexpi(txt,'TPWL{','end');
%     case {'AXES','TABLE'}
%         SEClCurly = regexpi(txt,'FIGURE{','end');
%     case 'VAR'
%         a = regexpi(txt,'VAR{','end');
%         b = regexpi(txt,'GVAR{','end');
%         SEClCurly = setdiff(a,b);
%     case 'GVAR'
%         SEClCurly = regexpi(txt,'GVAR{','end');
% end
% 
% if isempty(SEClCurly)
%     SECtext = [];
%     SECtextChild = [];
%     return;
% end
% 
% if strcmpi(sec,'TABLE')
%     sections = {'AXES','TABLE'};
% elseif strcmpi(sec,'AXES')
%     sections = {'TABLE','AXES'};
% end
% 
% if strcmpi(sec,'AXES') || strcmpi(sec,'TABLE')
%     leftEnd = [];
%     rightEnd = [];
%     for z = 1:2
%         SEClCurly_fg = regexpi(txt,'FIGURE{','end');
%         SEClCurly_ax = regexpi(txt,[sections{z},'{'],'end');
%         
%         %For each section beginning, we need to find the appropriate object end
%         %(i.e. right curly bracket)
%         %Allocate vector to contain locations of appropriate object ending
%         %specifiers
%         SECrCurly = zeros(size(SEClCurly_ax));
%         %Allocate cell to contain all text in a specific object
%         SECtextChild = cell(length(SEClCurly_fg),length(SEClCurly_ax));
%         cnt = [];
%         
%         for j = 1:length(SEClCurly_fg)
%             for i = 1:length(SEClCurly_ax) %Loop over all objects of type sec
%                 if ~isempty(intersect(cnt,i))
%                     continue;
%                 end
%                 if j ~= length(SEClCurly_fg) && SEClCurly_ax(i) > SEClCurly_fg(j+1)
%                     break;
%                 end
%                 %The first right curly bracket past the current left curly bracket is
%                 %the appropriate end
%                 SECrCurly(i) = min(rCurly(rCurly > SEClCurly_ax(i)));
%                 %Store the text in the appropriate location
%                 SECtextChild{j,i} = txt(SEClCurly_ax(i)+1:SECrCurly(i)-1);
%                 cnt = [cnt,i];
%                 
%                 leftEnd = [leftEnd; SEClCurly_ax(i)];
%                 rightEnd = [rightEnd; SECrCurly(i)];
%             end
%         end
%     end
%     
%     leftEnd = sort(leftEnd);
%     rightEnd = sort(rightEnd);
%     for i = length(leftEnd):-1:1 %Loop over all objects of type sec
% %     for i = length(SEClCurly_ax):-1:1 %Loop over all objects of type sec
%         %The first right curly bracket past the current left curly bracket is
%         %the appropriate end
%         if txt(leftEnd(i)-1) == 'E'
%             len = 5;
%         elseif txt(leftEnd(i)-1) == 'S'
%             len = 4;
%         end
%         txt(leftEnd(i)-len:rightEnd(i)) = [];
% %         SECrCurly(i) = min(rCurly(rCurly > SEClCurly_ax(i)));
%         %Store the text in the appropriate location
% %         txt(SEClCurly_ax(i)-4:SECrCurly(i)) = [];
%     end
%     
%     txt = removeAllWSpace(txt);
%     rCurly = regexpi(txt,'//}','end');
%     SEClCurly_fg = regexpi(txt,'FIGURE{','end');
%     SECrCurly = zeros(size(SEClCurly_fg));
%     SECtext = cell(length(SEClCurly_fg),1);
%     
%     for j = 1:length(SEClCurly_fg)
%         %The first right curly bracket past the current left curly bracket is
%         %the appropriate end
%         SECrCurly(j) = min(rCurly(rCurly > SEClCurly_fg(j)));
%         SECtext{j,1} = txt(SEClCurly_fg(j)+1:SECrCurly(j)-1);
%     end
%     return;
% end
% 
% SECtextChild = [];
% 
% %For each section beginning, we need to find the appropriate object end
% %(i.e. right curly bracket)
% %Allocate vector to contain locations of appropriate object ending
% %specifiers
% SECrCurly = zeros(size(SEClCurly));
% %Allocate cell to contain all text in a specific object
% SECtextCell = cell(length(SEClCurly),1);
% 
% for i = 1:length(SEClCurly) %Loop over all objects of type sec
%     %The first right curly bracket past the current left curly bracket is
%     %the appropriate end
%     SECrCurly(i) = min(rCurly(rCurly > SEClCurly(i)));
%     %Store the text in the appropriate location
%     SECtextCell{i,1} = txt(SEClCurly(i)+1:SECrCurly(i)-1);
% end
% % SECtext = SECtextCell{ind,1};
% 
% if ~strcmpi(sec,'CONFIG')
%     SECtext = SECtextCell{ind,1};
%     return;
% end
% 
% id = zeros(length(SECtextCell),1);
% for i = 1:length(SECtextCell)
%     delim = findstr(SECtextCell{i},'//');
%     id(i) = eval(determinePropText(SECtextCell{i},'id',delim));
% end
% index = find(id == ind);
% SECtext = SECtextCell{index,1};

end
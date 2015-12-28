function  [SECtext,SECtextChild] = determineSectionText(txt,sec,ind)
%This function determines all of the text in the section sec.
%--------------------------------------------------------------------------
%Inputs:
%-------
%txt          - string containing the text from the input file
%sec          - string specifying which section to read the objects from
%
%Outputs:
%--------
%SECtext      - string or cell array of strings that contains the text from
%               the objects of the section specified in sec
%SECtextChild - string or cell array of strings that contains the text from
%               the children objects of the section specified in sec
%--------------------------------------------------------------------------


%Determine the location of all right curly brackets in the input file
rCurly = regexpi(txt,'//(\s*)}','end');

%Determine the location of the beginning of the objects defined in sec (for
%example, if we are in the GROM section and we have specified a GROM1{} and
%GROM2{}, then SEClCurly will be a vector with two entries that specify the
%location in txt of the beginning of these objects.
switch upper(sec)
    case 'FOM'
        SEClCurly = regexpi(txt,'FOM{','end');
    case 'CONFIG'
        SEClCurly = regexpi(txt,'CONFIG{','end');
    case 'GROM'
        a = regexpi(txt,'GROM{','end');
        b = regexpi(txt,'PGROM{','end');
        SEClCurly = setdiff(a,b);
    case 'PGROM'
        SEClCurly = regexpi(txt,'PGROM{','end');
    case 'GNAT'
        SEClCurly = regexpi(txt,'GNAT{','end');
    case 'TPWL'
        SEClCurly = regexpi(txt,'TPWL{','end');
    case {'AXES','TABLE'}
        SEClCurly = regexpi(txt,'FIGURE{','end');
    case 'VAR'
        a = regexpi(txt,'VAR{','end');
        b = regexpi(txt,'GVAR{','end');
        SEClCurly = setdiff(a,b);
    case 'GVAR'
        SEClCurly = regexpi(txt,'GVAR{','end');
end

if isempty(SEClCurly)
    SECtext = [];
    SECtextChild = [];
    return;
end

if strcmpi(sec,'TABLE')
    sections = {'AXES','TABLE'};
elseif strcmpi(sec,'AXES')
    sections = {'TABLE','AXES'};
end

if strcmpi(sec,'AXES') || strcmpi(sec,'TABLE')
    leftEnd = [];
    rightEnd = [];
    for z = 1:2
        SEClCurly_fg = regexpi(txt,'FIGURE{','end');
        SEClCurly_ax = regexpi(txt,[sections{z},'{'],'end');
        
        %For each section beginning, we need to find the appropriate object end
        %(i.e. right curly bracket)
        %Allocate vector to contain locations of appropriate object ending
        %specifiers
        SECrCurly = zeros(size(SEClCurly_ax));
        %Allocate cell to contain all text in a specific object
        SECtextChild = cell(length(SEClCurly_fg),length(SEClCurly_ax));
        cnt = [];
        
        for j = 1:length(SEClCurly_fg)
            for i = 1:length(SEClCurly_ax) %Loop over all objects of type sec
                if ~isempty(intersect(cnt,i))
                    continue;
                end
                if j ~= length(SEClCurly_fg) && SEClCurly_ax(i) > SEClCurly_fg(j+1)
                    break;
                end
                %The first right curly bracket past the current left curly bracket is
                %the appropriate end
                SECrCurly(i) = min(rCurly(rCurly > SEClCurly_ax(i)));
                %Store the text in the appropriate location
                SECtextChild{j,i} = txt(SEClCurly_ax(i)+1:SECrCurly(i)-1);
                cnt = [cnt,i];
                
                leftEnd = [leftEnd; SEClCurly_ax(i)];
                rightEnd = [rightEnd; SECrCurly(i)];
            end
        end
    end
    
    leftEnd = sort(leftEnd);
    rightEnd = sort(rightEnd);
    for i = length(leftEnd):-1:1 %Loop over all objects of type sec
%     for i = length(SEClCurly_ax):-1:1 %Loop over all objects of type sec
        %The first right curly bracket past the current left curly bracket is
        %the appropriate end
        if txt(leftEnd(i)-1) == 'E'
            len = 5;
        elseif txt(leftEnd(i)-1) == 'S'
            len = 4;
        end
        txt(leftEnd(i)-len:rightEnd(i)) = [];
%         SECrCurly(i) = min(rCurly(rCurly > SEClCurly_ax(i)));
        %Store the text in the appropriate location
%         txt(SEClCurly_ax(i)-4:SECrCurly(i)) = [];
    end
    
    txt = removeAllWSpace(txt);
    rCurly = regexpi(txt,'//}','end');
    SEClCurly_fg = regexpi(txt,'FIGURE{','end');
    SECrCurly = zeros(size(SEClCurly_fg));
    SECtext = cell(length(SEClCurly_fg),1);
    
    for j = 1:length(SEClCurly_fg)
        %The first right curly bracket past the current left curly bracket is
        %the appropriate end
        SECrCurly(j) = min(rCurly(rCurly > SEClCurly_fg(j)));
        SECtext{j,1} = txt(SEClCurly_fg(j)+1:SECrCurly(j)-1);
    end
    return;
end

SECtextChild = [];

%For each section beginning, we need to find the appropriate object end
%(i.e. right curly bracket)
%Allocate vector to contain locations of appropriate object ending
%specifiers
SECrCurly = zeros(size(SEClCurly));
%Allocate cell to contain all text in a specific object
SECtextCell = cell(length(SEClCurly),1);

for i = 1:length(SEClCurly) %Loop over all objects of type sec
    %The first right curly bracket past the current left curly bracket is
    %the appropriate end
    SECrCurly(i) = min(rCurly(rCurly > SEClCurly(i)));
    %Store the text in the appropriate location
    SECtextCell{i,1} = txt(SEClCurly(i)+1:SECrCurly(i)-1);
end
% SECtext = SECtextCell{ind,1};

if ~strcmpi(sec,'CONFIG')
    SECtext = SECtextCell{ind,1};
    return;
end

id = zeros(length(SECtextCell),1);
for i = 1:length(SECtextCell)
    delim = findstr(SECtextCell{i},'//');
    id(i) = eval(determinePropText(SECtextCell{i},'id',delim));
end
index = find(id == ind);
SECtext = SECtextCell{index,1};

end
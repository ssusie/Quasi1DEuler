function  [lineOut] = removeComments(lineIn)

%Initialize arrays that will hold the cell location of a comment begin/end
%and the location in this cell where the comment begins/ends
startInd = [];
endInd = [];

for i = 1:length(lineIn)
    %Loop over all cells and determine if the current lines contains /*
    %and/or*/
    tempS = regexp(lineIn{i},'\/\*');
    tempE = regexp(lineIn{i},'\*\/');
    
    if ~isempty(tempS)
        %Add to the list if we find a comment beginning
        startInd = [startInd; i, tempS];
    end
    if ~isempty(tempE)
        %Add to the list if we find a comment end
        endInd = [endInd; i, tempE];
    end
    
end

%If the startInd and endInd variables do not have the same size, comments
%were not input properly.
if length(startInd) ~= length(endInd)
    error('All comments must start with /* and end with */');
end

for i = size(startInd,1):-1:1
    
    if startInd(i,1) == endInd(i,1)
        %The comment is contained in 1 line
        ind = startInd(i,1);
        cStart = startInd(i,2);
        cEnd = endInd(i,2)+1;
        lineIn{ind}(cStart:cEnd) = [];
        %If line is now empty, remove it.
        if isempty(deblank(lineIn{ind}))
            lineIn(ind) = [];
        end
    else
        %The comment spans multiple lineIn
        indS = startInd(i,1);
        indE = endInd(i,1);
        cStart = startInd(i,2);
        cEnd = endInd(i,2)+1;
        
        %Remove comments from line with /*
        lineIn{indS}(startInd(i,2):end) = [];
        %Empty lines between /* and */
        for j = (indS+1):indE-1
            lineIn(j) = '';
        end
        %Remove commenst from line with */
        lineIn{indE}(1:endInd(i,2)+1) = [];
        
        for k = indE:-1:indS
            if isempty(lineIn{k})
                lineIn(k) = [];
            end
        end
    end
end
lineOut = lineIn;
% %Remove comments
% ind = 1;
% for i = 1:length(lines)
%     %Find the start and end comment flags (if they exist)
%     tmpS = strfind(lines{ind},'/*');
%     tmpE = strfind(lines{ind},'*/');
%     
%     if ~isempty(tmpS) && isempty(tmpE)
%         for j = ind+1:length(lines)
%             tmpE2 = strfind(lines{ind},'*/');
%             if ~isempty(tmpE2)
%                 
%             end
%         end
%     elseif ~isempty(tmpS) && ~isempty(tmpE)
%         %If this line has a single comment
%         lines{ind}(tmpS:tmpE) = [];
%         ind = ind+1;
%     end
% end




end
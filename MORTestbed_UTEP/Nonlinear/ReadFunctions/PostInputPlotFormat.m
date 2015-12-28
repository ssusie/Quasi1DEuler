function  [out] = PostInputPlotFormat(cArray,N,flag)
%This function determines and format the plot object formatting fields from
%the PP file (i.e. plotspec and connCurveSpec)
%--------------------------------------------------------------------------
%Inputs:
%-------
%cArray - a cell array containing the cell array that results from the
%         checkPostInput command
%flag   - string used to determine proper error warning if an error should
%         occur
%
%Outputs:
%--------
%out    - cell array of cells containint the data to be store in either
%         plotspec or connCurveSpec (of the ppAXES class)
%--------------------------------------------------------------------------

if size(cArray,1) == 1
    tmp = cell(N,size(cArray,2));
    for i = 1:N
        [tmp{i,:}] = deal(cArray{1,:});
    end
elseif size(cArray,1) == N
    tmp = cArray;
elseif isempty(cArray)
    tmp = cell(N,1);
    [tmp{:}] = deal('');
else
    if strcmpi(flag,'plotspec')
        error('In the AXES field, plotspec must be either have 1 entry or numPlotObj entries');
    elseif strcmpi(flag,'connCurveSpec')
        error('In the AXES field, connCurveSpec must be either have 1 entry or the same number of entries (cells) as connWithCurve entries');
    end
end

out = cell((size(tmp,1)),1);
for i = 1:size(tmp,1)
    currCell = {};
    for j = 1:size(tmp,2)
        if ~isempty(tmp{i,j})
            currCell = [currCell, tmp{i,j}];
        else
            currCell = [currCell, {}];
        end
    end
    out{i} = currCell;
end
end
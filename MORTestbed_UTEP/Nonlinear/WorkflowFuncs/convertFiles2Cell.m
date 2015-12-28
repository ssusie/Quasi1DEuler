function  [varargout] = convertFiles2Cell(varargin)
%This function loads the indicated files and stores the specified variable
%in a cell array. 
%--------------------------------------------------------------------------
%Inputs:
%-------
%Inputs must come in pairs.  The first entry of a pair in the filename and
%the second is the variable name.  It is assumed that name files exist with
%the same filename core, with a number appended at the end to distinguish
%them.  For example, if there are files
%'ROM1.mat','ROM2.mat',...,'ROM100.mat' that each have a variable named
%rom, then this function will load each file individually and store the rom
%variables in a cell array.  The function call would look as follows:
%cellOfRoms = convertFiles2Cell('ROM','rom');
%This can be done for multiple sequences of files.  To do the sample
%operation for 'ROM1',...,'ROM100.mat' with the variable rom and for
%'GNAT1.mat',...'GNAT500.mat' with the variable gnat, use the command:
%[cellOfRoms,cellOfGnats] = convertFiles2Cell('ROM','rom','GNAT','gnat').
%
%Outputs:
%--------
%Each output is the cell array described above for each input pair.
%--------------------------------------------------------------------------

n = length(varargin)/2;

if round(n) ~= n
    error('There must be an even number of inputs... fName1, varName1,...');
end

for i = 1:n
    cnt = 1;
    varargout{i,1} = [];
    while exist([varargin{2*(i-1) + 1},num2str(cnt),'.mat'],'file')
        load([varargin{2*(i-1) + 1},num2str(cnt),'.mat']);
        varargout{i,1} = [varargout{i,1}, {eval(varargin{2*(i-1) + 2})}];
        cnt = cnt+1;
    end
end
end
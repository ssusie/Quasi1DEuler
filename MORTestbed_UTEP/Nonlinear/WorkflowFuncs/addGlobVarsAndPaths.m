function  [] = addGlobVarsAndPaths(CFGfname)

%Add the main directory to your path
folder = 'Nonlinear';
tmp = pwd;
ind = regexp(tmp,folder);
dir2add = tmp(1:ind+length(folder)-1);
addpath(dir2add);

%Extract parameters from CFG file
GVARtxt = readInFile('GVAR',[CFGfname,'.cfg'],1);

%Extract variables from VARtxt and evaluate them in the current
%workspace.
if ~isempty(GVARtxt)
    for i = 1:length(GVARtxt)
        evalin('caller',[GVARtxt{i},';']);
    end
end
end
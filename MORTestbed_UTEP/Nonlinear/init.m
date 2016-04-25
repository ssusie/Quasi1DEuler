%Script to initialize necessary folders
if ispc
    %If pc, need to use back slash notation
    addpath([pwd,'\ReadFunctions']);
    addpath([pwd,'\Classes']);
    addpath([pwd,'\Classes\Model']);
    addpath([pwd,'\Classes\PostProcess']);
    addpath([pwd,'\Classes\Problem']);
    addpath([pwd,'\Classes\Problem\zFEM']);
    addpath([pwd,'\Classes\Time']);
    addpath([pwd,'\Classes\Optimization']);
    addpath([pwd,'\Classes\OptimizationNew']);
    addpath([pwd,'\InputFileFuncs']);
    addpath([pwd,'\Miscellaneous']);
    addpath([pwd,'\Miscellaneous\distmesh']);
    addpath([pwd,'\Miscellaneous\Radial_Basis_Functions']);
    addpath([pwd,'\TimeIntegrate']);
    addpath([pwd,'\WorkflowFuncs']);
elseif isunix
    %If unix, need to use forward slash notation
    addpath([pwd,'/ReadFunctions']);
    addpath([pwd,'/Classes']);
    addpath([pwd,'/Classes/Model']);
    addpath([pwd,'/Classes/PostProcess']);
    addpath([pwd,'/Classes/Problem']);
    addpath([pwd,'/Classes/Problem/zFEM']);
    addpath([pwd,'/Classes/Time']);
    addpath([pwd,'/Classes/Optimization']);
    addpath([pwd,'/Classes/OptimizationNew']);
    addpath([pwd,'/InputFileFuncs']);
    addpath([pwd,'/Miscellaneous']);
%     addpath([pwd,'/Miscellaneous/distmesh']);
    addpath([pwd,'/Miscellaneous/Radial_Basis_Functions']);
    addpath([pwd,'/TimeIntegrate']);
    addpath([pwd,'/WorkflowFuncs']);
%     addpath('/Users/bokie89/Applications/Ipopt-3.10.3/build64/lib');
%     addpath('/afs/ir/software/matlab-2010b/tomlab');
%     addpath('/afs/ir/software/matlab-2010b/tomlab/mex');
end

clear; close all; clc;
fname = 'eulerExample300dof';
addGlobVarsAndPaths(fname);
%FOM training simulations
nTrain = 4;
targetIndex = 10;
[fom,prob] = initialize(fname,1:nTrain);

sampledParam = zeros(1,nTrain);
sampledInFunc = [];
for i=1:nTrain
    fom(i).executeModel;
    sampledParam(i) = fom(i).prob.pExitIncrease;
end
save([fname,'_FOM.mat'],'fom','prob');

% FOM target simulation
[targetFOM, targetProb] = initialize(fname,targetIndex);
tic;
targetFOM.executeModel;
walltime_fom = toc;
targetParam = targetFOM.prob.pExitIncrease;
targetFOM_sv = targetFOM.sv;

%% 
% Galerkin training simulations
methodROM=1; %solve with constraints  % =1  original,  =2 constrained
numCell=1; %number of domains for constraints
for i=1:nTrain
    galerkin(i) = ROM([fname,'.rom'], 'g', 1, prob(i), methodROM, numCell);
    trainIC = prob(i).ic;
    galerkin(i).clusterSnapsOnly(fom,trainIC);
    galerkin(i).computePOD('fromClustSnaps');
    galerkin(i).executeModel;
end

% Galerkin target simulation
methodROM=1; %solve with constraints  % =1  original,  =2 constrained
numCell=1; %number of domains for constraints
targetGalerkin = ROM([fname,'.rom'], 'g', 1, targetProb, methodROM, numCell);
targetIC = targetProb.ic;
targetGalerkin.clusterSnapsOnly(fom,targetIC);
targetGalerkin.computePOD('fromClustSnaps');
tic;
targetGalerkin.executeModel; 
walltime_galerkin = toc;

%%
% PG training simulations
for i=1:nTrain
    PG(i) = ROM([fname,'.rom'], 'pg', 1, prob(i), methodROM, numCell);
    trainIC = prob(i).ic;
    PG(i).clusterSnapsOnly(fom,trainIC);
    PG(i).computePOD('fromClustSnaps');
    PG(i).executeModel;
end

% PG target simulation
targetPG = ROM([fname,'.rom'], 'pg', 1, targetProb, methodROM, numCell);
targetPG.clusterSnapsOnly(fom,targetIC);
targetPG.computePOD('fromClustSnaps');
tic;
targetPG.executeModel; 
walltime_pg = toc;

%%
% PGcnstd training simulation 
methodROM=2; %solve with constraints  % =1  original,  =2 constrained
numCell=1; %number of domains for constraints
for i=1:nTrain
    PGC(i) = ROM([fname,'.rom'],'pg', 1, prob(i), methodROM, numCell);
    trainIC = prob(i).ic;
    PGC(i).clusterSnapsOnly(fom,trainIC);
    PGC(i).computePOD('fromClustSnaps');
    PGC(i).executeModel;
end

% PGcnstd target simulation
methodROM=2; %solve with constraints  % =1  original,  =2 constrained
numCell=1; %number of domains for constraints
targetPGC = ROM([fname,'.rom'],'pg', 1, targetProb, methodROM, numCell);
targetPGC.clusterSnapsOnly(fom,targetIC);
targetPGC.computePOD('fromClustSnaps');
tic;
targetPGC.executeModel;
walltime_pgc = toc;

%% 
% GNAT training simulations
fnameNLbase='NonlinBase';
NLSnapshot = 0;
methodGNAT = 1; % = 1 original, = 2 Real constraints, =3 Approximated constraints 
augment = false;
readResJacComputePODWriteFcn(fom,fnameNLbase,NLSnapshot,[]);
for i=1:nTrain
    disp(['GNAT training ',num2str(i)]);
    GNAT(i) = genGNAT([fname,'.rom'], PG(i), 1); % rom1 ~ ROM 'pg'
    [phiR,phiJ] = readNonlinBases(fnameNLbase,NLSnapshot,1);
    GNAT(i).createGNAT(prob(i), PG(i), phiR, phiJ, methodGNAT, augment);
    GNAT(i).executeModel;
    GNAT(i).associateFullProblem(prob(i));
    GNAT_sv(i).sv = GNAT(i).reconstructFullState(PG(i));
end

% GNAT target simulation
disp('GNAT target simulation');
targetGNAT = genGNAT([fname,'.rom'], targetPG, 1); % rom1 ~ ROM 'pg'
[phiR,phiJ] = readNonlinBases(fnameNLbase,NLSnapshot,1);
targetGNAT.createGNAT(targetProb(1), targetPG, phiR, phiJ, methodGNAT, augment);
tic;
targetGNAT.executeModel;
walltime_gnat = toc;
targetGNAT.associateFullProblem(targetProb(1));
targetGNAT_sv = targetGNAT.reconstructFullState(targetPG);

%%
% GNATcnstd training simulation
disp('GNAT with real constriants');
fnameNLbase='NonlinBase';
NLSnapshot = 0;
methodGNAT = 2; % = 1 original, = 2 Real constraints, =3 Approximated constraints 
augment = false;
readResJacComputePODWriteFcn(fom,fnameNLbase,NLSnapshot,[]);
for i=1:nTrain
    GC(i) = genGNAT([fname,'.rom'], PGC(i), 1); % ROM with constriants
    [phiR,phiJ] = readNonlinBases(fnameNLbase,NLSnapshot,1);
    GC(i).createGNAT(prob(i), PGC(i), phiR, phiJ, methodGNAT, augment);
    GC(i).executeModel;
    GC(i).associateFullProblem(prob(i));
    GC_sv(i).sv = GC(i).reconstructFullState(PGC(i));
end

% GNATcnstd target simulation
disp('GNAT with real constriants');
fnameNLbase='NonlinBase';
NLSnapshot = 0;
methodGNAT = 2; % = 1 original, = 2 Real constraints, =3 Approximated constraints 
augment = false;
readResJacComputePODWriteFcn(fom,fnameNLbase,NLSnapshot,[]);
targetGC = genGNAT([fname,'.rom'], targetPGC, 1); % ROM with constriants

[phiR,phiJ] = readNonlinBases(fnameNLbase,NLSnapshot,1);
targetGC.createGNAT(targetProb(1), targetPGC, phiR, phiJ, methodGNAT, augment);
tic;
targetGC.executeModel;
walltime_gc = toc;
targetGC.associateFullProblem(targetProb(1));
targetGC_sv = targetGC.reconstructFullState(targetPGC);

%%
% GNATcnstd(PGcnstd) target simulation
disp('GNAT with approxconstriants from PGcnstd snapshots');
targetPG.buildFluxBasisFromTranings(PGC); % build flux basis from PGcnstd snapshots.
fnameNLbase='NonlinBase';
NLSnapshot = 0;
methodGNAT = 3; % = 1 original, = 2 Real constraints, =3 Approximated constraints 
readResJacComputePODWriteFcn(fom,fnameNLbase,NLSnapshot,[]);
targetGCPGC = genGNAT([fname,'.rom'],targetPG,1);
[phiR,phiJ] = readNonlinBases(fnameNLbase,NLSnapshot,1);
targetGCPGC.createGNAT(targetProb(1),targetPG,phiR,phiJ,methodGNAT, augment);
targetGCPGC.preprocessGNATconstrained();
tic;
targetGCPGC.executeModel;
walltime_gcpgc = toc;

targetGCPGC.associateFullProblem(targetProb(1));
targetGCPGC_sv = targetGCPGC.reconstructFullState(targetPG);

%%
% GNATcnstd(GNATcnstd) target simulation
disp('GNAT with approxconstriants from GNATcnstd snapshots');
targetPG.buildFluxBasisFromTranings(GC_sv); % build flux basis from GNATcnstd snapshots.
fnameNLbase='NonlinBase';
NLSnapshot = 0;
methodGNAT = 3; % = 1 original, = 2 Real constraints, =3 Approximated constraints 
readResJacComputePODWriteFcn(fom,fnameNLbase,NLSnapshot,[]);
targetGCGC = genGNAT([fname,'.rom'], targetPG, 1);
[phiR,phiJ] = readNonlinBases(fnameNLbase,NLSnapshot,1);
targetGCGC.createGNAT(targetProb(1), targetPG, phiR, phiJ, methodGNAT, augment);
targetGCGC.preprocessGNATconstrained();
tic;
targetGCGC.executeModel;
walltime_gcgc = toc;

targetGCGC.associateFullProblem(targetProb(1));
targetGCGC_sv = targetGCGC.reconstructFullState(targetPG);

%%
% GNATcnstd(FOM) target simulation
disp('GNAT with approxconstriants from FOM snapshots');
targetPG.buildFluxBasisFromTranings(fom); % build flux basis from FOM snapshots
fnameNLbase='NonlinBase';
NLSnapshot = 0;
methodGNAT = 3; % = 1 original, = 2 Real constraints, =3 Approximated constraints 
readResJacComputePODWriteFcn(fom,fnameNLbase,NLSnapshot,[]);
targetGCF = genGNAT([fname,'.rom'],targetPG,1);
[phiR,phiJ] = readNonlinBases(fnameNLbase,NLSnapshot,1);
targetGCF.createGNAT(targetProb(1),targetPG,phiR,phiJ,methodGNAT,  augment);
targetGCF.preprocessGNATconstrained();
tic;
targetGCF.executeModel;
walltime_gcf = toc;

targetGCF.associateFullProblem(targetProb(1));
targetGCF_sv = targetGCF.reconstructFullState(targetPG);

%%
% GNATcnstd(GNAT) target simulation
disp('GNAT with approxconstriants from GNAT snapshots');
targetPG.buildFluxBasisFromTranings(GNAT_sv); % build flux basis from GNAT snapshots.
fnameNLbase='NonlinBase';
NLSnapshot = 0;
methodGNAT = 3; % = 1 original, = 2 Real constraints, =3 Approximated constraints
readResJacComputePODWriteFcn(fom,fnameNLbase,NLSnapshot,[]);
targetGCG = genGNAT([fname,'.rom'], targetPG, 1);
[phiR,phiJ] = readNonlinBases(fnameNLbase,NLSnapshot,1);
targetGCG.createGNAT(targetProb(1), targetPG, phiR, phiJ, methodGNAT, augment);
targetGCG.preprocessGNATconstrained();
tic;
targetGCG.executeModel;
walltime_gcg = toc;
targetGCG.associateFullProblem(targetProb(1));
targetGCG_sv = targetGCG.reconstructFullState(targetPG);

%%
% GNATcnstd(PG) target simulation
disp('GNAT with approxconstriants from PG snapshots');
targetPG.buildFluxBasisFromTranings(PG); % build flux basis from PG
fnameNLbase='NonlinBase';
NLSnapshot = 0;
methodGNAT = 3; % = 1 original, = 2 Real constraints, =3 Approximated constraints 
readResJacComputePODWriteFcn(fom,fnameNLbase,NLSnapshot,[]);
targetGCPG = genGNAT([fname,'.rom'], targetPG, 1);
[phiR,phiJ] = readNonlinBases(fnameNLbase,NLSnapshot,1);
targetGCPG.createGNAT(targetProb(1), targetPG, phiR, phiJ, methodGNAT, augment);
targetGCPG.preprocessGNATconstrained();
tic;
targetGCPG.executeModel;
walltime_gcpg = toc;
targetGCPG.associateFullProblem(targetProb(1));
targetGCPG_sv = targetGCPG.reconstructFullState(targetPG);

fprintf('Wall time for FOM is %e\n',walltime_fom);
fprintf('Wall time for G is %e\n',walltime_galerkin);
fprintf('Wall time for PG is %e\n',walltime_pg);
fprintf('Wall time for PGcnstd is %e\n',walltime_pgc);
fprintf('Wall time for GNAT is %e\n',walltime_gnat);
fprintf('Wall time for GNATcnstd is %e\n',walltime_gc);
fprintf('Wall time for GNATcnstd(FOM) is %e\n',walltime_gcf);
fprintf('Wall time for GNATcnstd(PG) is %e\n',walltime_gcpg);
fprintf('Wall time for GNATcnstd(GNAT) is %e\n',walltime_gcg);
fprintf('Wall time for GNATcnstd(PGcnstd) is %e\n',walltime_gcpgc);
fprintf('Wall time for GNATcnstd(GNATcnstd) is %e\n\n',walltime_gcgc);


fprintf('Speedup for G is %e\n',walltime_fom/walltime_galerkin);
fprintf('Speedup for PG is %e\n',walltime_fom/walltime_pg);
fprintf('Speedup for PGcnstd is %e\n',walltime_fom/walltime_pgc);
fprintf('Speedup for GNAT is %e\n',walltime_fom/walltime_gnat);
fprintf('Speedup for GNATcnstd is %e\n',walltime_fom/walltime_gc);
fprintf('Speedup for GNATcnstd(FOM) is %e\n',walltime_fom/walltime_gcf);
fprintf('Speedup for GNATcnstd(PG) is %e\n',walltime_fom/walltime_gcpg);
fprintf('Speedup for GNATcnstd(GNAT) is %e\n',walltime_fom/walltime_gcg);
fprintf('Speedup for GNATcnstd(PGcnstd) is %e\n',walltime_fom/walltime_gcpgc);
fprintf('Speedup for GNATcnstd(GNATcnstd) is %e\n',walltime_fom/walltime_gcgc);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('OneSubdomain/solutions','targetFOM', 'targetGalerkin', 'targetPG', 'targetPGC',...
        'targetGNAT', 'targetGC', 'targetGCPGC', 'targetGCGC',...
        'targetGCF', 'targetGCG', 'targetGCPG',...
        'targetGNAT_sv', 'targetGC_sv', 'targetGCPGC_sv', 'targetGCGC_sv',...
        'targetGCF_sv', 'targetGCG_sv', 'targetGCPG_sv', 'targetProb');
plotAll(targetFOM, targetGalerkin, targetPG, targetPGC,...
        targetGNAT, targetGC, targetGCPGC, targetGCGC,...
        targetGCF, targetGCG, targetGCPG,...
        targetGNAT_sv, targetGC_sv, targetGCPGC_sv, targetGCGC_sv,...
        targetGCF_sv, targetGCG_sv, targetGCPG_sv, targetProb);
% plotAll(targetFOM, targetProb);

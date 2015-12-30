clear; close all; clc;
fname = 'eulerExample300dof';
addGlobVarsAndPaths(fname);

%FOM
[fom,prob] = initialize(fname,1);
fom.executeModel;

%%
methodROM=2; % 1 = original,  2 = constrained
numCell=5; %number of domains for constraints
basisNumber=16; %truncate basis Phi; if don't want to truncate choose basisNumber=fom.cTimeIter
rom = ROM([fname,'.rom'],'pg',1,prob,methodROM,numCell,basisNumber);
trainIC = prob.ic;
rom.clusterSnapsOnly(fom,trainIC);
rom.computePOD('fromClustSnaps');
rom.executeModel;

%%
%GNAT
fnameNLbase='NonlinBase';
NLSnapshot=0;
methodGNAT=3; % =1 original, =2 Rom constrains; =3 Gnat constraints

readResJacComputePODWriteFcn(fom,fnameNLbase,NLSnapshot,[]);
gnat = genGNAT([fname,'.rom'],rom,1);
[phiR,phiJ] = readNonlinBases(fnameNLbase,NLSnapshot,1);
gnat.createGNAT(prob(1),rom,phiR,phiJ,methodGNAT,1);
gnat.executeModel;

gnat.associateFullProblem(prob(1));
svG = gnat.reconstructFullState(rom);
%%
%%%Postprocess%%%
RomErr = mean(ColumnwiseNorm(rom.sv(:,end)-fom.sv(:,end),2)./ColumnwiseNorm(fom.sv(:,end),2));
fprintf('Average Relative L2 ROM Error = %f %%\n',RomErr);
GnatErr = mean(ColumnwiseNorm(svG(:,end)-fom.sv(:,end),2)./ColumnwiseNorm(fom.sv(:,end),2));
fprintf('Average Relative L2 GNAT Error = %f %%\n',GnatErr);

[rhoF,uF,PF,cF,eF] = prob.getVariables(fom.sv(:,end));
[rhoR,uR,PR,cR,eR] = prob.getVariables(rom.sv(:,end));
[rhoG,uG,PG,cG,eG] = prob.getVariables(svG(:,end));

figure;
hfom  = plot(uF./cF,'k','linewidth',2); hold on;
hrom  = plot(uR./cR,'b--','linewidth',2);
hgnat = plot(uG./cG,'g');
xlabel('x')
ylabel('Mach')
legend([hfom,hrom,hgnat],'FOM','ROM','GNAT')

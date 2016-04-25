clear; close all; clc;
fname = 'eulerExample300dof';
addGlobVarsAndPaths(fname);

%FOM
[fom,prob] = initialize(fname,1);
fom.executeModel;
soln=fom.sv;
save fom_soln soln

%%
methodROM=2; %solve with constraints  % =1  original,  =2 constrained
numCell=1; %number of domains for constraints
basisNumber=10; %truncate basis Phi; if don't want to truncate choose basisNumber=fom.cTimeIter
rom = ROM([fname,'.rom'],'pg',1,prob,methodROM,numCell,basisNumber);
trainIC = prob.ic;
rom.clusterSnapsOnly(fom,trainIC);
rom.computePOD('fromClustSnaps');
rom.executeModel;


%%
%GNAT
fnameNLbase='NonlinBase';
NLSnapshot = 0;
methodGNAT = 4; % = 1 original, = 2 Rom constrains; = 3 Gnat constraints

readResJacComputePODWriteFcn(fom,fnameNLbase,NLSnapshot,[]);
gnat1 = genGNAT([fname,'.rom'],rom,1);
[phiR,phiJ] = readNonlinBases(fnameNLbase,NLSnapshot,1);
gnat1.createGNAT(prob(1),rom,phiR,phiJ,methodGNAT,1);
gnat1.executeModel;

gnat1.associateFullProblem(prob(1));
svG1 = gnat.reconstructFullState(rom);


%%
%%%Postprocess%%%
RomErr = mean(ColumnwiseNorm(rom.sv(:,end)-fom.sv(:,end),2)./ColumnwiseNorm(fom.sv(:,end),2));
fprintf('Average Relative L2 ROM Error = %f %%\n',RomErr);
GnatErr = mean(ColumnwiseNorm(svG(:,end)-fom.sv(:,end),2)./ColumnwiseNorm(fom.sv(:,end),2));
fprintf('Average Relative L2 GNAT Error = %f %%\n',GnatErr);
GnatErr1 = mean(ColumnwiseNorm(svG1(:,end)-fom.sv(:,end),2)./ColumnwiseNorm(fom.sv(:,end),2));
fprintf('Average Relative L2 GNAT Error = %f %%\n',GnatErr1);


[rhoF,uF,PF,cF,eF] = prob.getVariables(fom.sv(:,end));
[rhoR,uR,PR,cR,eR] = prob.getVariables(rom.sv(:,end));
[rhoG,uG,PG,cG,eG] = prob.getVariables(svG(:,end));
[rhoG1,uG1,PG1,cG1,eG1] = prob.getVariables(svG1(:,end));

figure;
hfom  = plot(uF./cF,'k','linewidth',2); hold on;
hrom  = plot(uR./cR,'b--','linewidth',2);
hgnat = plot(uG./cG,'c');
hgnat1 = plot(uG1./cG1,'g');

xlabel('x')
ylabel('Mach')
legend([hfom,hrom,hgnat, hgnat1],'FOM','ROM','GNATRom', 'GNATapprox')
%%
% norm of the real and approx constraints
figure(2)
% subplot(2,1,1)
plot(gnat.Anorm,'ko-', 'linewidth', 2), hold on
plot(gnat.Rnorm,'b*-', 'linewidth', 2)
legend('norm approx constr','norm real constr')
title('using Method 4')
% subplot(2,1,2)
% plot(gnat.Anorm(1:end-2),'ko-', 'linewidth', 2), hold on
% plot(gnat.Rnorm(1:end-2),'b*-', 'linewidth', 2)
% legend('norm approx constr','norm real constr')

%%
figure(3)
% subplot(2,1,1)
plot(gnat1.Anorm1,'ko-', 'linewidth', 2), hold on
plot(gnat1.Rnorm1,'b*-', 'linewidth', 2)
legend('norm approx constr','norm real constr')
title('using Method 5 with new snapshots')





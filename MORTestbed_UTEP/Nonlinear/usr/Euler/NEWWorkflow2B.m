clear; close all; clc;
fname = 'eulerExample300dof';
addGlobVarsAndPaths(fname);

%FOM
[fom,prob] = initialize(fname,1);
fom.executeModel;
fom_soln=fom.sv;
save fom_soln fom_soln

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
methodGNAT = 3; % = 1 original, = 2 Rom constrains; = 3 Gnat constraints
augment = false;

readResJacComputePODWriteFcn(fom,fnameNLbase,NLSnapshot,[]);
gnat_approx = genGNAT([fname,'.rom'],rom,1);
[phiR,phiJ] = readNonlinBases(fnameNLbase,NLSnapshot,1);
gnat_approx.createGNAT(prob(1),rom,phiR,phiJ,methodGNAT, augment, 1);
gnat_approx.executeModel;

gnat_approx.associateFullProblem(prob(1));
svG_approx = gnat_approx.reconstructFullState(rom);

%%
%GNAT
fnameNLbase='NonlinBase';
NLSnapshot = 0;
methodGNAT = 2; % = 1 original, = 2 Rom constrains; = 3 Gnat constraints
augment = false;

readResJacComputePODWriteFcn(fom,fnameNLbase,NLSnapshot,[]);
gnat_rom = genGNAT([fname,'.rom'],rom,1);
[phiR,phiJ] = readNonlinBases(fnameNLbase,NLSnapshot,1);
gnat_rom.createGNAT(prob(1),rom,phiR,phiJ,methodGNAT, augment, 1);
gnat_rom.executeModel;

gnat_rom.associateFullProblem(prob(1));
svG_rom = gnat_rom.reconstructFullState(rom);



%%
%GNAT
fnameNLbase='NonlinBase';
NLSnapshot = 0;
methodGNAT = 4; % = 1 original, = 2 Rom constrains; = 3 Gnat constraints

readResJacComputePODWriteFcn(fom,fnameNLbase,NLSnapshot,[]);
gnat_approx_rom = genGNAT([fname,'.rom'],rom,1);
[phiR,phiJ] = readNonlinBases(fnameNLbase,NLSnapshot,1);
gnat_approx_rom.createGNAT(prob(1),rom,phiR,phiJ,methodGNAT,augment, 1);
gnat_approx_rom.executeModel;

gnat_approx_rom.associateFullProblem(prob(1));
svGapprox_approx_rom = gnat_approx_rom.reconstructFullState(rom);


%%
%%%Postprocess%%%
RomErr = mean(ColumnwiseNorm(rom.sv(:,end)-fom.sv(:,end),2)./ColumnwiseNorm(fom.sv(:,end),2));
fprintf('Average Relative L2 ROM Error = %f %%\n',RomErr);
GnatErr_approx = mean(ColumnwiseNorm(svG_approx(:,end)-fom.sv(:,end),2)./ColumnwiseNorm(fom.sv(:,end),2));
fprintf('Average Relative L2 GNAT Error = %f %%\n',GnatErr_approx);
GnatErr_rom = mean(ColumnwiseNorm(svG_rom(:,end)-fom.sv(:,end),2)./ColumnwiseNorm(fom.sv(:,end),2));
fprintf('Average Relative L2 GNAT Error = %f %%\n',GnatErr_rom);
GnatErr_rom_approx = mean(ColumnwiseNorm(svGapprox_approx_rom(:,end)-fom.sv(:,end),2)./ColumnwiseNorm(fom.sv(:,end),2));
fprintf('Average Relative L2 GNAT Error = %f %%\n',GnatErr_rom_approx);


[rhoF,uF,PF,cF,eF] = prob.getVariables(fom.sv(:,end));
[rhoR,uR,PR,cR,eR] = prob.getVariables(rom.sv(:,end));
[rhoG_approx,uG_approx,PG_approx,cG_approx,eG_approx] = prob.getVariables(svG_approx(:,end));
[rhoG_rom,uG_rom,PG_rom,cG_rom,eG_rom] = prob.getVariables(svG_rom(:,end));
[rhoG_rom_approx,uG_rom_approx,PG_rom_approx,cG_rom_approx,eG_rom_approx] = ...
    prob.getVariables(svGapprox_approx_rom(:,end));
% 
figure(10);
hfom  = plot(uF./cF,'k','linewidth',2); hold on;
hrom  = plot(uR./cR,'b','linewidth',2);
hgnat = plot(uG_approx./cG_approx,'c--');
hgnat1 = plot(uG_rom./cG_rom,'g--');
hgnat2 = plot(uG_rom_approx./cG_rom_approx,'m--');

xlabel('x')
ylabel('Mach')
legend([hfom,hrom,hgnat,hgnat1, hgnat2],'FOM', 'ROM', 'GNAT ROM snapsh',...
    'GNAT_real', 'GNAT GNAT snapshots')

%%
% norm of the real and approx constraints
figure(2)
plot(gnat_approx.Anorm,'co-', 'linewidth', 2), hold on
plot(gnat_rom.Rnorm,'g*-', 'linewidth', 2), hold on
plot(gnat_approx_rom.Anorm1,'mo-', 'linewidth', 2)

legend('ROM snapsh','real constr','GNAT snapsh')
% legend('norm real constr','norm of approx with real snapsh')
% title('using Method 5')
% subplot(2,1,2)
% plot(gnat.Anorm(1:end-2),'ko-', 'linewidth', 2), hold on
% plot(gnat.Rnorm(1:end-2),'b*-', 'linewidth', 2)
% legend('norm approx constr','norm real constr')

%%
% figure(3)
% % subplot(2,1,1)
% plot(gnat1.Anorm1,'ko-', 'linewidth', 2), hold on
% plot(gnat1.Rnorm1,'b*-', 'linewidth', 2)
% legend('norm approx constr','norm real constr')
% title('using Method 5 with new snapshots')
% 
% 



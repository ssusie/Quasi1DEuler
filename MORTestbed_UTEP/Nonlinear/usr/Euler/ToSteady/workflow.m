clear; close all; clc;
fname = 'eulerExample300dof';
addGlobVarsAndPaths(fname);

%FOM
[fom,prob] = initialize(fname,1);
S = fom.prob.S;
% save('mySfom.mat', 'S');
% return;
fom.executeModel;
soln=fom.sv;
save soln soln
 
%   keyboard
%%
rom = ROM([fname,'.rom'],'pg',1,prob,1,10);
trainIC = prob.ic;
rom.clusterSnapsOnly(fom,trainIC);
rom.computePOD('fromClustSnaps');
rom.executeModel;

% save CZim5
% norm(fom.sv-rom.sv,'fro')/norm(fom.sv,'fro')
% keyboard
% break
%%
%ROM
% clear FnError
% for basNum=1:25
%     for cellNum=1:10
%         rom = ROM([fname,'.rom'],'pg',1,prob,cellNum, basNum);
%         trainIC = prob.ic;
%         rom.clusterSnapsOnly(fom,trainIC);
%         rom.computePOD('fromClustSnaps');
%         rom.executeModel;
%         Rom_sv(cellNum,basNum,:,:)=rom.sv;
%         save ROMsv Rom_sv
%         Z=fom.sv-rom.sv;
% %         nError(cellNum, :)=diag(Z'*Z)';
%         FnError(cellNum,basNum)=norm(Z,'fro')/norm(fom.sv,'fro')
%         save FnError FnError
%      end
% end
% % save nError nError
% 
%    break
% keyboard
%%
%GNAT
fnameNLbase='NonlinBase';
NLSnapshot=0;

readResJacComputePODWriteFcn(fom,fnameNLbase,NLSnapshot,[]);

gnat = genGNAT([fname,'.rom'],rom,1);
[phiR,phiJ] = readNonlinBases(fnameNLbase,NLSnapshot,1);
gnat.createGNAT(prob(1),rom,phiR,phiJ,1);
% keyboard
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
% legend([hfom,hrom,hgnat],'FOM','ROM','GNAT')
break
%%
% close all

original3=load('original50_3');
constrained3=load('constrained50_3');

original4=load('original50_4');
constrained4=load('constrained50_4');

original5=load('original50_5');
constrained5=load('constrained50_5');

[rhoO3,uO3,PO3,cO3,eO3] = prob.getVariables(original3.rom.sv(:,end));
[rhoC3,uC3,PC3,cC3,eC3] = prob.getVariables(constrained3.rom.sv(:,end));

[rhoO4,uO4,PO4,cO4,eO4] = prob.getVariables(original4.rom.sv(:,end));
[rhoC4,uC4,PC4,cC4,eC4] = prob.getVariables(constrained4.rom.sv(:,end));

[rhoO5,uO5,PO5,cO5,eO5] = prob.getVariables(original5.rom.sv(:,end));
[rhoC5,uC5,PC5,cC5,eC5] = prob.getVariables(constrained5.rom.sv(:,end));

h3=figure(3);
horig  = plot(uO3./cO3,'k--','linewidth',2); hold on;
hconstr  = plot(uC3./cC3,'b--','linewidth',2);
hfom  = plot(uF./cF,'r','linewidth',2); 
xlabel('x')
ylabel('Mach')
legend([horig,hconstr, hfom],'Original','Constrained', 'FOM')
title('3 basis vactors')

h4=figure(4);
horig  = plot(uO4./cO4,'k--','linewidth',2); hold on;
hconstr  = plot(uC4./cC4,'b--','linewidth',2);
hfom  = plot(uF./cF,'r','linewidth',2); hold on;
xlabel('x')
ylabel('Mach')
legend([horig,hconstr,hfom],'Original','Constrained', 'FOM')
title('4 basis vactors')

h5=figure(5);
horig  = plot(uO5./cO5,'k--','linewidth',2); hold on;
hconstr  = plot(uC5./cC5,'b--','linewidth',2);
hfom  = plot(uF./cF,'r','linewidth',2); 
xlabel('x')
ylabel('Mach')
legend([horig,hconstr, hfom],'Original','Constrained', 'FOM')
title('5 basis vactors')

% 
% saveas(h3, 'basis3', 'eps')
% saveas(h4, 'basis4', 'eps')
% saveas(h5, 'basis5', 'eps')
% 








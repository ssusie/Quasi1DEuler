%Workflow 1
% addpath('/Users/bokie89/Desktop/MORTestbed_main/ModularizedFolderNonlin');
fname = 'fhnexample1024dof';
addGlobVarsAndPaths(fname);

%%%FOM%%%
[fom,prob] = initialize(fname,1);
fom.executeModel;

%%%ROM%%%
rom = ROM([fname,'.rom'],'pg',1,prob);

trainIC = prob.ic;
rom.clusterSnapsOnly(fom,trainIC);
rom.computePOD('fromClustSnaps');
rom.executeModel();

%%%GNAT%%%
fnameNLbase='romNonlinBase';
NLSnapshot=1;
readResJacComputePODWriteFcn(rom,fnameNLbase,NLSnapshot,[]);

gnat = genGNAT([fname,'.rom'],rom,1);
[phiR,phiJ] = readNonlinBases(fnameNLbase,1,1);
gnat.createGNAT(prob,rom,phiR,phiJ,1);
gnat.executeModel;

gnat.associateFullProblem(prob);
svG = gnat.reconstructFullState(rom);

RomErr = mean(ColumnwiseNorm(rom.sv(:,2:end)-fom.sv(:,2:end),2)./ColumnwiseNorm(fom.sv(:,2:end),2));
fprintf('Average Relative L2 ROM Error = %f %%\n',RomErr);
GnatErr = mean(ColumnwiseNorm(svG(:,2:end)-fom.sv(:,2:end),2)./ColumnwiseNorm(fom.sv(:,2:end),2));
fprintf('Average Relative L2 GNAT Error = %f %%\n',GnatErr);


figure;
midv = fom.ndof/4;
midw = 3*midv;

subplot(2,1,1);
hfom  = plot(fom.sv(midv,:),'k','linewidth',2); hold on;
hrom  = plot(rom.sv(midv,:),'b--','linewidth',2);
hgnat = plot(svG(midv,:),'g');
xlabel('t')
ylabel('v(mid)')
legend([hfom,hrom,hgnat],'FOM','ROM','GNAT');

subplot(2,1,2);
hfom  = plot(fom.sv(midw,:),'k','linewidth',2); hold on;
hrom  = plot(rom.sv(midw,:),'b--','linewidth',2);
hgnat = plot(svG(midw,:),'g');
xlabel('t')
ylabel('w(mid)')
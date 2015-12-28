%Workflow 1
fname = 'nltEx100dof';
addGlobVarsAndPaths(fname);

%%%FOM%%%
% [fom,prob] = initialize(fname,1:4);
[fom,prob] = initialize(fname,1:4);
for i = 1:length(fom)
    fom(i).executeModel;
end

%%%ROM%%%
rom = ROM([fname,'.rom'],'pg',1,prob(4));

trainIC = prob(1).ic;
rom.clusterSnapsOnly(fom(1:3),trainIC);
rom.computePOD('fromClustSnaps');
rom.executeModel;

%%%GNAT%%%
fnameNLbase='romNonlinBase';
NLSnapshot=1;
readResJacComputePODWriteFcn(rom,fnameNLbase,NLSnapshot,[]);
gnat = genGNAT([fname,'.rom'],rom,1);
[phiR,phiJ] = readNonlinBases(fnameNLbase,1,1);
gnat.createGNAT(prob(4),rom,phiR,phiJ,1);
gnat.executeModel;
gnat.associateFullProblem(prob(4));
svG = gnat.reconstructFullState(rom);

%%%Postprocess%%%
RomErr = mean(ColumnwiseNorm(rom.sv(:,2:end)-fom(4).sv(:,2:end),2)./ColumnwiseNorm(fom(4).sv(:,2:end),2));
fprintf('Average Relative L2 ROM Error = %f %%\n',RomErr);
GnatErr = mean(ColumnwiseNorm(svG(:,2:end)-fom(4).sv(:,2:end),2)./ColumnwiseNorm(fom(4).sv(:,2:end),2));
fprintf('Average Relative L2 GNAT Error = %f %%\n',GnatErr);

figure;
hfom  = plot(fom(4).sv(1,:),'k','linewidth',2); hold on;
hrom  = plot(rom.sv(1,:),'b--','linewidth',2);
hgnat = plot(svG(1,:),'g');
xlabel('t')
ylabel('v(1)')
legend([hfom,hrom,hgnat],'FOM','ROM','GNAT')
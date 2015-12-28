% addpath('/Users/bokie89/Desktop/MORTestbed_main/ModularizedFolderNonlin');
fname = 'memsEx4050dof';
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


N = fom.prob.N;
figure;
hfom  = plot(fom.sv(1:N,end),'k','linewidth',2); hold on;
hrom  = plot(rom.sv(1:N,end),'b--','linewidth',2);
hgnat = plot(svG(1:N,end),'g');
xlabel('x')
ylabel('Beam position at final time')
legend([hfom,hrom,hgnat],'FOM','ROM','GNAT')


return;

save('FOM.mat','fom','prob');

for k = 1:3
    load('FOM.mat');
    
    %%%%%%%%%%ROM%%%%%%%%%%
    rom = ROM([fname,'.rom'],'pg',k,prob);
    
    trainIC = prob.ic;
    X = rom.determineSnapshots(fom,trainIC);
    rom.computePOD(X,trainIC);
    rom.executeModel;
    
    save(['ROM',num2str(k),'.mat'],'rom');
    
    %%%%%%%%%%Local GNAT%%%%%%%%%%
    [phiR,~] = pod(rom.res);
    [phiJ,~] = pod(rom.jac);
    
    gnat = genGNAT([fname,'.rom'],rom,k);
    gnat.createGNAT(prob,phiR,phiJ,rom.phi);
    clearAllExcept('fname','gnat','k');
    gnat.executeModel;
    load('FOM.mat');
    gnat.associateFullProblem(prob);
    save(['GNAT',num2str(k),'.mat'],'gnat');
    clearAllExcept('fname','k');
end

load FOM;
fomArray = {fom};
romArray = loadObjsIntoCell('ROM','rom',1:3);
gnatArray = loadObjsIntoCell('GNAT','gnat',1:3);

pp1 = PostProcess('MEMS.pp',fomArray,romArray,gnatArray,[]);
pp2 = PostProcess('MEMS2.pp',fomArray,romArray,gnatArray,[]);
pp3 = PostProcess('test.pp',fomArray,romArray,gnatArray,[]);

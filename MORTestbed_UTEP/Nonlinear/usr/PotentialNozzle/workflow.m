%Workflow 1
fname = 'SteadyNozzle';
addGlobVarsAndPaths(fname);

nnSp = 10;
xx = linspace(0,2,nnSp)';f0 = 1; f1 = 0.4; f2 = 1;
a = [1, 1, 1; 3, 2, 1; 8, 4, 2]\[f1-f0;0;f2-f0];
param = a(1)*xx.^3 + a(2)*xx.^2 + a(3)*xx + f0;

%Full Order Model
[fom,prob] = initialize(fname,1); 
prob.updateParameters(param);
fom.executeModel;

save('SteadyFOM.mat','fom','prob');
% load('SteadyFOM.mat');
xi = 0.5*(fom(1).prob.mesh.node(2:end,1) + fom(1).prob.mesh.node(1:end-1,1));

uF = fom.prob.velocity([prob.U0;fom.sv(:,end);prob.UL]);
rhoF = fom.prob.density(uF);
MF = fom.prob.mach(rhoF,uF);

figure(); ax=axes();
uF = fom.prob.velocity([prob.U0;fom.sv(:,end);prob.UL]);
rhoF = fom.prob.density(uF);
MF = fom.prob.mach(rhoF,uF);
set(gca,'nextplot','replacechildren','ylim',[0,max(MF)]);
for i = 1:fom.time.nstep
    uF = fom.prob.velocity([prob.U0;fom.sv(:,i);prob.UL]);
    rhoF = fom.prob.density(uF);
    MF = fom.prob.mach(rhoF,uF);
   
    plot(xi,MF,'linewidth',2); drawnow;
end


%Reduced Order Model - Galerkin
romG = ROM([fname,'.rom'],'g',1,prob(1));
trainIC = prob(1).ic;
romG.clusterSnapsOnly(fom,trainIC);
romG.computePOD('fromClustSnaps');
romG.executeModel;

figure(); ax=axes();
uG = fom.prob.velocity([prob.U0;romG.sv(:,end);prob.UL]);
rhoG = fom.prob.density(uG);
MG = fom.prob.mach(rhoG,uG);
set(gca,'nextplot','replacechildren','ylim',[0,max(MG)]);
for i = 1:romG.time.nstep
    uG = romG.prob.velocity([prob.U0;romG.sv(:,i);prob.UL]);
    rhoG = romG.prob.density(uG);
    MG = romG.prob.mach(rhoG,uG);
    plot(xi,MG,'r--','linewidth',2); drawnow;
end
save('SteadyGROM.mat','romG');
load('SteadyGROM.mat');

%GNAT - Galerkin
romGnlfile = 'GalROM_nonlinBase';
NLSnap=1;
readResJacComputePODWriteFcn(romG,romGnlfile,NLSnap,[]);
gnatG = genGNAT([fname,'.rom'],romG,1);
for j=1:romG.nBases
    [phiR,phiJ] = readNonlinBases(romGnlfile,NLSnap,j);
    gnatG.createGNAT(prob(1),romG,phiR,phiJ,j);
end
if romG.nBases > 1, gnatG.precomputeDistanceQuantities(romG); end;
gnatG.executeModel;
gnatG.associateFullProblem(prob(1));
svG_g = gnatG.reconstructFullState(romG);
save('SteadyGGNAT.mat','gnatG');

uGNATg = romG.prob.velocity([prob.U0;svG_g(:,end);prob.UL]);
rhoGNATg = romG.prob.density(uGNATg);
MGNATg = romG.prob.mach(rhoGNATg,uGNATg);
plot(xi,MGNATg,'go');

errGgnat = 100*mean(ColumnwiseNorm(fom.sv - svG_g,2)./ColumnwiseNorm(fom.sv,2))

figure;
plot(xi,MF,'k-','linewidth',2); hold on;
plot(xi,MG,'m--','linewidth',2);
plot(xi,MGNATg,'g:','linewidth',2);
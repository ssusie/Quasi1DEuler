%This should be the filename of the CFG file
fname = 'example1000dof';
addGlobVarsAndPaths(fname);
%This defines the simulation with 1000 dofs.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FOM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This line will initialize the 1st configuration in the Burger.cfg file (it
%will initialize configuration with id = 1.
%The output is an array of FOM objects (stored in fom) and OneDBurger
%objects (stored in prob).  In these case, these arrays only have one
%entry.  To initialize all four configurations, replace 1 with 1:4
[fom,prob] = initialize(fname,1);
% [fom,prob] = initialize(fname,1:4);

%This will loop over all foms created and execute the full order model
%simulation.
for i = 1:length(fom)
    fom(i).executeModel;
end
save([fname,'_FOM.mat'],'fom','prob');

% X = fom.sv;
% [U1,S1,V1] = svd(bsxfun(@minus,X(:,2:end),fom.sv(:,1)),0);
% [U2,S2,V2] = svd(X(:,2:end)-X(:,1:end-1),0);
% plot(1-cumsum(diag(S1))/sum(diag(S1))+1e-16,'ko'); hold on; plot(1-cumsum(diag(S2))/sum(diag(S2))+1e-16,'b*'); set(gca,'yscale','log');
% grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ROM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This line will initialize a single ROM object (put inside a loop for
%more).  It will read from the Burger.rom file and it will use the first
%id in the PGROM-block (the 'pg' means it will only look for the
%PGROM-block, where PG means Petrov-Galerkin; the 1 means it use the first
%one specified).  The prob(1) means that the rom that is returned from this
%function will be linked to the first entry in the prob array.
%To change to Galerkin, change 'pg' to 'g'.
rom = ROM([fname,'.rom'],'pg-reg',1,prob(1));

%Set the training initial condition. It is trivial in the case, but for
%local ROB with snapshots coming from multiple simulations, this is
%important.  If this is set to empty, the function rom.determineSnapshots
%knows to use the initial condition of the problem that rom is linked to
%(see previous comment for mention of this linking).
trainIC = prob(1).ic;
%This will compute the snapshots based of the state vectors of fom
%and the training initial condition.  This will involve (possibly)
%subtracting the initial condition from the state vectors OR for each
%state vector subtracting the state vector at the previous time step.
%This also involve choosing the times steps where we will take
%snapshots and the number of snapshots to take (specified in nsnaps,
%snapInt, and snapDist in the .rom file).  This function also clusters the
%snapshots (for the local case).
rom.clusterSnapsOnly(fom,trainIC);
%Compute the POD for each cluster
rom.computePOD('fromClustSnaps');
% rom.computePOD('fromClustSnaps',[],[],true,fom);

%Apriori error computations (projection error and residual ''bound'')
if rom.nBases>1
    rom.visualizeClustering();
end
rom.computeResNormBound(fom,'rel',2);
rB = rom.computeResNormBoundAtTimeVsBasisSize(500,fom,'rel',2);
rom.computeProjError(fom,true,2);

%This will execute the ROM online simulation.
rom.executeModel;

%Compute residual of rom solution at all times
rom.computeResNormFromSV('rel',2);

% return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%GNAT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%These lines will perform the pod of the residual and jacobian snapshots.
%Since no arguments were passed to POD, no truncation will take place.
%This is recommended by Kevin for computing sample nodes.  Truncation will
%take place inside createGNAT (although phiR, phiJ will not be affected
%because they are not stored anywhere) before computation of the online
%matrices.  You are free to truncate at this point if you wish.  Simply
%pass pod a second argument, which is the number of columns you want phiR
%or phiJ to have.  Note that if we are performing a local ROB analysis,
%rom.res and rom.jac will be a cell array of size nBases x 1, but pod knows
%how to handle this case, so there is no need to worry.
fnameNLbase='romNonlinBase';
NLSnapshot=2;
readResJacComputePODWriteFcn(rom,fnameNLbase,NLSnapshot,[]);

%This will initialize a single GNAT object (put inside loop if you want
%more).  It will read from the first id number (not id = 1, but the first
%id number, whatever it happends to be) (this is the meaning of the 1 in
%the 3rd argument) GNAT-block in the Burger.rom file. Furthermore, it will be 
%linked to the ROM object rom (i.e. the gnat will be an approximation built
%on rom - it will inherit nY and nBases from rom).
gnat = genGNAT([fname,'.rom'],rom,1);

%This will perform all of the offline GNAT computations.  Computation of
%sample nodes (or indices) and generation of online matrices (plus many
%more things that are not important from a user's standpoint).  phiR,
%phiJ, and rom.phi are used in these computations, but they are not stored
%(only rom.phi evaluated at the state mask is stored).  This also
%associates prob(1) with the gnat object, but it is not stored.  A small
%version of prob(1) is stored (internally named probGNAT).  The idea is
%that this will reduce the amount of memory being used during the
%computation and speed up the computation.
for i = 1:rom.nBases
    [phiR,phiJ] = readNonlinBases(fnameNLbase,1,i);
    gnat.createGNAT(prob(1),rom,phiR,phiJ,i);
end
if rom.nBases > 1
    gnat.precomputeDistanceQuantities(rom);
end

%This will perform the online simulation of the gnat object.    It is not
%done here, but the can clear all variables except gnat at this point (can
%be accomplished with the command clearAllExcept('gnat');). 
%This will reduce the computational resources needed and speed up the
%computation.  It is advised to save fom, prob, and rom to disk before
%doing this.  Otherwise you will have to run the code before this point
%again before postprocessing.
gnat.executeModel;

%This will reassociate the full scale problem with the gnat object.  This
%is necessary for postprocessing.  Now that the online simulation is
%complete, we are no concerned with storage optimization (as long as we
%don't need too much).
gnat.associateFullProblem(prob(1));

svG = gnat.reconstructFullState(rom);

RomErr = mean(ColumnwiseNorm(rom.sv-fom.sv,2)./ColumnwiseNorm(fom.sv,2));
fprintf('Average Relative L2 ROM Error = %f %%\n',RomErr);
GnatErr = mean(ColumnwiseNorm(svG-fom.sv,2)./ColumnwiseNorm(fom.sv,2));
fprintf('Average Relative L2 GNAT Error = %f %%\n',GnatErr);

figure;
for i = [100:100:1000]
    hfom=plot(fom.sv(:,i),'k','linewidth',2); hold on;
    hrom=plot(rom.sv(:,i),'b--','linewidth',2);
    hgnat=plot(svG(:,i),'g'); hold on;
end
xlabel('x')
ylabel('U')
legend([hfom,hrom,hgnat],'FOM','ROM','GNAT')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;
fname = 'eulerExample300dof';
addGlobVarsAndPaths(fname);
METHODS = {'FOM', 'Method1\_g', 'Method1\_pg', 'Method2', 'Method3',...
    'Method4', 'Method5\_ROM', 'Method5\_GNAT', 'Method5\_FOM'};
%FOM
[fom,prob] = initialize(fname,10);
fom.executeModel;
soln=fom.sv;
SOLN = fom.sv;
%% 
% Method 1; original ROM with 'g' or 'pg'
Gal_PetGal = {'g','pg'};

methodROM=1; %solve with constraints  % =1  original,  =2 constrained
numCell=1; %number of domains for constraints
basisNumber=10; %truncate basis Phi; if don't want to truncate choose basisNumber=fom.cTimeIter
rom1 = ROM([fname,'.rom'], Gal_PetGal{1}, 1, prob, methodROM, numCell, basisNumber);
trainIC = prob.ic;
rom1.clusterSnapsOnly(fom,trainIC);
rom1.computePOD('fromClustSnaps');
rom1.executeModel; 
SOLN = [SOLN, rom1.sv];
CONSTR = rom1.Cnorm';
%%
rom2 = ROM([fname,'.rom'], Gal_PetGal{2}, 1, prob, methodROM, numCell, basisNumber);
rom2.clusterSnapsOnly(fom,trainIC);
rom2.computePOD('fromClustSnaps');
rom2.executeModel; 
SOLN = [SOLN, rom2.sv];
CONSTR = [CONSTR, rom2.Cnorm'];

%%
% Method 2; Rom with constraints
methodROM=2; %solve with constraints  % =1  original,  =2 constrained
numCell=1; %number of domains for constraints
basisNumber=10; %truncate basis Phi; if don't want to truncate choose basisNumber=fom.cTimeIter
rom3 = ROM([fname,'.rom'],'pg', 1, prob, methodROM, numCell, basisNumber);
trainIC = prob.ic;
rom3.clusterSnapsOnly(fom,trainIC);
rom3.computePOD('fromClustSnaps');
rom3.executeModel;
SOLN = [SOLN, rom3.sv];
CONSTR = [CONSTR, rom3.Cnorm'];

%% 
% Method 3; GNAT original (no constriants) using ROM 'pg'

fnameNLbase='NonlinBase';
NLSnapshot = 0;
methodGNAT = 1; % = 1 original, = 2 Rom constrains; = 3 Gnat constraints
augment = false;
readResJacComputePODWriteFcn(fom,fnameNLbase,NLSnapshot,[]);

gnat0 = genGNAT([fname,'.rom'], rom1, 1); % rom1 ~ ROM 'pg'
[phiR,phiJ] = readNonlinBases(fnameNLbase,NLSnapshot,1);
gnat0.createGNAT(prob(1), rom1, phiR, phiJ, methodGNAT, augment, 1);
gnat0.executeModel;

gnat0.associateFullProblem(prob(1));
svG0 = gnat0.reconstructFullState(rom1);
SOLN = [SOLN, svG0];
CONSTR = [CONSTR, gnat0.Rnorm'];

%%
% Method 4; GNAT with real constraints (from ROM)
 
fnameNLbase='NonlinBase';
NLSnapshot = 0;
methodGNAT = 2; % = 1 original, = 2 Rom constrains; = 3 Gnat constraints
augment = false;

readResJacComputePODWriteFcn(fom,fnameNLbase,NLSnapshot,[]);
gnat1 = genGNAT([fname,'.rom'], rom3, 1); % ROM with constriants

[phiR,phiJ] = readNonlinBases(fnameNLbase,NLSnapshot,1);
gnat1.createGNAT(prob(1), rom3, phiR, phiJ, methodGNAT, augment, 1);
gnat1.executeModel;
gnat1.associateFullProblem(prob(1));
svG1 = gnat1.reconstructFullState(rom3);
SOLN = [SOLN, svG1];
CONSTR = [CONSTR, gnat1.Rnorm'];
Approx_CONSTR = [gnat1.Anorm_rom', gnat1.Anorm_fom'];
%%
% Method 5; GNAT with approx constraints using snapshots from ROM

rom3.buildFluxBasis(rom3.sv); % build flux basis from ROM snapshots.

fnameNLbase='NonlinBase';
NLSnapshot = 0;
methodGNAT = 3; % = 1 original, = 2 Rom constrains; = 3 Gnat constraints
augment = false;

readResJacComputePODWriteFcn(fom,fnameNLbase,NLSnapshot,[]);
gnat2 = genGNAT([fname,'.rom'],rom3,1);
[phiR,phiJ] = readNonlinBases(fnameNLbase,NLSnapshot,1);
gnat2.createGNAT(prob(1),rom3,phiR,phiJ,methodGNAT, augment, 1);
gnat2.executeModel;

gnat2.associateFullProblem(prob(1));
svG2 = gnat2.reconstructFullState(rom3);
SOLN = [SOLN, svG2];
CONSTR = [CONSTR, gnat2.Anorm'];
Real_CONSTR = gnat2.Rnorm';

%%
% Method 5; GNAT with approx constraints using snapshots from GNAT

fnameNLbase='NonlinBase';
NLSnapshot = 0;
methodGNAT = 4; % = 1 original, = 2 Rom constrains; = 3 Gnat constraints

readResJacComputePODWriteFcn(fom,fnameNLbase,NLSnapshot,[]);
gnat3 = genGNAT([fname,'.rom'], rom3, 1);
[phiR,phiJ] = readNonlinBases(fnameNLbase,NLSnapshot,1);
gnat3.createGNAT(prob(1), rom3, phiR, phiJ, methodGNAT, augment, 1);
gnat3.executeModel;

gnat3.associateFullProblem(prob(1));
svG3 = gnat3.reconstructFullState(rom3);
SOLN = [SOLN, svG3];
CONSTR = [CONSTR, gnat3.Anorm'];
Real_CONSTR = [Real_CONSTR ,gnat3.Rnorm'];
%%
% Method 5; GNAT with approx constraints using snapshots from FOM
rom3.buildFluxBasis(fom.sv);

fnameNLbase='NonlinBase';
NLSnapshot = 0;
methodGNAT = 5; % = 1 original, = 2 Rom constrains; = 3 Gnat constraints
augment = false;

readResJacComputePODWriteFcn(fom,fnameNLbase,NLSnapshot,[]);
gnat4 = genGNAT([fname,'.rom'],rom3,1);
[phiR,phiJ] = readNonlinBases(fnameNLbase,NLSnapshot,1);
gnat4.createGNAT(prob(1),rom3,phiR,phiJ,methodGNAT,  augment, 1);
gnat4.executeModel;

gnat4.associateFullProblem(prob(1));
svG4 = gnat4.reconstructFullState(rom3);
SOLN = [SOLN, svG4];
CONSTR = [CONSTR, gnat4.Anorm'];
Real_CONSTR = [Real_CONSTR ,gnat4.Rnorm'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

figure(1)
N = fom.cTimeIter+1;
thick = linspace(1,3, length(METHODS) - 1);
thick = sort(thick, 'descend');

for i = 1 : length(METHODS) - 1
    col = rand(1,3);
    figure(1)
    if i <=3
        mark = 'o';
        style = '-';
    else
        mark = '*';
        style = '--';
    end
    
    semilogy(CONSTR(:,i),'LineStyle', style, 'Marker', mark,...
        'color', col, 'linewidth', thick(i)), hold on
    xlabel('time step')
    ylabel('log norm(method specific constraints)')
    set(gca, 'Xtick', 1:N-1)
end
list = { 'Method1\_g', 'Method1\_pg', 'Method2', 'Method3', ...
    'Method4', 'Method5\_ROM', 'Method5\_GNAT', 'Method5\_FOM'};
legend(list)

%%
METHODS = {'FOM', 'Method1\_g', 'Method1\_pg', 'Method2', 'Method3',...
    'Method4', 'Method5\_ROM', 'Method5\_GNAT', 'Method5\_FOM'};

for i = 2 : length(METHODS)
    col = rand(1,3);
%     keyboard
    figure(2)
    z1 = abs(SOLN(1:3:end,N*i)- SOLN(1:3:end,N))./norm(SOLN(1:3:end,N),2);
    semilogy(z1,'-', 'color', col, 'linewidth', thick(i-1)), hold on
    title('rho')
    xlabel('cell number spatial domain')
    ylabel('log of Rel error at final time step for rho')

%     keyboard
    figure(3)
    z2 = abs(SOLN(2:3:end,N*i )- SOLN(2:3:end,N))./norm(SOLN(2:3:end,N),2);
    semilogy(z2,'-', 'color', col, 'linewidth', thick(i-1)), hold on
    title('rho\_u')
    xlabel('cell number spatial domain')
    ylabel('log of Rel error at final time step for rho u')

%     keyboard
    figure(4);
    z3 = abs(SOLN(3:3:end,N*i )- SOLN(3:3:end,N))./norm(SOLN(3:3:end,N),2);
    semilogy(z3,'-', 'color', col, 'linewidth', thick(i-1)); hold on;
    xlabel('cell number spatial domain')
    ylabel('log of Rel error at final time step for e')
    title('e')

%     keyboard
    if i < 5
        mark = 'o';
        style = '-';
    else
        mark = '+';
        style = '--';
    end
    figure(5);
    E =  ColumnwiseNorm(SOLN(:, N * (i - 1) + 1 : N * i )- SOLN(:,1:N)) ./ ColumnwiseNorm(SOLN(:,1:N),2);
    semilogy(E, 'LineStyle', style, 'Marker', mark, 'color', col, 'linewidth', thick(i-1)); 
    hold on % 2-norm of columns

    xlabel('time step')
    ylabel('log of Rel error at each time step')
    set(gca, 'Xtick', 1:N-1)
end

for i=2:5,  figure(i);  legend(METHODS{2:end});  end;


%%
figure(6)
for i = 1: length(METHODS)
    col = rand(1,3);
    [~,uF,~,cF,~] = prob.getVariables(SOLN(:,N*i));
    plot(uF./cF,'color', col, 'linewidth', 2); hold on;
    xlabel('spatial domain')
    ylabel('Mach')
end

legend(METHODS)

%%%%%%%%%%%%%
%%
for i = 2 : length(METHODS)
    col = rand(1,3);

    figure(7)
    z1 = abs(SOLN(1:3:end,N*i)- SOLN(1:3:end,N))./norm(SOLN(1:3:end,N),2);
    plot(z1,'-', 'color', col, 'linewidth', thick(i-1)), hold on
    title('rho')
    xlabel('cell number spatial domain')
    ylabel('Rel error at final time step for rho')

    figure(8)
    z2 = abs(SOLN(2:3:end,N*i )- SOLN(2:3:end,N))./norm(SOLN(2:3:end,N),2);
    plot(z2,'-', 'color', col, 'linewidth', thick(i-1)), hold on
    title('rho\_u')
    xlabel('cell number spatial domain')
    ylabel('Rel error at final time step for rho u')

    figure(9);
    z3 = abs(SOLN(3:3:end,N*i )- SOLN(3:3:end,N))./norm(SOLN(3:3:end,N),2);
    plot(z3,'-', 'color', col, 'linewidth', thick(i-1)); hold on;
    xlabel('cell number spatial domain')
    ylabel('Rel error at final time step for e')
    title('e')

    if i < 5
        mark = 'o';
        style = '-';
    else
        mark = '+';
        style = '--';
    end
  
end
for i=7:9, figure(i); legend(METHODS{2:end}); end;

%%
figure(10)
thick2 = [3,2.5,2];
for i = 1 : 3
    col = rand(1,3);
    
    semilogy(Real_CONSTR(:,i),'-o','color', col, 'linewidth', thick2(i)), 
    hold on
%     semilogy(CONSTR(:,4+i),'--+','color', col, 'linewidth', thick2(i)), 
    xlabel('time step')
    ylabel('log norm( REAL constr)')
end
set(gca, 'Xtick', 1:N-1)
% list = { 'Method5\_rom\_real', 'Method5\_rom', 'Method5\_GNAT\_real' ...
%     'Method5\_GNAT', 'Method5\_fom\_real', 'Method5\_fom'};
list = {  'Method5\_rom', 'Method5\_GNAT', 'Method5\_fom'};
legend(list)

%%
figure(11)
semilogy(CONSTR(:,5),'--+','color', 'r', 'linewidth', thick2(2)), hold on
semilogy(Approx_CONSTR(:,1),'-o','color', 'g', 'linewidth', thick2(2)), 
semilogy(Approx_CONSTR(:,2),'--v','color', 'b', 'linewidth', thick2(2)), 
xlabel('time step')
ylabel('log norm (constr)')
legend('Method4', 'Method4\_approx\_rom', 'Method4\_approx\_fom')



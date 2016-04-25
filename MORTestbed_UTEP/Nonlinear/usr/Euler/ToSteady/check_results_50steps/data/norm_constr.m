clear all; close all; clc;

G1=load('GconstrOrig');
G2=load('GconstrOrig_trunc4');
G3=load('GconstrOrig_trunc3');

PG1=load('PGconstrOrig');
PG2=load('PGconstrOrig_trunc4');
PG3=load('PGconstrOrig_trunc3');

Z1=load('constrZim');
Z2=load('constrZim_trunc4');
Z3=load('constrZim_trunc3');

%%
normZ_C=norm(Z1.C,'fro')
normZ4=norm(Z2.Czim_trunc4,'fro')
normZ3=norm(Z3.Czim_trunc3,'fro')

normPG=norm(PG1.Corig,'fro')
normPG4=norm(PG2.Corig_trunc4,'fro')
normPG3=norm(PG3.Corig_trunc3,'fro')

normG=norm(G1.Corig,'fro')
normG4=norm(G2.Corig_trunc4,'fro')
normG3=norm(G3.Corig_trunc3,'fro')

%%

for i=1:size(Z1.C,2)
    z(i)=norm(Z1.C,2);    
end

for i=1:size(PG1.Corig,2)
    pg(i)=norm(PG1.Corig,2);    
end

for i=1:size(G1.Corig,2)
    g(i)=norm(G1.Corig,2);    
end
figure(1)
plot(log(1+z),'r*', 'linewidth',2); hold on
plot(log(1+diag(sqrt(Z1.C'*Z1.C))), 'g*')
plot(log(1+pg),'bo', 'linewidth',2)
plot(log(1+g),'kv', 'linewidth',2)
legend('Constrained', 'PG', 'G')




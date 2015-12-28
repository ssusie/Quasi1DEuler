clear all; close all; clc;

G1=load('GconstrOrig');
G2=load('GconstrOrigTrunc4');
G3=load('GconstrOrigTrunc3');

PG1=load('PGconstrOrig');
PG2=load('PGconstrOrigTrunc4');
PG3=load('PGconstrOrigTrunc3');

Z1=load('constrZim');
Z2=load('constrZim_trunc4');
Z3=load('constrZim_trunc3');

%%
normZ_C=norm(Z1.Czim,'fro')
normZ4=norm(Z2.Czim,'fro')
normZ3=norm(Z3.Czim,'fro')

normPG=norm(PG1.Corig,'fro')
normPG4=norm(PG2.Corig_trunc4,'fro')
normPG3=norm(PG3.Corig_trunc3,'fro')

normG=norm(G1.Corig,'fro')
normG4=norm(G2.Corig_trunc4,'fro')
normG3=norm(G3.Corig_trunc3,'fro')

%%

for i=1:size(Z1.Czim,2)
    z(i)=norm(Z1.Czim(:,i),2);    
    pg(i)=norm(PG1.Corig(:,i),2);    
    g(i)=norm(G1.Corig(:,i),2);

    z3(i)=norm(Z3.Czim(:,i),2);    
    pg3(i)=norm(PG3.Corig_trunc3(:,i),2);    
    g3(i)=norm(G3.Corig_trunc3(:,i),2);

    z4(i)=norm(Z2.Czim(:,i),2);    
    pg4(i)=norm(PG2.Corig_trunc4(:,i),2);    
    g4(i)=norm(G2.Corig_trunc4(:,i),2);

end

figure(1)
plot(z,'r*', 'linewidth',2); hold on
plot(pg,'bo', 'linewidth',2)
plot(g,'kv', 'linewidth',2)
legend('Constrained', 'PG', 'G')
title('norm of constraints on each timestep with 5 basis vectors')

figure(2)
plot(z4,'r*', 'linewidth',2); hold on
plot(pg4,'bo', 'linewidth',2)
plot(g4,'kv', 'linewidth',2)
legend('Constrained', 'PG', 'G')
title('norm of constraints on each timestep with 4 basis vectors')

figure(3)
plot(z3,'r*', 'linewidth',2); hold on
plot(pg3,'bo', 'linewidth',2)
plot(g3,'kv', 'linewidth',2)
legend('Constrained', 'PG', 'G')
title('norm of constraints on each timestep with 3 basis vectors')



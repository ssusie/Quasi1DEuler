clear all; close all; clc;

G1=load('Gconstr5');
G2=load('Gconstr999');
G3=load('Gconstr9999');

PG1=load('PGconstr5');
PG2=load('PGconstr999');
PG3=load('PGconstr9999');

Z1=load('constrZim5');
Z2=load('constrZim999');
Z3=load('constrZim9999');

%%
normZ_C=norm(Z1.Czim,'fro')
normZ4=norm(Z2.Czim,'fro')
normZ3=norm(Z3.Czim,'fro')

normPG=norm(PG1.Corig,'fro')
normPG4=norm(PG2.Corig,'fro')
normPG3=norm(PG3.Corig,'fro')

normG=norm(G1.Corig,'fro')
normG4=norm(G2.Corig,'fro')
normG3=norm(G3.Corig,'fro')

%%

for i=1:size(Z1.Czim,2)
    z(i)=norm(Z1.Czim(:,i),2);    
    pg(i)=norm(PG1.Corig(:,i),2);    
    g(i)=norm(G1.Corig(:,i),2);

    z3(i)=norm(Z3.Czim(:,i),2);    
    pg3(i)=norm(PG3.Corig(:,i),2);    
    g3(i)=norm(G3.Corig(:,i),2);

    z4(i)=norm(Z2.Czim(:,i),2);    
    pg4(i)=norm(PG2.Corig(:,i),2);    
%     g4(i)=norm(G2.Corig(:,i),2);

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
% plot(g4,'kv', 'linewidth',2)
legend('Constrained', 'PG')
title('norm of constraints 99.9%')

figure(3)
plot(z3,'r*', 'linewidth',2); hold on
plot(pg3,'bo', 'linewidth',2)
plot(g3,'kv', 'linewidth',2)
legend('Constrained', 'PG', 'G')
title('norm of constraints 99.99%')



clear all; close all; clc;
original=load('originalROM');
constrained=load('constrainedROM');


[rhoO,uO,PO,cO,eO] = prob.getVariables(original.rom.sv(:,end));
[rhoC,uC,PC,cC,eC] = prob.getVariables(constrained.rom.sv(:,end));

figure;
horig  = plot(uF./cF,'k','linewidth',2); hold on;
hconstr  = plot(uR./cR,'b--','linewidth',2);
xlabel('x')
ylabel('Mach')
legend([hfom,hrom,hgnat],'Original','Constrained')









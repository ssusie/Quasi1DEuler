function plotAll(fom, rom1, rom2, rom3,...
                  gnat0, gnat1, gnat2, gnat3, gnat4, gnat5, gnat6,...
                  svG0, svG1, svG2, svG3, svG4, svG5, svG6, prob)

close all;
N = fom.cTimeIter+1;
%METHODS = {'FOM','G','PG','PGcnstd','GNAT','GNATcnstd','GNATcnstd(PGcnstd)',...
%    'GNATcnstd(GNATcnstd)','GNATcnstd(FOM)','GNATcnstd(GNAT)'};
%list = { 'G','PG','PGcnstd','GNAT','GNATcnstd','GNATcnstd(PGcnstd)',...
%    'GNATcnstd(GNATcnstd)','GNATcnstd(FOM)','GNATcnstd(GNAT)'};
 METHODS = {'FOM','G','PG','PGcnstd','GNAT','GNATcnstd','GNATcnstd(PGcnstd)',...
     'GNATcnstd(GNATcnstd)','GNATcnstd(FOM)','GNATcnstd(GNAT)','GNATcnstd(PG)'};
 list = { 'G','PG','PGcnstd','GNAT','GNATcnstd','GNATcnstd(PGcnstd)',...
     'GNATcnstd(GNATcnstd)','GNATcnstd(FOM)','GNATcnstd(GNAT)','GNATcnstd(PG)'};
% METHODS = {'FOM','G','PG','PGcnstd','GNAT','GNATcnstd','GNATcnstd(PGcnstd)',...
%     'GNATcnstd(GNATcnstd)','GNATcnstd(FOM)','GNATcnstd(PG)'};
% list = { 'G','PG','PGcnstd','GNAT','GNATcnstd','GNATcnstd(PGcnstd)',...
%     'GNATcnstd(GNATcnstd)','GNATcnstd(FOM)','GNATcnstd(PG)'};

%% 
figure(1)
hfom = semilogy(fom.Cnorm,'k-*','linewidth',5,'markersize',13); hold on;
hrom1 = semilogy(rom1.Cnorm,'LineStyle','-','Marker','o','color','b','linewidth',5,'markersize',10);
hrom2 = semilogy(rom2.Cnorm,'LineStyle','-','Marker','x','color','g','linewidth',5,'markersize',9);
hrom3 = semilogy(rom3.Cnorm,'LineStyle','-','Marker','v','color','c','linewidth',5,'markersize',10);
hgnat0 = semilogy(gnat0.Rnorm,'LineStyle','-','Marker','o','color',[0.7,0.7,0.9],'linewidth',4,'markersize',4);
hgnat1 = semilogy(gnat1.Rnorm,'LineStyle','-','Marker','o','color',[0.5,0.2,0.1],'linewidth',4,'markersize',8);
hgnat2 = semilogy(gnat2.Anorm,'LineStyle','--','Marker','o','color',[0.1,0.8,0.2],'linewidth',3.5,'markersize',6);
hgnat3 = semilogy(gnat3.Anorm,'LineStyle','-','Marker','o','color','m','linewidth',3,'markersize',4);
hgnat4 = semilogy(gnat4.Anorm,'LineStyle','--','Marker','o','color',[0.1,0.1,0.1],'linewidth',2.5,'markersize',3);
hgnat5 = semilogy(gnat5.Anorm,'LineStyle','-','Marker','o','color',[0.1,0.55,0.1],'linewidth',2,'markersize',2);
hgnat6 = semilogy(gnat6.Anorm,'LineStyle','-','Marker','o','color',[0.1,0.55,1],'linewidth',1,'markersize',1);
legend([hfom,hrom1,hrom2,hrom3,hgnat0,hgnat1,hgnat2,hgnat3,hgnat4,hgnat5,hgnat6],METHODS);

xlabel('time step')
ylabel('norm of method specific constraints')

%%
figure(5)
Erom1 =  ColumnwiseNorm(rom1.sv-fom.sv)./ColumnwiseNorm(fom.sv,2);
Erom2 =  ColumnwiseNorm(rom2.sv-fom.sv)./ColumnwiseNorm(fom.sv,2);
Erom3 =  ColumnwiseNorm(rom3.sv-fom.sv)./ColumnwiseNorm(fom.sv,2);
Egnat0 = ColumnwiseNorm(svG0-fom.sv)./ColumnwiseNorm(fom.sv,2);
Egnat1 = ColumnwiseNorm(svG1-fom.sv)./ColumnwiseNorm(fom.sv,2);
Egnat2 = ColumnwiseNorm(svG2-fom.sv)./ColumnwiseNorm(fom.sv,2);
Egnat3 = ColumnwiseNorm(svG3-fom.sv)./ColumnwiseNorm(fom.sv,2);
Egnat4 = ColumnwiseNorm(svG4-fom.sv)./ColumnwiseNorm(fom.sv,2);
Egnat5 = ColumnwiseNorm(svG5-fom.sv)./ColumnwiseNorm(fom.sv,2);
Egnat6 = ColumnwiseNorm(svG6-fom.sv)./ColumnwiseNorm(fom.sv,2);
hrom1 = semilogy(Erom1, 'LineStyle','-','Marker','o','color','b','linewidth',5,'markersize',10); hold on;
hrom2 = semilogy(Erom2, 'LineStyle','-','Marker','x','color','g','linewidth',5,'markersize',9);
hrom3 = semilogy(Erom3, 'LineStyle','-','Marker','v','color','c','linewidth',5,'markersize',10);
hgnat0 = semilogy(Egnat0, 'LineStyle','-','Marker','o','color',[0.7,0.7,0.9],'linewidth',4,'markersize',4);
hgnat1 = semilogy(Egnat1, 'LineStyle','-','Marker','o','color',[0.5,0.2,0.1],'linewidth',4,'markersize',8);
hgnat2 = semilogy(Egnat2, 'LineStyle','--','Marker','o','color',[0.1,0.8,0.2],'linewidth',3.5,'markersize',6);
hgnat3 = semilogy(Egnat3, 'LineStyle','-','Marker','o','color','m','linewidth',3,'markersize',4);
hgnat4 = semilogy(Egnat4, 'LineStyle','--','Marker','o','color',[0.1,0.1,0.1],'linewidth',2.5,'markersize',3);
hgnat5 = semilogy(Egnat5, 'LineStyle','-','Marker','o','color',[0.1,0.55,0.1],'linewidth',2,'markersize',2);
hgnat6 = semilogy(Egnat6, 'LineStyle','-','Marker','o','color',[0.1,0.55,1],'linewidth',1,'markersize',1);
xlabel('time step')
ylabel('rel error with FOM solution')
legend([hrom1,hrom2,hrom3,hgnat0,hgnat1,hgnat2,hgnat3,hgnat4,hgnat5,hgnat6],list);
% legend([hrom1,hrom2,hrom3,hgnat0,hgnat1,hgnat2,hgnat3,hgnat4,hgnat5],list);

%%
figure(6)
loops = size(fom.sv,2);
F(loops) = struct('cdata',[],'colormap',[]);
gap = 1;
for i=[loops]
    [rhoF,uF,PF,cF,eF] = prob.getVariables(fom.sv(:,i));
    [rhoR1,uR1,PR1,cR1,eR1] = prob.getVariables(rom1.sv(:,i));
    [rhoR2,uR2,PR2,cR2,eR2] = prob.getVariables(rom2.sv(:,i));
    [rhoR3,uR3,PR3,cR3,eR3] = prob.getVariables(rom3.sv(:,i));
    [rhoGN0,uGN0,PGN0,cGN0,eGN0] = prob.getVariables(real(svG0(:,i)));
    [rhoGN1,uGN1,PGN1,cGN1,eGN1] = prob.getVariables(svG1(:,i));
    [rhoGN2,uGN2,PGN2,cGN2,eGN2] = prob.getVariables(svG2(:,i));
    [rhoGN3,uGN3,PGN3,cGN3,eGN3] = prob.getVariables(svG3(:,i));
    [rhoGN4,uGN4,PGN4,cGN4,eGN4] = prob.getVariables(svG4(:,i));    
    [rhoGN5,uGN5,PGN5,cGN5,eGN5] = prob.getVariables(svG5(:,i));  
    [rhoGN6,uGN6,PGN6,cGN6,eGN6] = prob.getVariables(svG6(:,i));  
    hfom  = plot(uF(1:gap:end)./cF(1:gap:end),'k-*','linewidth',5,'markersize',13); hold on;
    hrom1  = plot(uR1(1:gap:end)./cR1(1:gap:end),'b-o','linewidth',5,'markersize',10);
    hrom2  = plot(real(uR2(1:gap:end)./cR2(1:gap:end)),'g-x','linewidth',5,'markersize',9);
    hrom3  = plot(uR3(1:gap:end)./cR3(1:gap:end),'c-v','linewidth',5,'markersize',10);
    hgnat0 = plot(real(uGN0(1:gap:end)./cGN0(1:gap:end)),'-o','linewidth',4,'markersize',4,'color',[0.7,0.7,0.9]);
    hgnat1 = plot(uGN1(1:gap:end)./cGN1(1:gap:end),'-o','linewidth',4,'markersize',8,'color',[0.5,0.2,0.1]);
    hgnat2 = plot(uGN2(1:gap:end)./cGN2(1:gap:end),'--o','linewidth',3.5,'markersize',6,'color',[0.1,0.8,0.2]);
    hgnat3 = plot(uGN3(1:gap:end)./cGN3(1:gap:end),'m-o','linewidth',3,'markersize',4);
    hgnat4 = plot(uGN4(1:gap:end)./cGN4(1:gap:end),'--o','linewidth',2.5,'markersize',3,'color',[0.1,0.1,0.1]);
    hgnat5 = plot(uGN5(1:gap:end)./cGN5(1:gap:end),'-o','linewidth',2,'markersize',2,'color',[0.1,0.55,0.1]);
    hgnat6 = plot(real(uGN6(1:gap:end)./cGN6(1:gap:end)),'-o','linewidth',1,'markersize',1,'color',[0.1,0.55,1]);
    legend([hfom,hrom1,hrom2,hrom3,hgnat0,hgnat1,hgnat2,hgnat3,hgnat4,hgnat5,hgnat6],METHODS);
%     legend([hfom,hrom1,hrom2,hrom3,hgnat0,hgnat1,hgnat2,hgnat3,hgnat4,hgnat5],METHODS);
    xlabel('spatial domain')
    ylabel('Mach')
    ylim([0 5]);
end

%% Movie (either comment or do not comment)   
% loops = size(fom.sv,2);
% figure(7)
% F(loops) = struct('cdata',[],'colormap',[]);
% gap = 2;
% for i=1:1:loops
%     [rhoF,uF,PF,cF,eF] = prob.getVariables(fom.sv(:,i));
%     [rhoR1,uR1,PR1,cR1,eR1] = prob.getVariables(rom1.sv(:,i));
%     [rhoR2,uR2,PR2,cR2,eR2] = prob.getVariables(rom2.sv(:,i));
%     [rhoR3,uR3,PR3,cR3,eR3] = prob.getVariables(rom3.sv(:,i));
%     [rhoGN0,uGN0,PGN0,cGN0,eGN0] = prob.getVariables(svG0(:,i));
%     [rhoGN1,uGN1,PGN1,cGN1,eGN1] = prob.getVariables(svG1(:,i));
%     [rhoGN2,uGN2,PGN2,cGN2,eGN2] = prob.getVariables(svG2(:,i));
%     [rhoGN3,uGN3,PGN3,cGN3,eGN3] = prob.getVariables(svG3(:,i));
%     [rhoGN4,uGN4,PGN4,cGN4,eGN4] = prob.getVariables(svG4(:,i));  
%     [rhoGN5,uGN5,PGN5,cGN5,eGN5] = prob.getVariables(svG5(:,i)); 
%     [rhoGN6,uGN6,PGN6,cGN6,eGN6] = prob.getVariables(svG6(:,i));
%     hfom  = plot(uF(1:gap:end)./cF(1:gap:end),'k-*','linewidth',2,'markersize',13); hold on;
%     hrom1  = plot(uR1(1:gap:end)./cR1(1:gap:end),'b-o','linewidth',2,'markersize',10);
%     hrom2  = plot(uR2(1:gap:end)./cR2(1:gap:end),'g-x','linewidth',2,'markersize',7);
%     hrom3  = plot(uR3(1:gap:end)./cR3(1:gap:end),'c-v','linewidth',2,'markersize',7);
%     hgnat0 = plot(uGN0(1:gap:end)./cGN0(1:gap:end),'-.','linewidth',2,'markersize',7,'color',[0.7,0.7,0.9]);
%     hgnat1 = plot(uGN1(1:gap:end)./cGN1(1:gap:end),'-.','linewidth',2,'markersize',7,'color',[0.5,0.2,0.1]);
%     hgnat2 = plot(uGN2(1:gap:end)./cGN2(1:gap:end),'--v','linewidth',2,'markersize',7,'color',[0.1,0.8,0.2]);
%     hgnat3 = plot(uGN3(1:gap:end)./cGN3(1:gap:end),'m-.','linewidth',2,'markersize',3);
%     hgnat4 = plot(uGN4(1:gap:end)./cGN4(1:gap:end),'--.','linewidth',2,'markersize',7,'color',[0.1,0.1,0.1]);
%     hgnat5 = plot(uGN5(1:gap:end)./cGN5(1:gap:end),'-*','linewidth',2,'markersize',7,'color',[0.1,0.55,0.1]);
%     hgnat6 = plot(uGN6(1:gap:end)./cGN6(1:gap:end),'-+','linewidth',2,'markersize',7,'color',[0.1,0.55,1]);
%     legend([hfom,hrom1,hrom2,hrom3,hgnat0,hgnat1,hgnat2,hgnat3,hgnat4,hgnat5,hgnat6],METHODS);
% %     axis([0 length(uF(:,i)) -0.5 3])
%     xlabel('spatial domain')
%     ylabel('Mach')
%     drawnow
%     F(i) = getframe(gcf);
%     hold off;
% end


%%
N = fom.cTimeIter+1;
figure(10)
hgnat2 = semilogy(gnat2.Rnorm,'LineStyle','--','Marker','o','color',[0.1,0.8,0.2],'linewidth',3.5,'markersize',6); hold on;
hgnat3 = semilogy(gnat3.Rnorm,'LineStyle','-','Marker','o','color','m','linewidth',3,'markersize',4);
hgnat4 = semilogy(gnat4.Rnorm,'LineStyle','--','Marker','o','color',[0.1,0.1,0.1],'linewidth',2.5,'markersize',3);
hgnat5 = semilogy(gnat5.Rnorm,'LineStyle','-','Marker','o','color',[0.1,0.55,0.1],'linewidth',2,'markersize',2);
hgnat6 = semilogy(gnat5.Rnorm,'LineStyle','-','Marker','o','color',[0.1,0.55,1],'linewidth',1,'markersize',1);
% list = {'GNATcnstd(PGcnstd)', 'GNATcnstd(GNATcnstd)','GNATcnstd(FOM)','GNATcnstd(GNAT)'};
list = {'GNATcnstd(PGcnstd)', 'GNATcnstd(GNATcnstd)','GNATcnstd(FOM)','GNATcnstd(GNAT)','GNATcnstd(PG)'};
legend([hgnat2,hgnat3,hgnat4,hgnat5,hgnat6],list);
% legend([hgnat2,hgnat3,hgnat4,hgnat5],list);
xlabel('time step')
ylabel('norm of REAL constr')

end


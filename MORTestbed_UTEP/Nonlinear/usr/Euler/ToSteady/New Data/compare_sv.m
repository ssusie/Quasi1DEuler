clear all; close all; clc;

 extensions = {'5','999','9999'};
% extensions = {'999'};

figureIndex = 1;
for iExtension = 1:length(extensions)
	Gal = load(['G',extensions{iExtension}]);
	PetGal = load(['PG',extensions{iExtension}]);
	LSC = load(['Czim',extensions{iExtension}]);
    Gcons=load(['Gconstr', extensions{iExtension}]);
    PGcons=load(['PGconstr', extensions{iExtension}]);
    LSQcons=load(['constrZim', extensions{iExtension}]);
    
    FrobNormLSQ=norm(LSQcons.Czim,'fro');
    FrobNormPG=norm(PGcons.Corig,'fro');
    FrobNormG=norm(Gcons.Corig,'fro');
    normLSQconst_EachTime=zeros(1,size(LSQcons.Czim,2));
    normPGconst_EachTime=zeros(1,size(LSQcons.Czim,2));
    normGconst_EachTime=zeros(1,size(LSQcons.Czim,2));
    for i=1:size(LSQcons.Czim,2)
        normLSQconst_EachTime(i)=norm(LSQcons.Czim(:,i),2);
        normPGconst_EachTime(i)=norm(PGcons.Corig(:,i),2);
%         if i==23, keyboard, end
        if iExtension == 2 
%             keyboard
            normGconst_EachTime=[];
        else
           normGconst_EachTime(i)=norm(Gcons.Corig(:,i),2);   
        end
    end
    
    h100=figure(figureIndex);
    figureIndex = figureIndex + 1;
    plot(normLSQconst_EachTime,'b*', 'linewidth',2); hold on
    plot(normPGconst_EachTime,'ro', 'linewidth',2)
%     if extensions{iExtension} ~= 999 
     plot(normGconst_EachTime,'gv', 'linewidth',2)
%     end
    legend('Constrained', 'PG', 'G')
    xlabel('timestep')
    ylabel('norm of constraint')
    title(['norm of constraints on each timestep for energy ', extensions{iExtension}])
    
    h1=figure(figureIndex);
    figureIndex = figureIndex + 1;
	semilogy(abs(Gal.rom.sv(1:3:end,end)- PetGal.fom.sv(1:3:end,end)), 'g--', 'linewidth',2); hold on
	semilogy(abs(PetGal.rom.sv(1:3:end,end)- PetGal.fom.sv(1:3:end,end)), 'r--','linewidth',2); hold on
	semilogy(abs(LSC.rom.sv(1:3:end,end)- PetGal.fom.sv(1:3:end,end)), 'b--', 'linewidth',2);
	%plot(PetGal.fom.sv(1:3:end,end), 'k', 'linewidth',2);
	legend('G','PG','LSC')
	title(['rho at final time step for case ',extensions{iExtension}])


	h2=figure(figureIndex);
	figureIndex = figureIndex + 1;
	semilogy(abs(Gal.rom.sv(2:3:end,end)- PetGal.fom.sv(2:3:end,end)), 'g--', 'linewidth',2); hold on
	semilogy(abs(PetGal.rom.sv(2:3:end,end)- PetGal.fom.sv(2:3:end,end)), 'r--','linewidth',2); hold on
	semilogy(abs(LSC.rom.sv(2:3:end,end)- PetGal.fom.sv(2:3:end,end)), 'b--', 'linewidth',2);
	%plot(PetGal.fom.sv(2:3:end,end), 'k', 'linewidth',2);
	legend('G','PG','LSC')
	title(['rho*u at final time step',extensions{iExtension}])

	h3=figure(figureIndex);
	figureIndex = figureIndex + 1;
	semilogy(abs(PetGal.fom.sv(3:3:end,end) - Gal.rom.sv(3:3:end,end)), 'g--', 'linewidth',2); hold on
	semilogy(abs(PetGal.fom.sv(3:3:end,end) - PetGal.rom.sv(3:3:end,end)), 'r--','linewidth',2); hold on
	semilogy(abs(PetGal.fom.sv(3:3:end,end) - LSC.rom.sv(3:3:end,end)), 'b--', 'linewidth',2);
	%plot(PetGal.fom.sv(3:3:end,end), 'k', 'linewidth',2);
	legend('G','PG','LSC')
	title(['e at final time step',extensions{iExtension}])

	h10=figure(figureIndex);
	figureIndex = figureIndex + 1;
	A1=PetGal.fom.sv(:,:)-LSC.rom.sv(:,:);
	A2=PetGal.fom.sv-Gal.rom.sv;
	A3=PetGal.fom.sv-PetGal.rom.sv;
	for iTime = 1:size(A1,2)
			A1_time_errors(iTime) = norm(A1(:,iTime))/norm(PetGal.fom.sv(:,iTime));
			A2_time_errors(iTime) = norm(A2(:,iTime))/norm(PetGal.fom.sv(:,iTime));
			A3_time_errors(iTime) = norm(A3(:,iTime))/norm(PetGal.fom.sv(:,iTime));
			fom_time_norm(iTime) = norm(PetGal.fom.sv(:,iTime));
	end
	
	plot(A2_time_errors, 'g', 'linewidth',2); hold on
	plot(A3_time_errors, 'r', 'linewidth',2);
	plot(A1_time_errors, 'b', 'linewidth',2); 
	title(['norm of state-space error between FOM and ROM' ,extensions{iExtension}])
	xlabel('time iteration')
	ylabel('relative error')
	legend('G','PG','LSC')
	
	disp('LSC')
	norm(A1_time_errors,1)/norm(fom_time_norm,1)
	norm(A1_time_errors,2)/norm(fom_time_norm,2)
	norm(A1_time_errors,'inf')/norm(fom_time_norm,'inf')
	
	disp('Gal')
	norm(A2_time_errors,1)/norm(fom_time_norm,1)
	norm(A2_time_errors,2)/norm(fom_time_norm,2)
	norm(A2_time_errors,'inf')/norm(fom_time_norm,'inf')
	
	disp('LSPG')
	norm(A3_time_errors,1)/norm(fom_time_norm,1)
	norm(A3_time_errors,2)/norm(fom_time_norm,2)
	norm(A3_time_errors,'inf')/norm(fom_time_norm,'inf')
	
	clear A1 A2 A3

end
break

clear all; close all; clc;

Gal=load('Goriginal');
GalTrunc4=load('Goriginal_trunc4');
GalTrunc3=load('Goriginal_trunc3');

PetGal=load('PGoriginal50');
PetGalTrunc4=load('PGoriginal50_trunc4');
PetGalTrunc3=load('PGoriginal50_trunc3');

LSC=load('LSC50_5');
LSCTrunc4=load('LSC50_4');
LSC3=load('LSC50_3');

%figure(100)
%plot(gal.rom)



h1=figure(1);
plot(Gal.rom.sv(1:3:end,end), 'g--', 'linewidth',2); hold on
plot(PetGal.rom.sv(1:3:end,end), 'r--','linewidth',2); hold on
plot(LSC.rom.sv(1:3:end,end), 'b--', 'linewidth',2);
plot(PetGal.fom.sv(1:3:end,end), 'k', 'linewidth',2);
legend('G','PG','LSC','FOM')
title('rho with 5 basis vectors at final time step')


h2=figure(2);
plot(Gal.rom.sv(2:3:end,end), 'g--', 'linewidth',2); hold on
plot(PetGal.rom.sv(2:3:end,end), 'r--','linewidth',2); hold on
plot(LSC.rom.sv(2:3:end,end), 'b--', 'linewidth',2);
plot(PetGal.fom.sv(2:3:end,end), 'k', 'linewidth',2);
legend('G','PG','LSC','FOM')
title('rho*u with 5 basis vectors at final time step')


h3=figure(3);
plot(Gal.rom.sv(3:3:end,end), 'g--', 'linewidth',2); hold on
plot(PetGal.rom.sv(3:3:end,end), 'r--','linewidth',2); hold on
plot(LSC.rom.sv(3:3:end,end), 'b--', 'linewidth',2);
plot(PetGal.fom.sv(3:3:end,end), 'k', 'linewidth',2);
legend('G','PG','LSC','FOM')
title('e with 5 basis vectors at final time step')


h4=figure(4);
plot(GalTrunc4.rom.sv(1:3:end,end), 'g--', 'linewidth',2); hold on
plot(PetGalTrunc4.rom.sv(1:3:end,end), 'r--','linewidth',2); hold on
plot(LSCTrunc4.rom.sv(1:3:end,end), 'b--', 'linewidth',2);
plot(PetGalTrunc4.fom.sv(1:3:end,end), 'k', 'linewidth',2);
legend('G','PG','LSC','FOM')
title('rho with 4 basis vectors at final time step')



h5=figure(5);
plot(GalTrunc4.rom.sv(2:3:end,end), 'g--', 'linewidth',2); hold on
plot(PetGalTrunc4.rom.sv(2:3:end,end), 'r--','linewidth',2); hold on
plot(LSCTrunc4.rom.sv(2:3:end,end), 'b--', 'linewidth',2);
plot(PetGalTrunc4.fom.sv(2:3:end,end), 'k', 'linewidth',2);
legend('G','PG','LSC','FOM')
title('rho*u with 4 basis vectors at final time step')


h6=figure(6);
plot(GalTrunc4.rom.sv(3:3:end,end), 'g--', 'linewidth',2); hold on
plot(PetGalTrunc4.rom.sv(3:3:end,end), 'r--','linewidth',2); hold on
plot(LSCTrunc4.rom.sv(3:3:end,end), 'b--', 'linewidth',2);
plot(PetGalTrunc4.fom.sv(3:3:end,end), 'k', 'linewidth',2);
legend('G','PG','LSC','FOM')
title('e with 4 basis vectors at final time step')




h7=figure(7);
plot(GalTrunc3.rom.sv(1:3:end,end), 'g--', 'linewidth',2); hold on
plot(PetGalTrunc3.rom.sv(1:3:end,end), 'r--','linewidth',2); hold on
plot(LSC3.rom.sv(1:3:end,end), 'b--', 'linewidth',2);
plot(PetGalTrunc3.fom.sv(1:3:end,end), 'k', 'linewidth',2);
legend('G','PG','LSC','FOM')
title('rho with 3 basis vectors at final time step')



h8=figure(8);
plot(GalTrunc3.rom.sv(2:3:end,end), 'g--', 'linewidth',2); hold on
plot(PetGalTrunc3.rom.sv(2:3:end,end), 'r--','linewidth',2); hold on
plot(LSC3.rom.sv(2:3:end,end), 'b--', 'linewidth',2);
plot(PetGalTrunc3.fom.sv(2:3:end,end), 'k', 'linewidth',2);
legend('G','PG','LSC','FOM')
title('rho*u with 3 basis vectors at final time step')


h9=figure(9);
plot(GalTrunc3.rom.sv(3:3:end,end), 'g--', 'linewidth',2); hold on
plot(PetGalTrunc3.rom.sv(3:3:end,end), 'r--','linewidth',2); hold on
plot(LSC3.rom.sv(3:3:end,end), 'b--', 'linewidth',2);
plot(PetGalTrunc3.fom.sv(3:3:end,end), 'k', 'linewidth',2);
legend('G','PG','LSC','FOM')
title('e with 3 basis vectors at final time step')
break
%%%%%%%%%%
h10=figure(10);
A1=LSC3.fom.sv(1:3:end,:)-LSC.rom.sv(1:3:end,:);
A2=(LSC3.fom.sv(1:3:end,:)-Gal.rom.sv(1:3:end,:));
A3=(LSC3.fom.sv(1:3:end,:)-PetGal.rom.sv(1:3:end,:));

plot(diag(sqrt(A1'*A1)), 'r', 'linewidth',2); hold on
plot(diag(sqrt(A2'*A2)), 'b', 'linewidth',2);
plot(diag(sqrt(A3'*A3)), 'k', 'linewidth',2);
title('norm of rho_fom-rho_rom')
legend('LSC','G','PG')
clear A1 A2 A3
% break

h11=figure(11);
A1=LSC3.fom.sv(2:3:end,:)-LSC.rom.sv(2:3:end,:);
A2=LSC3.fom.sv(2:3:end,:)-Gal.rom.sv(2:3:end,:);
A3=LSC3.fom.sv(2:3:end,:)-PetGal.rom.sv(2:3:end,:);

plot(diag(sqrt(A1'*A1)), 'r', 'linewidth',2); hold on
plot(diag(sqrt(A2'*A2)), 'b', 'linewidth',2);
plot(diag(sqrt(A3'*A3)), 'k', 'linewidth',2);
title('norm of rho*u_fom-rho*u_rom')
legend('LSC','G','PG')
clear A1 A2 A3

% break
h12=figure(12);
A1=LSC3.fom.sv(3:3:end,:)-LSC.rom.sv(3:3:end,:);
A2=LSC3.fom.sv(3:3:end,:)-Gal.rom.sv(3:3:end,:);
A3=LSC3.fom.sv(3:3:end,:)-PetGal.rom.sv(3:3:end,:);

plot(diag(sqrt(A1'*A1)), 'r', 'linewidth',2); hold on
plot(diag(sqrt(A2'*A2)), 'b', 'linewidth',2);
plot(diag(sqrt(A3'*A3)), 'k', 'linewidth',2);
title('norm of e_fom-e_rom')
legend('LSC','G','PG')
clear A1 A2 A3

%%%%%%%%%%
h13=figure(13);
A1=LSC3.fom.sv(1:3:end,:)-LSCTrunc4.rom.sv(1:3:end,:);
A2=LSC3.fom.sv(1:3:end,:)-GalTrunc4.rom.sv(1:3:end,:);
A3=LSC3.fom.sv(1:3:end,:)-PetGalTrunc4.rom.sv(1:3:end,:);

plot(diag(sqrt(A1'*A1)), 'r', 'linewidth',2); hold on
plot(diag(sqrt(A2'*A2)), 'b', 'linewidth',2);
plot(diag(sqrt(A3'*A3)), 'k', 'linewidth',2);
title('norm of rho_fom-rho_rom with 4 basis vectors')
legend('LSC','G','PG')
clear A1 A2 A3


h14=figure(14);
A1=LSC3.fom.sv(2:3:end,:)-LSCTrunc4.rom.sv(2:3:end,:);
A2=LSC3.fom.sv(2:3:end,:)-GalTrunc4.rom.sv(2:3:end,:);
A3=LSC3.fom.sv(2:3:end,:)-PetGalTrunc4.rom.sv(2:3:end,:);

plot(diag(sqrt(A1'*A1)), 'r', 'linewidth',2); hold on
plot(diag(sqrt(A2'*A2)), 'b', 'linewidth',2);
plot(diag(sqrt(A3'*A3)), 'k', 'linewidth',2);
title('norm of rho*u_fom-rho*u_rom with 4 basis vectors')
legend('LSC','G','PG')
clear A1 A2 A3

h15=figure(15);
A1=LSC3.fom.sv(3:3:end,:)-LSCTrunc4.rom.sv(3:3:end,:);
A2=LSC3.fom.sv(3:3:end,:)-GalTrunc4.rom.sv(3:3:end,:);
A3=LSC3.fom.sv(3:3:end,:)-PetGalTrunc4.rom.sv(3:3:end,:);

plot(diag(sqrt(A1'*A1)), 'r', 'linewidth',2); hold on
plot(diag(sqrt(A2'*A2)), 'b', 'linewidth',2);
plot(diag(sqrt(A3'*A3)), 'k', 'linewidth',2);
title('norm of e_fom-e_rom with 4 basis vectors')
legend('LSC','G','PG')
clear A1 A2 A3

%%%%%%%%%%
h16=figure(16);
A1=LSC3.fom.sv(1:3:end,:)-LSC3.rom.sv(1:3:end,:);
A2=LSC3.fom.sv(1:3:end,:)-GalTrunc3.rom.sv(1:3:end,:);
A3=LSC3.fom.sv(1:3:end,:)-PetGalTrunc3.rom.sv(1:3:end,:);

plot(diag(sqrt(A1'*A1)), 'r', 'linewidth',2); hold on
plot(diag(sqrt(A2'*A2)), 'b', 'linewidth',2);
plot(diag(sqrt(A3'*A3)), 'k', 'linewidth',2);
title('norm of rho_fom-rho_rom with 3 basis vectors')
legend('LSC','G','PG')
clear A1 A2 A3

h17=figure(17);
A1=LSC3.fom.sv(2:3:end,:)-LSC3.rom.sv(2:3:end,:);
A2=LSC3.fom.sv(2:3:end,:)-GalTrunc3.rom.sv(2:3:end,:);
A3=LSC3.fom.sv(2:3:end,:)-PetGalTrunc3.rom.sv(2:3:end,:);

plot(diag(sqrt(A1'*A1)), 'r', 'linewidth',2); hold on
plot(diag(sqrt(A2'*A2)), 'b', 'linewidth',2);
plot(diag(sqrt(A3'*A3)), 'k', 'linewidth',2);
title('norm of rho*u_fom-rho*u_rom with 3 basis vectors')
legend('LSC','G','PG')
clear A1 A2 A3

h18=figure(18);
A1=LSC3.fom.sv(3:3:end,:)-LSC3.rom.sv(3:3:end,:);
A2=LSC3.fom.sv(3:3:end,:)-GalTrunc3.rom.sv(3:3:end,:);
A3=LSC3.fom.sv(3:3:end,:)-PetGalTrunc3.rom.sv(3:3:end,:);

plot(diag(sqrt(A1'*A1)), 'r', 'linewidth',2); hold on
plot(diag(sqrt(A2'*A2)), 'b', 'linewidth',2);
plot(diag(sqrt(A3'*A3)), 'k', 'linewidth',2);
title('norm of e_fom-e_rom with 3 basis vectors')
legend('LSC','G','PG')
clear A1 A2 A3






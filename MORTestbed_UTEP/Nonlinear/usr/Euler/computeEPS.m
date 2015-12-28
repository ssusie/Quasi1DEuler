U=fom.sv(:,125);
obj=prob;
U = reshape(U,3,obj.nVol);
%Extract primitive variables.  Take into account that the first
%three and last three variables are PRIMITIVE VARIABLES while
%all other are conservative variables.
[rho,u,P,c] = obj.conservativeToPrimitive(U(:,2:end-1));
rho = [U(1,1),rho,U(1,end)]; %Density
u   = [U(2,1),u,U(2,end)]; %Velocity
P   = [U(3,1),P,U(3,end)]; %Pressure
c   = [sqrt(obj.gamma*P(1)/rho(1)),c,sqrt(obj.gamma*P(end)/rho(end))]; %Speed of sound
e   = [P(1)/(obj.gamma-1)+rho(1)*u(1)^2/2,U(3,2:end-1),P(end)/(obj.gamma-1)+rho(end)*u(end)^2/2]; %Energy

[rhoH,uH,cH,drhoH,duH,dcH] = roeTerms(obj,rho,u,P,e);

lamH=uH'; lam1 = u(1:end-1)'; lam2 = u(2:end)';
eps1=obj.sigma0*max([zeros(99,1),lamH-lam1,lam2-lamH],[],2);
abs(lamH)<eps1


lamH=uH'+cH'; lam1 = u(1:end-1)'+c(1:end-1)'; lam2 = u(2:end)'+c(2:end)';
eps2=obj.sigma0*max([zeros(99,1),lamH-lam1,lam2-lamH],[],2);
abs(lamH)<eps2

lamH=uH'-cH'; lam1 = u(1:end-1)'-c(1:end-1)'; lam2 = u(2:end)'-c(2:end)';
eps3=obj.sigma0*max([zeros(99,1),lamH-lam1,lam2-lamH],[],2);
abs(lamH)<eps3

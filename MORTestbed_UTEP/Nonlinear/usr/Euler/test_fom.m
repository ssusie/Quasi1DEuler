clc

susie= load('Susie');
w_guess=obj.phi'*(susie.fom.sv(:,2)-susie.fom.sv(:,1));
save  true_guess w_guess
%%
obj.sv(:,itnump1)=obj.sv(:,obj.cTimeIter)+obj.phi*w_guess;
[Res,Jres]=obj.TimeScheme.TimeIntNLFunc(obj.prob,obj.sv(:,itnump1),obj.sv(:,itnump1-1),t);
Df=Jres*obj.phi;

[rho, u, P,c,e,dc_p,dc]=obj.prob.getVariables(obj.sv(:,itnump1));
[Q,dQ]=obj.prob.forceTerm(u,P);
[roeF, droeF]=obj.prob.roeFlux(rho,u,P,c,e,dc);

Rright= bsxfun(@times,roeF(:,2:end),obj.prob.S(3:end-1));
Rleft= -bsxfun(@times,roeF(:,1:end-1),obj.prob.S(2:end-2));

g(1:3,:)=[(obj.prob.SVol(2:end-1).*obj.prob.dx(2:end-1))*Phi1(2:end-1,:)*w_guess;(obj.prob.SVol(2:end-1).*obj.prob.dx(2:end-1))*Phi2(2:end-1,:)*w_guess; (obj.prob.SVol(2:end-1).*obj.prob.dx(2:end-1))*Phi3(2:end-1,:)*w_guess ]+...
     obj.time.dt*(Rright(:,end)+Rleft(:,1)-Q(:,2:end-1)*(obj.prob.SVol(2:end-1).*obj.prob.dx(2:end-1))')
     
     derivQ=zeros(3,5);
 for i=1:obj.prob.nVol
     derivQ=derivQ+squeeze(dQ(:,:,i))*[Phi1(i,:);Phi2(i,:);Phi3(i,:)]*(obj.prob.SVol(i).*obj.prob.dx(i));
 end

 dUdV=[1,0,0;u(1),rho(1),0;0.5*u(1)*u(1),rho(1)*u(1),1/(obj.prob.gamma-1)];
 droeF(:,1:3,1)=droeF(:,1:3,1)*dUdV;
 dUdV=[1,0,0;u(end),rho(end),0;0.5*u(end)*u(end),rho(end)*u(end),1/(obj.prob.gamma-1)];
 droeF(:,4:6,end)=droeF(:,4:6,end)*dUdV;
 
 J2L = -obj.prob.S(2)*droeF(:,1:3,1);
 J2R =  obj.prob.S(obj.prob.nVol)*droeF(:,4:6,obj.prob.nVol-1);
 
 Dg= [(obj.prob.SVol(2:end-1).*obj.prob.dx(2:end-1))*Phi1(2:end-1,:); (obj.prob.SVol(2:end-1).*obj.prob.dx(2:end-1))*Phi2(2:end-1,:); (obj.prob.SVol(2:end-1).*obj.prob.dx(2:end-1))*Phi3(2:end-1,:)]+...
     obj.time.dt*(J2L*[Phi1(1,:); Phi2(1,:); Phi3(1,:)] + ...
     J2R*[Phi1(end,:); Phi2(end,:); Phi3(end,:)]- derivQ);

 
%  J2L = [-obj.prob.S(2)*droeF(:,1:3,1),- obj.prob.S(k)*droeF(:,4:6,1)];
%  J2R = [obj.prob.S(obj.prob.nVol)*droeF(:,1:3,obj.prob.nVol-1),obj.prob.S(obj.prob.nVol)*droeF(:,4:6,obj.prob.nVol-1)];
% 
%  
%   Dg=[ [(obj.prob.SVol(2:end-1).*obj.prob.dx(2:end-1))*Phi1(2:end-1,:); (obj.prob.SVol(2:end-1).*obj.prob.dx(2:end-1))*Phi2(2:end-1,:); (obj.prob.SVol(2:end-1).*obj.prob.dx(2:end-1))*Phi3(2:end-1,:)]+...
%          obj.time.dt*([J2L(1,:)*(Phi1(1:6,:)*w_guess)*Phi1(1,:); J2L(2,:)*(Phi2(1:6,:)*w_guess)*Phi2(1,:); J2L(3,:)*(Phi3(1:6,:)*w_guess)*Phi3(1,:)] + ...
%                       [J2R(1,:)*Phi1(end-5:end,:)*w_guess*Phi1(end,:); J2R(2,:)*Phi2(end-5:end,:)*w_guess*Phi2(end,:); J2R(3,:)*Phi3(end-5:end,:)*w_guess*Phi3(end,:)]- derivQ)]
                              
                              
 H=Df'*Df; h=Df'*Res;
 P=H\Dg';
 x=H\h;
 S=-inv(Dg*P);
 Qzim=P*S;
 del_w=Qzim*g-(x+Qzim*Dg*x)
 
 
 
 
 
 
 
 
 
 
 
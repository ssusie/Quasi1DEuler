function [r1, dr1]=myroe(SV)
         
%                   SV=obj.partialUprev+obj.phiYhat*w_increment;
                  U = reshape(SV,3,obj.probGNAT.nVolSamp);
                  startInd = 1 + obj.probGNAT.ind1;
                  endInd   = obj.probGNAT.nVolSamp - obj.probGNAT.indN;
                  [rho,u,P,c] = obj.probGNAT.conservativeToPrimitive(U(:,startInd:endInd));
                  e=U(3,startInd:endInd);
                  if obj.probGNAT.ind1
                      rho = [U(1,1),rho]; %Density
                      u   = [U(2,1),u]; %Velocity
                      P   = [U(3,1),P]; %Pressure
                      c   = [sqrt(obj.probGNAT.gamma*P(1)/rho(1)),c]; %Speed of sound
                      e   = [P(1)/(obj.probGNAT.gamma-1)+rho(1)*u(1)^2/2,e]; %Energy
                  end
                  if obj.probGNAT.indN
                     rho = [rho,U(1,end)]; %Density
                      u   = [u,U(2,end)]; %Velocity
                      P   = [P,U(3,end)]; %Pressure
                      c   = [c,sqrt(obj.probGNAT.gamma*P(end)/rho(end))]; %Speed of sound
                      e   = [e,P(end)/(obj.probGNAT.gamma-1)+rho(end)*u(end)^2/2]; %Energy
                  end
     
                  dc_cons = [0.5*obj.probGNAT.gamma./(c.*rho).*(0.5*(obj.probGNAT.gamma-1)*u.*u - P./rho);...
                      -0.5*obj.probGNAT.gamma*(obj.probGNAT.gamma-1)*u./(rho.*c);...
                      0.5*obj.probGNAT.gamma*(obj.probGNAT.gamma-1)./(rho.*c)]';

                  %[rho, u, P,c,e,~,dc]=obj.prob.getVariables(obj.sv(:,itnump1));
                  [Q,dQ]=obj.probGNAT.forceTermGNAT(u,P);
                  %keyboard
                  [roeF, droeF]=obj.probGNAT.roeFluxGNAT(rho,u,P,c,e,dc_cons);

                 r1=roeF(1,1);
                 dr1=zeros(
                 
end
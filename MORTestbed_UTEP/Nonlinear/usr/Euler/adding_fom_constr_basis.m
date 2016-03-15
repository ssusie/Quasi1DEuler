
soln=fom.sv;

for jj=1:size(soln,2)
    
    U = reshape(soln(:,jj),3,obj.prob.nVol);
    [rho,u,P,c] = obj.prob.conservativeToPrimitive(U(:,2:end-1));
    rho = [U(1,1),rho,U(1,end)]; %Density
    u   = [U(2,1),u,U(2,end)]; %Velocity
    P   = [U(3,1),P,U(3,end)]; %Pressure
    c   = [sqrt(obj.prob.gamma*P(1)/rho(1)),c,sqrt(obj.prob.gamma*P(end)/rho(end))]; %Speed of sound
    e   = [P(1)/(obj.prob.gamma-1)+rho(1)*u(1)^2/2,U(3,2:end-1),P(end)/(obj.prob.gamma-1)+rho(end)*u(end)^2/2];
    dc_cons = [0.5*obj.prob.gamma./(c.*rho).*(0.5*(obj.prob.gamma-1)*u.*u - P./rho);...
        -0.5*obj.prob.gamma*(obj.prob.gamma-1)*u./(rho.*c);...
        0.5*obj.prob.gamma*(obj.prob.gamma-1)./(rho.*c)]';
    [roeF,droeF] = obj.prob.roeFlux(rho,u,P,c,e,dc_cons);
    [Q,dQ]=obj.prob.forceTerm(u,P);
    
    dUdV=[1,0,0;u(1),rho(1),0;0.5*u(1)*u(1),rho(1)*u(1),1/(obj.prob.gamma-1)];
    droeF(:,1:3,1)=droeF(:,1:3,1)*dUdV;
    dUdV=[1,0,0;u(end),rho(end),0;0.5*u(end)*u(end),rho(end)*u(end),1/(obj.prob.gamma-1)];
    droeF(:,4:6,end)=droeF(:,4:6,end)*dUdV;
    
    roeF1(:,jj) = roeF(1,:)';
    roeF2(:,jj) = roeF(2,:)';
    roeF3(:,jj) = roeF(3,:)';
    
    Fright(:,jj)= reshape(bsxfun(@times,roeF(:,2:end),obj.prob.S(3:end-1)), size(roeF,1)*size(roeF(:, 2:end),2),1);
    Fleft(:,jj) = reshape(-bsxfun(@times,roeF(:,1:end-1),obj.prob.S(2:end-2)),size(roeF,1)*size(roeF(:,2:end),2),1);
    forceQ(:,jj)= reshape(Q, size(Q,1)*size(Q,2), 1);
    J2L=zeros(3*(obj.prob.nVol-2), 3*obj.prob.nVol);
    J2R=J2L;
    for kl= 2:obj.prob.nVol-1
        i=kl-1;
        J2L(3*(i-1)+1:3*i,3*(kl-2)+1:3*kl) = [-obj.prob.S(kl)*droeF(:,1:3,kl-1), -obj.prob.S(kl)*droeF(:,4:6,kl-1)];
        J2R(3*(i-1)+1:3*i,3*(kl-1)+1:3*(kl+1)) = [obj.prob.S(kl+1)*droeF(:,1:3,kl), obj.prob.S(kl+1)*droeF(:,4:6,kl)];
    end
    
    
    if jj==1
        
        sq=Q(:,2:end-1)*(obj.prob.SVol(2:end-1).*obj.prob.dx(2:end-1))';
        sdq=zeros(3,obj.trunc);
        Phi = obj.phi(:,1:obj.trunc);
        for i=2:obj.prob.nVol-1
            sdq=sdq+squeeze(dQ(:,:,i))*Phi(3*i-2:3*i,:)*(obj.prob.SVol(i).*obj.prob.dx(i));
        end
        sq_true = sq;
        sdq_true = sdq;
        save fom_sq_true sq_true
        save fom_sdq_true sdq_true
        
        J2L_true=J2L;
        J2R_true=J2R;
        DforceQ=dQ;
        romdroeF= droeF;
        Q2_true=Q(2,:);
        
%         save romroeF roeF
        save fom_Fright Fright
        save fom_Fleft Fleft
        save fom_forceQ forceQ
        
        save fom_DforceQ DforceQ
        save fom_J2R_true J2R_true
        save fom_Q2_true Q2_true
        save fom_J2L_true J2L_true
        save fom_romdroeF romdroeF
    end
%     
%     NewFrBasis = [NewFrBasis, Fright(:,jj), J2R*obj.phi];
%     NewFlBasis = [NewFlBasis, Fleft(:,jj), J2L*obj.phi];
%     NewQBasis  = [NewQBasis, forceQ(:,jj), dQ*obj.phi];
    romroeF(:,jj)=roeF(:);
    if imag(Fright(:,jj))~=0
        disp(['righ flux is complex for snapshot number', num2str(jj)])
    elseif imag(Fleft(:,jj))~=0
        disp(['left flux is complex for snapshot number', num2str(jj)])
    elseif imag(forceQ)~=0
        disp(['force term is complex for snapshot number', num2str(jj)])
    end
end
\

% [ur2,sr2,~] = svd(NewFrBasis,0);
% [ul2,sl2,~] = svd(NewFlBasis,0);
[ur,sr,~]   = svd(Fright,0);
[ul, sl, ~] = svd(Fleft,0);
[uQ,sQ,~]   = svd(forceQ,0);
[uroe1,s1,~] = svd(roeF1,0);
[uroe2,s2,~] = svd(roeF2,0);
[uroe3,s3,~] = svd(roeF3,0);

% 
% obj.phiFr_fom = ur2;%(:,1:mFr);
% obj.phiFl_fom = ul2;%(:,1:mFr);
obj.phiFright_fom = ur;%(:,1:mF); %99.9%
obj.phiFleft_fom  = ul;%(:,1:mF);  %99.9%
% obj.phiRoeF1 = uroe1;%(:,1:mr);
% obj.phiRoeF2 = uroe2;%(:,1:mr);
% obj.phiRoeF3 = uroe3;%(:,1:mr);
obj.phiQ_fom = uQ;%(:,1:mQ);


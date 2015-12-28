function  [d2Rdu_r,d2Rdudz_r] = SecondDerivsPG(rom,R,J,dRdp,phi,h,t,w,order)
%Compute second order term in dcdw and dcdz for PG-ROM

if nargin < 9, order = 2; end;

r = size(phi,2);
np = length(rom.prob.p);

d2Rdu_r = zeros(r,r);
if nargout==2, d2Rdudz_r = zeros(r,np); else d2Rdudz_r=[]; end

switch order
    case 2
%INCREIBLY SLOW, INEFFICIENT WAY TO COMPUTE SECOND DERIVATIVES (used as a check)
%         d2Rdu_r2 = zeros(r,r);
%         d2Rdudz_r2 = zeros(r,np);
%         N=length(w);
%         d2Rdu=zeros(N,N);
%         d2Rdudz = zeros(N,np);
%         for i = 1:N
%             e=zeros(size(w)); e(i)=1;
%             [~,Jp,dRdp_p] = prob.ResSens(w + h*e,t);
%             [~,Jm,dRdp_m] = prob.ResSens(w - h*e,t);
%             for k = 1:N
%                 for j = 1:N
%                     d2Rdu(k,i,j) = (0.5/h)*(Jp(k,j)-Jm(k,j));
%                 end
%                 for j = 1:np
%                     d2Rdudz(k,i,j) = (0.5/h)*(dRdp_p(k,j) - dRdp_m(k,j));
%                 end
%             end
%         end
% 
%         for i = 1:length(w)
%             d2Rdu_r2 = d2Rdu_r2 + R(i)*phi'*squeeze(d2Rdu(i,:,:))*phi;
%             d2Rdudz_r2 = d2Rdudz_r2 + R(i)*phi'*squeeze(d2Rdudz(i,:,:));
%         end
        for j = 1:r
%             [~,Jp,dRdp_p] = rom.TimeScheme.TimeIntNLFuncSens(rom.prob,w + h*phi(:,j),w + h*phi(:,j),t);
%             [~,Jm,dRdp_m] = rom.TimeScheme.TimeIntNLFuncSens(rom.prob,w - h*phi(:,j),w - h*phi(:,j),t);
            [~,Jp,dRdp_p] = rom.prob.ResSens(w + h*phi(:,j),t);
            [~,Jm,dRdp_m] = rom.prob.ResSens(w - h*phi(:,j),t);
            
            Jr_p = Jp*phi; clear Jp;
            Jr_m = Jm*phi; clear Jm;
            
            d2Rdu_r(:,j) = (0.5/h)*R'*(Jr_p - Jr_m);
            if nargout==2, d2Rdudz_r(j,:) = (0.5/h)*R'*(dRdp_p - dRdp_m); end;
        end
    case 1
        Jr_m = J*phi; clear J;
        if nargout==2, dRdp_m = dRdp; clear dRdp; end
        for j = 1:r
            [~,Jp,dRdp_p] = rom.prob.ResSens(w + h*phi(:,j),t);
            
            Jr_p = Jp*phi; clear Jp;
            
            d2Rdu_r(:,j) = (1/h)*R'*(Jr_p - Jr_m);
            if nargout==2, d2Rdudz_r(j,:) = (1/h)*R'*(dRdp_p - dRdp_m); end;
        end
end
end
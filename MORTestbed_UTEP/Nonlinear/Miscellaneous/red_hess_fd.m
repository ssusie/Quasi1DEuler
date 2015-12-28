function  [H,H_gamma] = red_hess_fd(prob,nY,phi,y,t,gamma)
%H(j,k,i) = d^2 f_i/dy_j dy_k

h = 1e-4;
H=zeros(nY,nY,nY);
for k = 1:nY
    [~,Jr] = prob.ResJac(phi*y + h*phi(:,k),t);
    Jr_red = phi'*Jr*phi; clear Jr;
    
    [~,Jl] = prob.ResJac(phi*y - h*phi(:,k),t);
    Jl_red = phi'*Jl*phi;  clear Jl;
    
    for j = 1:nY
        for i = 1:nY
            H(j,k,i) = (0.5/h)*(Jr_red(i,j) - Jl_red(i,j));
        end
    end
end

H_gamma=[];
if nargin == 6
   H_gamma=zeros(nY,nY,nY);
   
   for i = 1:nY
       H_gamma(:,:,i) = gamma(:,:,i) + gamma(:,:,i)';
   end
end

end

%Slow!
% for i = 1:nY
%     for j = 1:nY
%         for k = 1:nY
%             [~,Jr] = fom.prob.ResJac(phi*y + h*phi(:,k),t);
%             Jr = Jr*phi(:,j);
%             
%             [~,Jl] = fom.prob.ResJac(phi*y - h*phi(:,k),t);
%             Jl = Jl*phi(:,j);
%             
%             H(j,k,i) = (0.5/h)*(phi(:,i)'*(Jr - Jl));
%         end
%     end
% end
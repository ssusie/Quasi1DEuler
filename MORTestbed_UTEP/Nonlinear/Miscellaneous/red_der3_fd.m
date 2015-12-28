function  [T,T_omega,err] = red_der3_fd(prob,nY,phi,y,t,omega)
%T(j,k,l,i) = d^3 f_i/dy_j dy_k dy_l

h = 5e-1;
T=zeros(nY,nY,nY,nY);

[~,J_red] = prob.ResJac(phi*y,t);
J_red=phi'*J_red*phi;

Jp_red = zeros(nY,nY,nY);
Jm_red = zeros(nY,nY,nY);
for t = 1:nY
    [~,Jp] = prob.ResJac(phi*y + h*phi(:,t),t);
    [~,Jm] = prob.ResJac(phi*y - h*phi(:,t),t);
    
    Jp_red(:,:,t) = phi'*Jp*phi;
    Jm_red(:,:,t) = phi'*Jm*phi;
end
clear Jp Jm;

for l = 1:nY
    for k = l:nY
        [~,Jp_c] = prob.ResJac(phi*y + h*(phi(:,k) + phi(:,l)),t);
        [~,Jm_c] = prob.ResJac(phi*y - h*(phi(:,k) + phi(:,l)),t);

        Jp_c_red = phi'*Jp_c*phi; clear Jp_c;
        Jm_c_red = phi'*Jm_c*phi; clear Jm_c;
            
        for j = k:nY
%             if (k == l)
%                 for i = 1:nY
%                     T(j,k,l,i) = (1/h/h)*(Jp_red(i,j,l) - 2*J_red(i,j) + Jm_red(i,j,l));
%                 end
%             elseif (j == k)
%                 for i = 1:nY
%                     T(j,k,l,i) = (1/h/h)*(Jp_red(i,l,j) - 2*J_red(i,l) + Jm_red(i,l,j));
%                 end
%             elseif (j == l)
%                 for i = 1:nY
%                     T(j,k,l,i) = (1/h/h)*(Jp_red(i,k,j) - 2*J_red(i,k) + Jm_red(i,k,j));
%                 end
%             else
                for i = 1:nY
                    T(j,k,l,i) = (0.5/h/h)*((Jp_c_red(i,j) - 2*J_red(i,j) + Jm_c_red(i,j)) ...
                                            - (Jp_red(i,j,k) - 2*J_red(i,j) + Jm_red(i,j,k)) ...
                                            - (Jp_red(i,j,l) - 2*J_red(i,j) + Jm_red(i,j,l)));
                end
%             end
        end
    end
end

for l = 1:nY
    for k = 1:nY
        for j = 1:nY
            ind = sort([j,k,l],'descend');
            if (j >= k) && (k >=l)
                continue;
            end
            T(j,k,l,:) = T(ind(1),ind(2),ind(3),:);
        end
    end
end

%Check for needed symmetry
err=0;
for i = 1:nY
    for j = 1:nY
        for k = 1:nY
            for l = 1:nY
                err = err + abs(T(j,k,l,i) - T(k,j,l,i)) + ...
                    abs(T(j,k,l,i) - T(j,l,k,i)) + ...
                    abs(T(j,k,l,i) - T(l,k,j,i));
            end
        end
    end
end

T_omega=[];
if nargin == 6
    T_omega=zeros(nY,nY,nY,nY);
    
    for i = 1:nY
        for j = 1:nY
            for k = 1:nY
                for l = 1:nY
                    T_omega(j,k,l,i) = omega(j,k,l,i) + omega(j,l,k,i) + ...
                                       omega(k,j,l,i) + omega(k,l,j,i) + ...
                                       omega(l,j,k,i) + omega(l,k,j,i);
                end
            end
        end
    end
end
end
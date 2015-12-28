function  [Zsamp] = lhs(zlow,zup,k,stream)

if size(zlow,1) ~= size(zup,1)
    error('Upper and lower bounds must be the same length');
end

if nargin == 3 || isempty(stream)
    stream=RandStream('mt19937ar','seed',71489);
end

n = size(zlow,1);

dz    = (1/k)*(zup - zlow);
Zpoints = zeros(n,k);
for i = 1:n
    for j = 1:k
        Zpoints(i,j) = zlow(i) + (j-1)*dz(i) + dz(i)*stream.rand();
    end
end

Zsamp = zeros(n,k);
Zsamp(1,:) = Zpoints(1,:);
for i = 2:n
    perm = findperm([1:k]',stream);
    Zsamp(i,:) = Zpoints(i,perm);
end

% if plotflag && n == 2
%     figure; axes; hold on;
%     for i = 1:k
%         plot([i,i]*dz(1),[zlow(2),zup(2)],'k--');
%     end
%     for i = 1:k
%         plot([zlow(1),zup(1)],[i,i]*dz(2),'k--');
%     end
%     plot(Zsamp(1,:),Zsamp(2,:),'bo')
% end
% 
% if plotflag && n == 3
%     figure; axes; hold on;
%     for i = 1:k
%         for j = 1:k
%             plot3([zlow(1),zup(1)],[i,i]*dz(2),[j,j]*dz(3),'k--');
%         end
%     end
%     for i = 1:k
%         for j = 1:k
%             plot3([i,i]*dz(1),[zlow(2),zup(2)],[j,j]*dz(3),'k--');
%         end
%     end
%     for i = 1:k
%         for j = 1:k
%             plot3([i,i]*dz(1),[j,j]*dz(2),[zlow(3),zup(3)],'k--');
%         end
%     end
%     plot3(Zsamp(1,:),Zsamp(2,:),Zsamp(3,:),'bo')
% end

end

function  [perm] = findperm(vec,stream)

n = size(vec,1);
perm = zeros(n,1);

for i = 1:n
    ind = stream.randi(length(vec));
    perm(i) = vec(ind);
    vec(ind,:) = [];
end

end
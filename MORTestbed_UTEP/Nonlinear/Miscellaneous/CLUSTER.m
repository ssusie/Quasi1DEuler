function  [Component] = CLUSTER(X)
%An automatic shape independent clustering technique

n = size(X,2);

[E,W] = determineRNG(X);
[SortedW,I] = sort(W);
SortedE = E(I,:);
[SortedE,SortedW] = deleteDuplicate(SortedE,SortedW);

m = size(SortedE,1);

if isempty(SortedW)
    Min = inf;
    Max = 0;
else
    Min = SortedW(1);
    Max = SortedW(m);
end

if Max < 2*Min
    fprintf('Cannot Cluster Data Further!\n');
    Component = X;
    return;
end

FirstOrder = diff(SortedW);
FirstOrder = sort(FirstOrder);
t = 0.5*(FirstOrder(1) + FirstOrder(m-1));
thresh = SortedW(diff(SortedW) >= t & SortedW(1:m-1,1) >= 2*Min);
thresh = max(thresh);

if isempty(thresh)
    fprintf('Cannot Cluster Data Further!\n');
    Component = X;
    return;
end

ind = (W <= thresh);
Ep = E(ind,:);
Wp = W(ind,:);

Component = connectedComponents(Ep,n);
k = length(Component);

if k == 1 
    fprintf('Cannot Cluster Data Further!\n');
    Component = X;
    return;
end

if k > sqrt(n)
    return;
else
    for i = 1:k
        Component{i} = CLUSTER(X(:,Component{i}));
    end
end

%Plot the clusters
% plotClusters(Component,1);

end

function  [E,W] = determineRNG(X)

n = size(X,2);
D = zeros(n);
%Compute the distance between each pair of points
for i = 1:n
    for j = 1:n
        if i > j
            continue;
        end
        D(i,j) = dist(X(:,i), X(:,j));
        D(j,i) = D(i,j);
    end
end

%if d(xi,xj) <= max(d(xi,xk),d(xj,xk)) for all k ~= i,j
%then there is an edge between xi and xj (i.e. they are
%relative neighbors)
E = []; W = [];
for i = 1:n
    for j = 1:n
        if i >= j
            continue;
        end
        flag = true;
        for k = 1:n
            if k == j || k == i
                continue;
            end
                        
            if D(i,j) > max(D(i,k),D(j,k))
                flag = false;
                break;
            end
        end
        if flag
            E = [E;[i,j]];
            W = [W;D(i,j)];
        end
    end
end
end

function  [d] = dist(x,y)

d = norm(x-y,1);

end

function [E,W] = deleteDuplicate(E,W)

m = size(E,1);
for i = m:-1:2
   if abs(W(i) - W(i-1)) < 1e-9
       if length(intersect(E(i,:),E(i-1,:))) == 2
           E(i,:) = [];
           W(i,:) = [];
       end
   end
end

end

function  [Comp] = connectedComponents(E,n)

G = spalloc(n,n,2*n);
for i = 1:size(E,1)
    G(E(i,1),E(i,2)) = 1;
    G(E(i,2),E(i,1)) = 1;
end

remainNodes = [1:n]';

Comp = [];
while ~isempty(remainNodes)
    C = findConnected(G,remainNodes(1),[]);
    Comp = [Comp; {C}];
    remainNodes = setdiff(remainNodes,C);
end
end

function  [visited] = findConnected(G,i,visited)

tmp = find(G(i,:)); tmp = tmp(:);
k = length(tmp);
visited = [visited;i];

% C = [];

for j = 1:k
    if isempty(intersect(visited,tmp(j)))
        unvisited = setdiff(tmp,visited);
        [tmpVis] = findConnected(G,unvisited(1),visited);
%         [tmpC,tmpVis] = findConnected(G,unvisited(1),visited);
%         C = [setdiff(tmp,visited);tmpC];
        visited = unique([visited;tmpVis]);
    end
end

end

function  [] = plotGraph(G,nodes)

h = figure;
ax = axes;

n = length(nodes);
X = rand(n,2);
plot(X(:,1),X(:,2),'ko','MarkerSize',16,'linewidth',2); hold on;

for i = 1:n
    text(X(i,1)-0.01,X(i,2),num2str(nodes(i)));
end

for i = 1:n
   for j = 1:n
      if i == j
          continue;
      end
      if G(nodes(i),nodes(j)) == 1
          x = [X(i,1);X(j,1)];
          y = [X(i,2);X(j,2)];
          v = [diff(x);diff(y)];
          plot(x + v(1)/norm(v)*0.015*[1;-1],y + v(2)/norm(v)*0.015*[1;-1],'r-','linewidth',2);
      end
   end
end


end
function  [Inc] = plotClusters(Component,inc)
Inc = inc;

if ~iscell(Component)
    for j = 1:size(Component,2)
        text(Component(1,j),Component(2,j),num2str(inc));
    end
    
    return;
end

for i = 1:length(Component)
    inc = inc + 1;
    Inc = plotClusters(Component{i},inc);
end
end
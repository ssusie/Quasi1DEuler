function  [out] = createModelField(model,n)

out = cell(n,1);
for i = 1:n
    out{i,1} = [model,num2str(i)];
end

end
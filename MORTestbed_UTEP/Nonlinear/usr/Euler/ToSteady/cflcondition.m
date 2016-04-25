function  [cfl] = cflcondition(n,obj)

if strcmpi(class(obj),'FOM')
    cfl = min(1e8,max(0.001*n^4,500/norm(obj.prob.ResJac(obj.sv(:,obj.cTimeIter)))^(0.7)));
elseif strcmpi(class(obj),'ROM')
    cfl = min(1e8,max(0.001*n^4,500/norm(obj.prob.ResJac(obj.sv(:,obj.cTimeIter)))^(0.7)));
else
    cfl = min(1e8,max(0.001*n^4,1/norm(obj.probGNAT.ResJacGNAT(obj.probGNAT.ic + obj.phiYhat*sum(obj.sv(:,1:obj.cTimeIter),2)))^(0.7)));
end


end
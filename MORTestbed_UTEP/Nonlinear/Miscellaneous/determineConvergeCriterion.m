function  [tol,tolIt] = determineConvergeCriterion(obj,convCrit)

if obj.newt.converge == 1
    tol = obj.newt.eps*convCrit;
    tolIt = [];
elseif obj.newt.converge == 2
    tol = obj.newt.eps(1);
    tolIt = obj.newt.eps(2);
else
    tol = max(obj.newt.eps(1)*convCrit,obj.newt.eps(2));
    tolIt = obj.newt.eps(3);
end

end

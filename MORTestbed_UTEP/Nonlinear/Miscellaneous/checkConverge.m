function  [conv] = checkConverge(obj,convCrit,du,tol,tolIt)


if convCrit < tol
    %If were are using convergence criteria 2 or 3, compute the
    %difference between two iterates and compute to absolute
    %tolerance
    if obj.newt.converge ~= 1
        itDiffNorm = norm(du,2);
        if itDiffNorm<tolIt
            conv=true;
            return;
        end
    else
        conv=true;
        return;
    end
end

if norm(du,2) == 0
    conv=true;
    return;
end

conv=false;

end
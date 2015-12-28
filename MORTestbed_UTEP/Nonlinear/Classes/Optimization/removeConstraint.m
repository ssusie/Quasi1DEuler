function  [Inew] = removeConstraint(Iold,lam)

Inew = Iold;
removeInd = lam < 0;
indInI = find(Iold);
Inew(indInI(removeInd)) = false;

end

function  [Inew] = addConstraint(Iold,gamma,alpha)

%
% [~,minGind] = min(gamma);
% newInd = ind(minGind(1));
%
% Inew = determineActiveSet(x,A,b);
% Inew = (Iold & Inew);
% Inew(newInd) = true;

Inew = Iold;

% ind = find(~Iold);
[alpha_max,minGind] = min(gamma);
if alpha == alpha_max
    Inew(minGind(1)) = true;
    %newInd = ind(minGind(1));
    %Inew(newInd) = true;
end
end

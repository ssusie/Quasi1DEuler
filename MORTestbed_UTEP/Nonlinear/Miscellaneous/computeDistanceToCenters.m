function [newCenter] = computeDistanceToCenters(obj,U,it)
%This function computes the distance from the local state vector U to the
%centers of the clusters using the precomputed quantities.
%--------------------------------------------------------------------------
%Inputs:
%-------
%obj   - ROM or locGNAT object
%U     - state vector that we wish to compute the distance between its
%        location in state space and the cluster centers
%
%Outputs:
%--------
%distC - nBases x 1 vector containing the distance between
%        UrLoc and the center of each cluster.
%--------------------------------------------------------------------------

if ((strcmpi(obj.basisUpdate,'border_exact') || strcmpi(obj.basisUpdate,'border') || (strcmpi(obj.basisUpdate,'border_fast_approx'))) && strcmpi(class(obj),'ROM'))
    z = zeros(obj.nBases,1);
    for i = 1:obj.nBases
        z(i) = norm(obj.sv(:,obj.cTimeIter) - obj.clusCenter(:,i),2);
    end
    [~,I] = sort(z);
    newCenter = I(1);
    return;
end

if it == 1
    z = reshape([obj.preCompDistQuant.d],obj.nBases,obj.nBases);
else
    z=zeros(obj.nBases);
    if strcmpi(obj.basisUpdate,'border_fast')
        for i1 = 1:obj.nBases
            for i2 = 1:obj.nBases
                if (i2 > i1)
                    z(i1,i2) = obj.preCompDistQuant(i1,i2).d;
                    for j = 1:obj.numswitch
                        z(i1,i2) = z(i1,i2) + obj.preCompDistQuant(i1,i2).e*obj.ROBcompon(j).alpha'*U{j};
                        for k = 1:obj.nBases
                            z(i1,i2) = z(i1,i2) + (obj.preCompDistQuant(i1,i2).f{k}*obj.ROBcompon(j).beta{k}' + ...
                                obj.preCompDistQuant(i1,i2).g{k}*obj.ROBcompon(j).N{k})*U{j};
                        end
                    end
                end
            end
        end
    elseif strcmpi(obj.basisUpdate,'border_fast_approx')
        for i1 = 1:obj.nBases
            for i2 = 1:obj.nBases
                if (i2 > i1)
                    z(i1,i2) = obj.preCompDistQuant(i1,i2).d;
                    for j = 1:obj.numswitch
                        z(i1,i2) = z(i1,i2) + obj.preCompDistQuant(i1,i2).h{j}*U{j};
                    end
                end
            end
        end
    elseif strcmpi(obj.basisUpdate,'none')
        for i1 = 1:obj.nBases
            for i2 = 1:obj.nBases
                if (i2 > i1)
                    z(i1,i2) = obj.preCompDistQuant(i1,i2).d;
                    for j = 1:obj.nBases
                        z(i1,i2) = z(i1,i2) + obj.preCompDistQuant(i1,i2).g{j}*U{j};
                    end
                end
            end
        end
    end
end

z = z-z';
for iLocBasis = 1:obj.nBases
    if sum(z(iLocBasis,:)<=0) == obj.nBases
        newCenter = iLocBasis;
        return;
    end
end
end
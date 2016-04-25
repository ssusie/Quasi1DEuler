

reconDD = (obj.reconQ2 * squeeze(dQ(1,:,2:10))')';

sumDerivQ2=zeros(1, length(w_increment));

for i=2:obj.problem.prob.nVol-1
    sumDerivQ2=sumDerivQ2+reconDD(:,i) * obj.problem.phi(3*i-2:3*i,:) * (obj.problem.prob.SVol(i).*obj.problem.prob.dx(i));
end

sdQ=sumDerivQ2;
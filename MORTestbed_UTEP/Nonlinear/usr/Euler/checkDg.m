points=1:floor(obj.prob.nVol/(obj.ncell)):obj.prob.nVol+1;

for k=1:10
D1(k)=(obj.prob.SVol(2:points(2)-1).*obj.prob.dx(2:points(2)-1))*Phi1(2:points(2)-1,k)+...
    obj.prob.S(points(2))*Phi2(points(2)-1,k)-obj.prob.S(2)*Phi2(2,k);
end



for k=1:10
D1(k)=(obj.prob.SVol(2:end-1).*obj.prob.dx(2:end-1))*Phi1(2:end-1,k)+...
    obj.prob.S(end-1)*Phi2(end-1,k)-obj.prob.S(2)*Phi2(2,k);
end


 for jj=1:298
if J2L_true(3*(jj-1)+1:3*jj,:)==JnewL(1:3,:)
jj
end
 end

 
 for jj=1:298
if J2L_true(3*(jj-1)+1:3*jj,:)==obj.JhatFluxL(4:6,1:6)
jj
end
end
function  ts=defineSurface(p,t,expr,plotflag)

tri=surftri(p,t);

p1=reshape(p(tri,1),size(tri));
p2=reshape(p(tri,2),size(tri));
p3=reshape(p(tri,3),size(tri));

ind=feval(expr,p1,p2,p3);

ts = tri(ind,:);

if plotflag
    simpplot(p,ts);
%     for i = 1:size(ts,1)
%         plot([p(ts(i,1),1),p(ts(i,2),1)],[p(ts(i,1),2),p(ts(i,2),2)],'k'); hold on;
%         plot([p(ts(i,1),1),p(ts(i,3),1)],[p(ts(i,1),2),p(ts(i,3),2)],'k');
%         plot([p(ts(i,2),1),p(ts(i,3),1)],[p(ts(i,2),2),p(ts(i,3),2)],'k');
%     end
%     axis equal;
end

end
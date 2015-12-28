function  [samp_feas] = MakeSampSatisfySplineConstrs(fname,dfom,samp_infeas,AlowBnd,AupBnd)

nsample = size(samp_infeas,2);
samp_feas=zeros(size(samp_infeas));

switch class(dfom.prob)
    case 'SteadyNozzle'
        xi = 0.5*(dfom.prob.mesh.node(1:end-1) + dfom.prob.mesh.node(2:end));
    case 'quasi1dEuler'
        xi = dfom.prob.mesh.gridpt;
end

if nargin < 4 || isempty(AlowBnd), AlowBnd = -1e9*ones(size(xi)); end;
if nargin < 5 || isempty(AupBnd), AupBnd  =  1e9*ones(size(xi)); end;

nlconstr = @(p) dfom.prob.NozzleConstr_Splines(p,xi,AlowBnd,AupBnd);
optionslhs = optimset('Algorithm','sqp',...
                    'Display','none',...
                    'GradObj','on',...
                    'GradConstr','on',...
                    'Hessian','bfgs',...
                    'TolFun',1e-6,...
                    'TolCon',1e-12,...
                    'TolX',1e-6);
                
optObj = PDEopt([fname,'.opt'],1,dfom);
for i = 1:nsample
    objoverw.objective='mindist';%objoverw.objective='feasibility'
    objoverw.point = samp_infeas(:,i);
    samp_feas(:,i) = optObj.fmincon_nand(samp_infeas(:,i),optionslhs,[],[],nlconstr,objoverw);
end
end

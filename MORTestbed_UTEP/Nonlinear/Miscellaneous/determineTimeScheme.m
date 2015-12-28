function  [obj] = determineTimeScheme(txt,model)


%if strcmpi(class(model.prob),'quasi1dEuler')
%    obj = quasiEuler1DImplicitPseudoTime(model);
%    return;
%end

if isempty(txt)
%    if strcmpi(class(model.prob),'HNL2dSS')
%        obj = SteadyTime(model);
%    else
        obj = BackwardEuler(model);
%    end
    return;
end

switch txt
    case 'quasiEuler1DImplicitPseudoTime'
        obj = quasiEuler1DImplicitPseudoTime(model);
    case 'Steady'
        %if strcmpi(class(model.prob),'HNL2dSS') || strcmpi(class(model.prob),'SteadyNozzle') || strcmpi(class(model.prob),'quasi1dEuler')
        obj = SteadyTime(model);
        %else
        %    fprintf('Cannot use Steady time stepping scheme with unsteady problem.  Using BackwardEule time stepping.');
        %end
    case 'BackwardEuler'
%         if ~(strcmpi(class(model.prob),'HNL2dSS') || strcmpi(class(model.prob),'SteadyNozzle') || strcmpi(class(model.prob),'quasiEuler1D'))
            obj = BackwardEuler(model);
%         else
%             fprintf('Cannot use Unsteady time stepping scheme with steady problem.  Using Steady.');
%             obj = SteadyTime(model);
%         end
    case 'ImplicitNewmark'
        obj = ImplicitNewmark(model);
    case 'specialPGsteady'
        obj = specialPGsteady(model);
    case 'GenAlpha'
        obj = GenAlpha(model);
end
end

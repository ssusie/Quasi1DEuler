function  [fomobj,probobj] = initialize(fname,cfgid)

for i = 1:length(cfgid)
    cfgobj(i) = CONFIG([pwd,SlashBS,fname,'.cfg'],cfgid(i));    
    switch cfgobj(i).problem
        case 1
            probobj(i) = OneDBurgers(cfgobj(i));
        case 2
            probobj(i) = NLTransLine(cfgobj(i));
        case 3
            probobj(i) = FHN(cfgobj(i));
        case 4
            probobj(i) = HNL2dSS(cfgobj(i));
        case 5
            probobj(i) = ConvDiffReact(cfgobj(i));
        case 6
            probobj(i) = MEMS(cfgobj(i));
        case 7
            probobj(i) = SteadyNozzle(cfgobj(i));
        case 8
            probobj(i) = quasi1dEuler(cfgobj(i));
        case 9
            probobj(i) = OneDFKPP(cfgobj(i));
        case 10
            probobj(i) = Lorenz(cfgobj(i));
        case 11
            probobj(i) = Contrived(cfgobj(i));
        case 12
            probobj(i) = Inverter(cfgobj(i));
        case 13
            probobj(i) = Rossler(cfgobj(i));
        case 14
            probobj(i) = OneDKdV(cfgobj(i));
        case 15
            probobj(i) = structuralFEM(cfgobj(i));
        case 16
            probobj(i) = lidINS(cfgobj(i));
        otherwise
            error(['In ',[pwd,SlashBS,fname,'.cfg'],': The problem number was not in the correct bounds.  See documentation for mapping between number and problem.']);
    end
    fomobj(i) = FOM([pwd,SlashBS,fname,'.cfg'],cfgid(i),probobj(i));
end

end
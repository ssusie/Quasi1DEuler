%This should be the filename of the CFG file
fname = 'example1000dof';
addGlobVarsAndPaths(fname);
%This defines the simulation with 1000 dofs.

% Create PROBLEM objects (define parameters for the problem in CFG file)
[fom,prob] = initialize(fname,1);

% Create wavelet ROM objects (WROM does not exist, it must be implemented)
% and solve at each parameter configuration
for i = 1:length(prob)
    wrom(i) = WROM([fname,'_wrom.rom'],'pg',i,prob(i));
    wrom(i).executeModel(); % This should generate the snapshots in wrom(i).sv
                            % using the wavelet ROM instead of the FOM
end

% Create snapshot-based ROM object
for i = 1:length(prob)
    rom(i) = ROM([fname,'.rom'],'pg',1,prob(i));
    rom(i).clusterSnapsOnly(wrom,wrom(1).sv(:,1)); % Pass your worm object, which has the snapshots
    rom(i).computePOD('fromClustSnaps');
    rom(i).executeModel
end

% Run FOM to check error
for i = 1:length(prob)
    fom(i).executeModel();
    romErr(i) = mean(ColumnwiseNorm(rom(i).sv-fom(i).sv,2)./ColumnwiseNorm(fom(i).sv,2));
end
disp(romErr)

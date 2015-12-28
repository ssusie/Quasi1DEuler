function  [SnapsLoc] = computeSnapshotFromDist(nsnap,param,dist,seedNum,lowBND,upBND)
%This function computes snapshot location from a given distribution
%--------------------------------------------------------------------------
%Inputs:
%-------
%nsnap    - scalar indicating the number of snapshots to collect
%param    - vector (size depends on the distribution) containing the
%           parameters of the distribution (max. time index for
%           uniform, mean and standard deviation for normal, )
%dist     - string indicating the distribution to use ('uniform', 'normal',
%           'beta', 'binomial', 'chi2', 'exponential', 'gamma',
%           'lognormal', 'negative binomial', 'f', 'poisson', 'weibull',
%           'rayleigh')
%seed     - integer specifying the seed number to use when generating
%           random numbers from the dist distribution
%lowBND   - integer specifying the smallest allowed time step
%upBND    - integer specifying the largest allowed time step
%
%Outputs:
%--------
%SnapsLoc - vector containing the time step numbers where snapshots are to
%           be taken
%--------------------------------------------------------------------------

%Create random stream
s = RandStream.create('mt19937ar','seed',seedNum);

%Initialize vector to store snapshot locations
SnapsLoc = zeros(1,nsnap);
%Initialize counter to track the number of snapshots we have taken
i = 1;

while i <= nsnap
    %Generate a random number from the appropriate distribution and make it
    %a whole number (since we are interested in time step numbers)
    switch dist
        case 'uniform'
            SnapsLoc = round(linspace(1,param,nsnap));
            return;
        case 'normal'
            temp = round(normrndSeed(param(1),param(2),s));
        case 'lognormal'
            temp = round(lognrndSeed(param(1),param(2),s));
        case 'exponential'
            temp = round(exprndSeed(param,s));
        %The distributions below are not supported yet.    
        case 'beta'
            temp = round(betarnd(param(1),param(2)));
        case 'binomial'
            temp = round(binornd(param(1),param(2)));
        case 'chi2'
            temp = round(chi2rnd(param));
        case 'gamma'
            temp = round(gamrnd(param(1),param(2)));
        case 'negative binomial'
            temp = round(nbinrnd(param(1),param(2)));
        case 'f'
            temp = round(frnd(param(1),param(2)));
        case 'poisson'
            temp = round(poissrnd(param));
        case 'weibull'
            temp = round(wblrnd(param));
        case 'rayleigh'
            temp = round(raylrnd(param(1),param(2)));
    end
    if (sum(SnapsLoc == temp) == 0) && (temp >= lowBND) && (temp <= upBND)
        %If we do not already have this location, store it and increment
        %the counter.  Otherwise, try again.
        SnapsLoc(i) = temp;
        i = i+1;
    end
end
%Sort the snapshots in ascending order
SnapsLoc = sort(SnapsLoc);
end
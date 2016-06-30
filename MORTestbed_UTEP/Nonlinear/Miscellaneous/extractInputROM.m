function  [out] = extractInputROM(txt,N,id,flag,def_val)
%This code extracts the model properties when using the multiple input
%notation (i.e. in the ROM file).
%--------------------------------------------------------------------------
%Inputs:
%-------
%txt      - string containing the text from the input file
%N        - number of id's for this model (i.e. number of instances of this
%           model)
%id       - index of model array to extract 
%flag     - string declaring which property to extract
%def_val  - the default value of the property
%
%Outputs:
%--------
%out - value of the field at the ind_th row
%--------------------------------------------------------------------------

%Determine the baseunit and basecols for each property
switch flag
    case {'initsnapref','basisUpdate','fileRoot','basisSelect','clustering','TimeScheme','lastTimeScheme','waveName'}
        baseunit = {'char'};
        basecols = nan;
    case {'svdUpdateBuffer','dt','nstep','maxsteps','maxIter','addElemTol','nBases','nSampleProd','printLevel'}
        baseunit = {'double'};
        basecols = 1;
    case 'linesrchProp'
        baseunit = {'double'};
        basecols = 8;
    case {'saveNL','steadyconverge'}
        baseunit = {'double','logical'};
        basecols = 1;
    case {'timeQuiet','newtQuiet','linesrch','i1','ie'}
        baseunit = {'logical'};
        basecols = 1;
    case 'T'
        baseunit = {'double'};
        basecols = 2;
    case 'addInd'
        baseunit = {'cell','char'};
        basecols = nan;
    case 'snapDist'
        baseunit = {'cell'};
        basecols = 3;
    case {'nY','nYMin','nYMax','nYrelEnergy','nI','nR','nJ','nSample','nFlux','nQ','nGreed','waveLvl','wavenY','sigma','mu'}
        baseunit = {'double'};
        basecols = nan;
    case 'snapInt'
        baseunit = {'double'};
        basecols = 2;
    case 'nsnap'
        baseunit = {'double','char'};
        basecols = [1,3];
    case 'eps'
        baseunit = {'double'};
        basecols = [1,2,3];
end

%Extract the line of code from txt that is going to be evaluted.
%Note that it should already be in MATLAB syntax.
tmp = determinePropText(txt,flag);
if isempty(tmp)
    out = def_val;
    return;
else
    %Evaluate the MATLAB expression in tmp to get a (possibly) vector or cell
    %array of entries.  We need to determine which row to use.
    prop  = evalin('caller',tmp);
end
% prop = extractInputRobust(txt,flag,def_val);

n = size(prop,1);

%If the property is empty, just return;
if n == 0
    out = prop;
    return;
end

%If the property is not of the appropriate size, throw error.
if n ~= N && n ~= 1
    error(['All properties in ROM file must have either 1 ROW or length(id) ROWs.  The ',flag,' property has incorrect number of ROWs. Cannot continue.']);
end

%If the input is a 'scalar', set id to 1 if it is not already
if n == 1
    id = 1;
end

%Determine if prop is in the appropriate baseunit or if it is wrapped in a cell
flagClass     = false;
flagClassCell = false;
for i = 1:length(baseunit)
    flagClass = flagClass || strcmpi(class(prop),baseunit{i});
    flagClassCell = flagClassCell || strcmpi(baseunit{i},'cell');
end

%Extract the appropriate quantity
if flagClass
    out = prop(id,:);
elseif iscell(prop) && ~flagClassCell;
    out = prop{id};
end

%Determine if any column sizes are violated (if basecols = nan, we don't
%know the appropriate number of columns apriori so it isn't checked)
flagColSize = false;
for i = 1:length(basecols)
    flagColSize = flagColSize || size(out,2) == basecols(i);
    if isnan(basecols(i))
        flagColSize = true;
    end
end

%Make sure out has the appropriate number of columns
if ~flagColSize
    error(['Not the appropriate number of columns in the ',flag,' field in the ROM file']);
end

end

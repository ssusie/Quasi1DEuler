function   [sValsR,sValsJ] = readResJacComputePODWriteFcn(obj,fileRoot,NLsnapcoll,svdOptions)
% function   [sValsR,sValsJ] =
% readResJacComputePODWriteFcn(obj,fileRoot,NLsnapcoll,lowRnkSVDapprox)

if nargin < 4, svdOptions=[]; end

sValsR=[]; sValsJ=[];

if strcmpi(class(obj(1)),'FOM') && NLsnapcoll == 0
    %Make sure all FOM objects are sorted into the same number of bases
    nbases = length(obj(1).res);
    for i = 1:length(obj)
        if length(obj(i).res) ~= nbases
            error('All FOM objects in the OBJ array must be separated into the same number of bases');
        end
    end
    
    %Don't include object if residuals were not saved
    ind2include = true(length(obj),1);
    for i = 1:length(obj)
        if obj(i).saveNL == 0
            ind2include(i) = false;
            fprintf('Residual Snapshots not saved because saveNL = 0 for object %i.  Continuing...\n',i);
            continue;
        end
    end
       
    %For each basis, extract residual snapshots and compute POD
    for j = 1:nbases
        [phiR,~,sValsR,~] = PODnonlin(obj(ind2include),j,'phiR',NLsnapcoll,svdOptions);
        writeNonlinBases(fileRoot,phiR,[],j); clear phiR;
    end
    
elseif strcmpi(class(obj(1)),'FOM') && NLsnapcoll == 0.5
    
    %Make sure all FOM objects are separated into the same number of bases
    nbases = length(obj(1).res);
    for i = 1:length(obj)
        if length(obj(i).res) ~= nbases
            error('All ROM objects in the OBJ array must have the same number of bases');
        end
    end
    
    %Don't include object if residuals were not saved
    ind2include = true(length(obj),1);
    for i = 1:length(obj)
        if obj(i).saveAllIt == 0
            ind2include(i) = false;
            fprintf('Residual Snapshots not saved because saveNL = 0 for object %i.  Continuing...\n',i);
            continue;
        end
    end
    
    %For each basis, extract residual snapshots and compute POD
    for j = 1:nbases
        [phiR,~,sValsR,~] = PODnonlin(obj(ind2include),j,'phiR',NLsnapcoll,svdOptions);
        writeNonlinBases(fileRoot,phiR,[],j); clear phiR;
    end
    
elseif strcmpi(class(obj(1)),'ROM')
    
    %Make sure all ROM objects have the same number of bases
    nbases = obj(1).nBases;
    for i = 1:length(obj)
        if obj(i).nBases ~= nbases
            error('All ROM objects in the OBJ array must have the same number of bases');
        end
    end
    
    for j = 1:nbases
        %For each basis, extract residual snapshots and compute POD
        if NLsnapcoll == 1 || NLsnapcoll == 1.5 || NLsnapcoll == 2
            [phiR,~,sValsR,~] = PODnonlin(obj,j,'phiR',NLsnapcoll,svdOptions);
            writeNonlinBases(fileRoot,phiR,[],j); clear phiR;
        end
        
        %For each basis, extract jacobian snapshots and compute POD
        if NLsnapcoll == 2
            [~,phiJ,~,sValsJ] = PODnonlin(obj,j,'phiJ',NLsnapcoll,svdOptions);
            writeNonlinBases(fileRoot,[],phiJ,j); clear phiJ;
        end
    end
end
end

function  [] = writeNonlinBases(fileRoot,phiR,phiJ,basisNum)
%Write phiR
if ~isempty(phiR)
    fid = fopen([fileRoot,'phiR',num2str(basisNum),'.bin'],'wb');
    fwrite(fid,size(phiR)','double');
    fwrite(fid,phiR,'double');
    fclose(fid);
end
%Write phiJ
if ~isempty(phiJ)
    fid = fopen([fileRoot,'phiJ',num2str(basisNum),'.bin'],'wb');
    fwrite(fid,size(phiJ)','double');
    fwrite(fid,phiJ,'double');
    fclose(fid);
end
end

function  [phiR,phiJ,sValsR,sValsJ] = PODnonlin(obj,basisNum,flag,NLSnapColl,svdOptions)
%This function extracts the nonlinear snapshots from basis number basisNum
%of each object in obj (based on NLSnapColl), computes the POD for each
%basis individually (uses svdOptions to determine how the POD is
%performed), and flag determines whether to do this for phiR or phiJ. 

%Initialize variables
sValsR=[]; sValsJ=[];
phiR=[]; phiJ=[];
RES=[]; JAC=[];

%%%%%PhiR%%%%%
for i = 1:length(obj)
    if strcmpi(flag,'phiR')
        %For each basis, extract residual snapshots and append to the
        %snapshot matrix
        [tmp,~] = obj(i).readResJacFiles(NLSnapColl,basisNum);
        RES=[RES,tmp]; clear tmp;
    end
end
if strcmpi(flag,'phiR')
    %Compute POD of residual snapshots
    [phiR,sValsR] = pod(RES,[],svdOptions); clear RES;
end

%%%%%PhiJ%%%%%
for i = 1:length(obj)
    if strcmpi(flag,'phiJ')
        %For each basis, extract jacobian snapshots and append to the
        %snapshot matrix
        [~,tmp] = obj(i).readResJacFiles(NLSnapColl,basisNum);
        JAC=[JAC,tmp]; clear tmp;
    end
end
if strcmpi(flag,'phiJ')
    %Compute POD of jacobian snapshots
    [phiJ,sValsJ] = pod(JAC,[],svdOptions); clear JAC;
end
end
function  [phiR,phiJ] = readNonlinBases(fileRoot,NLsnapcoll,basisNum)

phiR = [];
phiJ = [];

if isempty(basisNum)
    basisNum=1;
end

if NLsnapcoll == 0
    fid = fopen([fileRoot,'phiR1.bin'],'r');
    if fid == -1
        disp(['Could not open ',fileRoot,'phiR1.bin']);
        return;
    end
    sizePhiR = fread(fid,2,'double');
    phiR = reshape(fread(fid,prod(sizePhiR),'double'),sizePhiR(1),sizePhiR(2));
    fclose(fid);
end

if NLsnapcoll > 0
    fid = fopen([fileRoot,'phiR',num2str(basisNum),'.bin'],'r');
    if fid == -1
        disp(['Could not open ',fileRoot,'phiR',num2str(basisNum),'.bin']);
        return;
    end
    sizePhiR = fread(fid,2,'double');
    phiR = reshape(fread(fid,prod(sizePhiR),'double'),sizePhiR(1),sizePhiR(2));
    fclose(fid);
end

if NLsnapcoll == 2 || NLsnapcoll == 3
    fid = fopen([fileRoot,'phiJ',num2str(basisNum),'.bin'],'r');
    if fid == -1
        disp(['Could not open ',fileRoot,'phiJ',num2str(basisNum),'.bin']);
        return;
    end
    sizePhiJ = fread(fid,2,'double');
    phiJ = reshape(fread(fid,prod(sizePhiJ),'double'),sizePhiJ(1),sizePhiJ(2));
    fclose(fid);
end

end
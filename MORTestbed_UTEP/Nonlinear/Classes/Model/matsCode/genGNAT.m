function  [obj] = genGNAT(ROMfile,romobj,id)
%This is the constructor of the GNAT class.  It reads the appropriate
%entries from the input file and stores them in the class instance. 
%--------------------------------------------------------------------------
%Inputs:
%-------
%ROMfile - string indicating the filename of the rom file to use
%romobj  - ROM object that will be associated with this GNAT instance
%
%Outputs:
%--------
%obj     - instance of the GNAT object that was constructed
%--------------------------------------------------------------------------

if romobj.nBases == 1
    obj = GNAT(ROMfile,romobj,id);
else
    obj = locGNAT(ROMfile,romobj,id);
end

end
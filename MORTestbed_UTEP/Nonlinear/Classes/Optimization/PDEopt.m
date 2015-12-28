function  [optobj] = PDEopt(OPTfile,optid,modelobj)
%This function calls the constructor of the appropriate optimization
%function and returns the object instance.
%--------------------------------------------------------------------------
%Inputs:
%-------
%OPTfile  - string containing filename of OPT file
%optid    - id(s) from the OPT file to construct
%modelobj - model Object (FOM, ROM, GNAT, locGNAT)
%
%Outputs:
%--------
%optobj   - OPT object
%--------------------------------------------------------------

for i = optid(:)'
    %Extract parameters from CFG file
    GVARtxt = modelobj.GVARtxt;
    VARtxt     = [GVARtxt; readInFile('VAR',OPTfile,1)];
    %Extract variables from VARtxt and evaluate them in the current
    %workspace.
    if ~isempty(VARtxt)
        for j = 1:length(VARtxt)
            eval([VARtxt{j},';']);
        end
    end
    
    %Extract the OPT text from the cfg file
    OPTtext     = readInFile('OPT',OPTfile,optid);
    
    %Determine the OPT properties based on the text in the OPT file
    type   = extractModelPropMultInput(OPTtext,1,1,'type',[]);
    method = extractModelPropMultInput(OPTtext,1,1,'method',[]);
    
    switch type
        case 'bound constrained'
            %Box-Constrained Linesearch Minimization
            optobj = BoxConstMin(OPTfile,optid,modelobj);
        case 'nonlinearly constrained'
            switch method
                case 'SQP'
                    %SQP
                    %optobj = SQPold2(OPTfile,optid,modelobj);
                    optobj = SQP(OPTfile,optid,modelobj);
                case 'IntPt'
                case 'augLagPenFunc'
            end
    end
end
end
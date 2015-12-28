function  [err] = computeError(UROM,UFOM,nType,flag)
%This function computes the maximum relative error over all time steps
%(using the norm specified in nType)
%--------------------------------------------------------------------------
%Inputs:
%-------
%UROM  - The NdofsOfFOM x nstep (reconstructed) solution matrix from the ROM
%        analysis. Each column corresponds to the state vector at a
%        particular time step
%UFOM  - The NdofsOfFOM x nstep (reconstructed) solution matrix from the FOM
%        analysis. Each column corresponds to the state vector at a
%        particular time step
%nType - positive double (or inf) specifying the type of norm to use in
%        calculating the maximum relative error
%type  - string indicating whether to compute the relative or absolute
%        error ('rel' = relative, 'abs' = absolute).
%
%Outputs:
%--------
%err   - non-negative real number indicating the maximum relative 2-norm
%        error in the ROM approximation of the FOM over all time steps
%--------------------------------------------------------------------------

%Compute absolute error vectors for each time step
UErr = UROM - UFOM; 

%Compute the ratio (i.e. the relative error) of the nType-norm of the
%absolute error vector at each time step over the nType-norm of the FOM
%state vector at each time step.  Then find the maximum error over all time
%steps.
if lower(flag(1)) == 'r'
    err = ColumnwiseNorm(UErr,nType)./ColumnwiseNorm(UFOM,nType);
elseif lower(flag(1)) == 'a'
    err = ColumnwiseNorm(UErr,nType);
end

end
function  [phi,newk] = PODCompression(modelobj,k,type)
%This function calculates the Proper Orthogonal Decomposition for a problem
%in the time domain.
%--------------------------------------------------------------------------
%Inputs:
%-------
%modelobj - a n x Nsnap matrix (n is the full order model size)
%k        - scalar specifying the size of the reduced basis
%type     - string indicating the collection to compress ('sv' to compress
%           the state vector snapshot matrix, 'res' to compress the
%           residual snapshot matrix,'jac' to compress the jacobian
%           snapshot matrix
%
%Outputs:
%--------
%phi  - the POD compressed basis
%newk - the new model order after truncation based on rank of PHI
%--------------------------------------------------------------------------

if strcmpi(type,'sv')
    [phi,newk] = pod(modelobj.sv,k);
elseif strcmpi(type,'res')
    [phi,newk] = pod(modelobj.res,k);
elseif strcmpi(type,'jac')
    [phi,newk] = pod(modelobj.jac,k);    
end

end
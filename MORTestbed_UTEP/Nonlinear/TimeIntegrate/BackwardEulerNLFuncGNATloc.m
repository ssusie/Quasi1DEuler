function  [Rhat,JVhat] = BackwardEulerNLFuncGNATloc(modelobj)
%This function sets up parameters to be used in Backward Euler calculation
%dX/dt = R(X) where R is the residual of the nonlinear function
%--------------------------------------------------------------------------
%Inputs:
%-------
%modelobj     - Model object
%t            - 1 x (nstep+1) vector containing the times at which the
%               solution is to be computed (first entry is the time of the
%               initial condition)
%dt           - scalar indicating the time step size
%Outputs:
%--------
%R            - the residual of the nonlinear function, adjusted for the
%               Backward Euler scheme (wrapped in the nonlinear solver)
%J            - the jacobian of the nonlinear function, adjusted for the
%               Backward Euler scheme (wrapped in the nonlinear solver)
%--------------------------------------------------------------------------

%Determine the Residual and Jacobian (function handles) of the continuous
%system of equations
itnum = modelobj.cTimeIter;
t = modelobj.time.T(1) + modelobj.time.dt*itnum;
[Rhat,Jhat] = modelobj.probGNAT(modelobj.cLocBasis).ResJacGNAT(modelobj.partialULoc(modelobj.uniqueJROWind{modelobj.cLocBasis}),t);
% [Rhat,Jhat] =
% probobj.ResJacGNAT(modelobj.partialU,t,modelobj.sampleInd,modelobj.irstart,modelobj.jrow);
%Add contribution to the residual for the Backward Euler scheme)
Rhat = Rhat - (modelobj.partialULoc(modelobj.jdiagHat{modelobj.cLocBasis} ,1)-modelobj.partialUprevLoc(modelobj.jdiagHat{modelobj.cLocBasis} ,1))./modelobj.time.dt;
% Rhat = Rhat - (modelobj.partialU(logical(modelobj.jdiag),1)-modelobj.partialUprev(logical(modelobj.jdiag),1))./modelobj.time.dt;
%Set up the appropriate jacobian function for the
%backward euler integration scheme used (see notes for
%details)
JVhat = modelobj.computeJPHIhat(Jhat - modelobj.jdiag{modelobj.cLocBasis}./modelobj.time.dt);
end
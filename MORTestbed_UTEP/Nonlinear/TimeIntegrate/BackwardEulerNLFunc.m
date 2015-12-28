function  [R,J] = BackwardEulerNLFunc(modelobj)
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
[R,J] = modelobj.prob.ResJac(modelobj.sv(:,itnum+1),t);
%Add contribution to the residual for the Backward Euler scheme)
R = R - (modelobj.sv(:,itnum+1)-modelobj.sv(:,itnum))./modelobj.time.dt;
%Set up the appropriate jacobian function for the
%backward euler integration scheme used (see notes for
%details)
J = J - speye(modelobj.prob.config.ndof)./modelobj.time.dt;
end

function [v,dvdr] = phi4 ( r, r0 )

%*****************************************************************************80
%
%% PHI4 evaluates the gaussian radial basis function.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    29 June 2012
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
%    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
%    Third Edition,
%    Cambridge University Press, 2007,
%    ISBN13: 978-0-521-88068-8,
%    LC: QA297.N866.
%
%  Parameters:
%
%    Input, real R(), the radial separation.
%    0 < R.
%
%    Input, real R0, a scale factor.
%
%    Output, real V(), the value of the radial basis function.
%
  v = exp ( - 0.5 * r.^2 / r0^2 );
  dvdr = -1/r0^2*diag(r.*v);
  return
end

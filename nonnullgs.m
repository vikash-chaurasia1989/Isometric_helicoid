function [XN, info] = nonnullgs (X, x, b, M)

% Perform one step of the interval gauss seidel method to solve the problem
% M(z - x) = b inside a box X. This function works on lines where the
% diagonal element M(i,i) does not contain zero.
%
% It is a internal Function. To be used with intmincon, intminunc and
% intsolve functions
%
% This is an INTLAB file. It requires to have INTLAB installed under
% MATLAB to function properly.
%
% This algorithm is based on section 5.7 of
% Global Optimization Using Interval Analysis. Hansen E., Walster G. W.
% 2004. Marcel Dekker inc.
%
% Arguments with * are needed.
%
% Argument    i/o   description
%   X*         i    Initial Box.
%
%   x*         i    A numerical (non interval) vector inside the box X. The
%                   expansion point of the problem above
%
%   b*         i    An interval vector. The left side in equation above.
%
%   M*         i    An interval preconditioned matrix.
%
%   XN         o    The final box.
%
%   info       o    Information parameter.
%                     -1 The Gauss Seidel method produces a NaN entry on
%                     XN and the initial box can be deleted. See the
%                     reference above.
%                     0  The Gauss Seidel method produces a new box XN.
%
%    WARRANTY
%
%    Because the program is licensed free of charge, there is
%    no warranty for the program, to the extent permitted by applicable
%    law. Except when otherwise stated in writing the copyright holder
%    and/or other parties provide the program "as is" without warranty
%    of any kind, either expressed or implied, including, but not
%    limited to, the implied warranties of merchantability and fitness
%    for a particular purpose. The entire risk as to the quality and
%    performance of the program is with you. Should the program prove
%    defective, you assume the cost of all necessary servicing, repair
%    or correction.
%
%    History
%
%    02-13-2009   first version

  XN   = X;
  n    = length(XN);
  info = 0;  

  if ~in(0,M(1,1))
   N     = x(1) + (b(1) - M(1,2:n) * (XN(2:n)-x(2:n))') / M(1,1);
   XN(1) = intersect(N,X(1)); 
   
   if isNaN(XN(1)) || isempty(XN(1))
    info = -1;   
    return;
   end
  end
  
  
  for i = 2:n
    if ~in(0, M(i,i))
     N     = x(i) + (b(i) - M(i,1:i-1) * (XN(1:i-1)-x(1:i-1))' - ...
            M(i,i+1:n) * (XN(i+1:n) - x(i+1:n))') / M(i,i);
     XN(i) = intersect(N,X(i));
     
     if isNaN(XN(i)) || isempty(XN(i))
      info = -1;   
      return;
     end
    end
  end
end
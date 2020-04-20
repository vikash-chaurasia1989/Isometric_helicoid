function [XN, gap, info] = nullgs (X, x, b, M)

% Perform one step of the interval gauss seidel method to solve the problem
% M(z - x) = b inside a box X. This function works on lines where the
% diagonal element M(i,i) contains zero. 
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
%   gap        o    The widest gap obtained during the process. See the
%                   reference
%
%   info       o    Information parameter.
%                     -1 The Gauss Seidel method produces a NaN entry on
%                     XN and the initial box can be deleted. See the
%                     reference.
%                     0  The Gauss Seidel method produces a new box XN
%                     and no gap has been founded.
%                     > 0  The Gauss Seidel method produces a new box XN. 
%                     at least one gap has been founded during the process
%                     and the index of the widest is given by info.
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

  XN     = X;
  n      = length(XN);  
  gap    = [];
  maxgap = 0;
  info   = 0;
  
  if in(0,M(1,1))
   R   = b(1) - M(1,2:n) * (XN(2:n)-x(2:n))';
   div = intdiv(R, M(1,1));
   y   = intersect(x(1) + div, X(1));
   
   if all(isNaN(y))
    info = -1;
    return;
   end
   
   if length(y) == 1
    XN(1) = y;
   else
    if ~isNaN(intersect(y(1), y(2)))
     XN(1) = hull(y(1),y(2));
    else
     if y(1).inf < y(2).inf
      gw = y(2).inf - y(1).sup;
     else
      gw = y(1).inf - y(2).sup;
     end
     
     if gw > maxgap
      info   = 1;
      gap    = [y(1) y(2)];
     end
    end
   end
  end 
 
  for i = 2:n
    if in(0,M(i,i))
     R   = b(i) - M(i,1:i-1) * (XN(1:i-1)-x(1:i-1))' - ...
                M(i,i+1:n) * (XN(i+1:n)-x(i+1:n))';
     div = intdiv(R, M(i,i));
     y   = intersect(x(i) + div, X(i));
   
     if all(isNaN(y))
      info = -1;
      return;
     end
   
     if length(y) == 1
      XN(i) = y;
     else
      if ~isNaN(intersect(y(1), y(2)))
       XN(i) = hull(y(1),y(2));
      else
       if y(1).inf < y(2).inf
        gw = y(2).inf - y(1).sup;
       else
        gw = y(1).inf - y(2).sup;
       end
       
       if gw > maxgap
        info   = i;
        gap    = [y(1) y(2)]; 
       end
      end
     end
    end
  end
end
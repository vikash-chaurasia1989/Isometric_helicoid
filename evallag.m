function [y, dy] = evallag (f, c, x, lambda)

% Evaluate the lagrangian system and it's first derivative in a pair 
% (x, lambda) for the objective function f and equality constraints c.
%
% It is a internal Function. To be used with intmincon function.
%
% This is an INTLAB file. It requires to have INTLAB installed under
% MATLAB to function properly.
%
% Arguments with * are needed.
%
% Argument    i/o   description
%   f*         i    Objective function.
%
%   c*         i    Equality constraints.
%
%   x*         i    A numerical or interval vector where we want to
%                   evaluate de function
%
%   lambda*    i    Numerical or interval vector. The Lagrange multipliers
%                   which we use on evaluation.   
%
%   y          o    Numerical or interval vector. Lagrangian system
%                   applied in (x,lambda).
%
%   dy         o    Numerical or interval matrix. The Jacobian matrix of
%                   the Lagrangian system applied in (x,lambda).
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

  if nargout == 1
   dc  = feval(c, gradientinit(x));
   df  = feval(f, gradientinit(x)); 
   lag = df(:).dx' - sum((diag(lambda) * dc(:).dx),1)';
   y   = [lag; dc(:).x];
  end

  if nargout == 2
   dc  = feval(c, hessianinit(x));
   df  = feval(f, hessianinit(x));
   lag = df(:).dx' - sum((diag(lambda) * dc(:).dx),1)';    
   y   = [lag; dc(:).x];
   hl  =  df(:).hx;
   
   for i = 1: length(lambda)
     hl = hl - lambda(i) * dc(i).hx;
   end
   
   dy = [ hl, -dc(:).dx'; dc(:).dx, zeros(length(lambda))];
  end
end
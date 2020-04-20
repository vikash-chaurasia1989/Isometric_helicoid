function [x, ival, exitflag, output, Kraw, res] = intksolve(f, x0, options)

% Apply the Krawczyk operator to test existence of zeros of a nonlinear
% system of equations inside a box x0. If the second Moore's condition
% is satisfied the algorithm returns the unique zero of f inside x0.
%
% Example of Usage:
% >> f  = @(x) [x(1)^2+ x(2)^2 - 1, x(1)^2 - x(2)];
% >> x0 = [infsup(.7, .9), infsup(.5, .7)];
% >> x  = intksolve(f, x0)
% intval x = 
%  [0.7861, 0.7862] [0.6180, 0.6181] 
%
% This is an INTLAB file. It requires to have INTLAB installed under
% MATLAB to function properly.
%
% This method uses automatic differentiation to evaluate the jacobian 
% matrix in a box X. The user just need to supply the objective function f
% but he/she must to be sure that this function is differentiable in an open
% set wich contains the initial box x0.
%
% The algorithm presented here is based on the classical Moore's Paper
% A Test for Existence of Solutions to Nonlinear Systems. Moore R. E., 1977
% Siam Journal of Numerical Analysis.Vol 14, No 4. 611 - 615.
%
% Users can set the following parameters: maximum number of function eva
% luations, maximum number of outer iterates and the level of display.
% See the intoptimset function.
%
% Arguments with * are needed.
%
% Argument    i/o   description
%   f*         i    Objective Function.
%
%   X0*        i    Initial Interval(s). A row vector containing the
%                   initial box
%
%   options    i    Vector of options. See intoptimset and intoptimget.
%
%   x          o    Vector with interval Solution. If it is unique.
%
%   ival       o    Vector with function evaluated at solution. If it is
%                   unique.
%
%   exitflag   o    Integer identifying the reason the algorithm terminated
%                     -1 Algorithm finished with success. There are at 
%                        least one solution of f inside x0.  
%                     0  Algorithm finished with success. There is a unique 
%                        solution founded with Krawczyk operator. 
%                     1  Algorithm finished with maximum iterates.
%                     2  Algorithm finished with maximum function
%                     evaluations.
%                     3  Moore's Condition does not hold. Verify the
%                     reference above.
%
%   output     o    Information vector
%                     output(1) = number of outer iterations.
%                     output(2) = number of function evaluations.
%                     output(3) = number of jacobian evaluations.
%                     output(4) = number of matrix invertions.
%                     output(5) = percentage of the domain reached by the
%                     method.
%
%   Kraw       o    Krawczyk operator applied to initial box x0. See the
%                   reference above.
%
%   res        o    The value of ||eye(n) - YJ(X_{0})||. See the reference
%                   above
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

  if nargin < 3
   opt = intoptimset();
   
   if nargin < 2
    error('Missing Arguments. You must specify an objective function and an initial box.');
   end
  else
   opt = options;
  end
  
  maxIte   = opt(5);
  maxFunEv = opt(6);
  display  = opt(8);
  
  if opt(9) == 1
   maxIte   = inf;
   maxFunEv = inf;
  end
  
  n        = length(x0);
  exitflag = -1;
  x        = x0;
  ival     = intval(zeros(1,n));
  total    = prod(diam(x0));
  ite      = 0;
  evf      = 0;
  evj      = 0;
  evi      = 0;
  pct      = 0;
  y        = mid(x);
  ef       = feval(f, intval(y));
  eF       = feval(f, gradientinit(x));
  [r,c]    = size(ef);  
  
  if c > r
   ef = ef';
  end
  
  evf      = evf + 2;
  evj      = evj + 1;
   
  if evf > maxFunEv
   exitflag = 2;
   output   = [ite evf evj evi pct];
   return;
  end
  
  Y    = inv(mid(eF(:).dx));
  evi  = evi + 1;
  K    = y' - Y * ef + (eye(n) - Y * eF(:).dx)* (x - y)';
  Kraw = K;
  
  if ~all(in(K',x))
   fprintf('Moore`s Conditions does not hold.\n');
   exitflag = 3;
   output   = [ite evf evj evi pct];
   return;
  end
  
  res = norm(eye(n) - Y * eF(:).dx,inf);
  res = res.sup;
  
  if res < 1
   exitflag = 0;   
  else
   output   = [ite evf evj evi pct];
   return;
  end
  
  vx  = prod(diam(x));
  xn  = intersect(x, K');
  ite = ite + 1;
        
  if ite > maxIte
   exitflag = 1;
   output   = [ite evf evj evi pct];
   return;
  end
 
  while xn ~= x      
    vnx   = prod(diam(xn));
    pct   = pct + (vx - vnx) / total; 
    
    if display > 0
     if mod(ite, 20 * display) == 0
      clc;
      fprintf('\n ite \t evf \t evj \t evi \t reached \n');
     end
     
     if mod(ite, display) == 0
      fprintf('\n %d \t %d \t %d \t %d %.3f\n', ite, evf, evj, evi, pct);
     end
    end
    
    x     = xn;
    vx    = vnx;
    y     = mid(x);
    ef    = feval(f, intval(y));
    eF    = feval(f, gradientinit(x));
    [r,c] = size(ef);  
  
    if c > r
     ef = ef';
    end
  
    evf      = evf + 2;
    evj      = evj + 1;
   
    if evf > maxFunEv
     exitflag = 2;     
     break;
    end
  
    Y   = inv(mid(eF(:).dx));
    evi = evi + 1;
    K   = y' - Y * ef + (eye(n) - Y * eF(:).dx)* (x - y)';    
    xn  = intersect(x, K');
        
    ite = ite + 1;
        
    if ite > maxIte
     exitflag = 1;
     break;
    end
  end
  
  pct = pct + vx / total;
  
  if nargout >= 4
   output = [ite evf evj evi pct];
  end
end
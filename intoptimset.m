function options = intoptimset()

% Create an optimization options vector. To change default settings on
% Intesolve methods you must run this function specifying an output
% vector and change the needed settings. Rum intoptimget function to see 
% the parameters setted in an options vector
%
% Example of Usage: 
% >> opt = intoptimset();
% >> opt(1) = 10^-10; % Set the function tolerance to 10^-10.
% >> opt(8) = 1;      % Change the display level to 1.
% >> x = intminunc(func, x0, opt);
%
% This is an INTLAB file. It requires to have INTLAB installed under
% MATLAB to works properly.
%
% Arguments with * are needed.
%
% Argument    i/o   description
%   options*   o    Options vector
%                    options(1) Function tolerance. epsF
%
%                    options(2) Box tolerance.epsX
%
%                    options(3) Equality constraints tolerance. epsEC 
%                    (Constrained optimization problems only).
%
%                    options(4) Initial guess for global minimum. (Optimiza
%                    tion problems only).
%
%                    options(5) Maximum number of iterations allowed.
%
%                    options(6) Maximum number of function evaluations
%                    allowed.
%
%                    options(7) Maximum number of equality constraint
%                    evaluations allowed. Constrained optimization problems
%                    only.
%
%                    options(8) Level of display 
%                            0 = no information.
%                            n = display informations at each n iterations.
%
%                    options(9) Iterate until convergence. If true ignores
%                    any values in options 5, 6 and 7. 
%
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

  if nargout < 1
   error('Missing arquments. you must specify an output vector (For example opt = intoptimset())');
  end
  
  options(1)  = 10^-8;    % Function Tolerance (epsF).
  
  options(2)  = 10^-8;    % Box Tolerance (epsX).
  
  options(3)  = 10^-8;    % Equality constraints tolerance (epsEC). 
                          % Constrained Problems only.
                          
  options(4)  = inf;      % Initial guess for the global minimum. 
                          % Optimization problems only.
  
  options(5)  = 1000;     % Maximum number of iterations allowed.
  
  options(6)  = 10000;    % Maximum number of function evaluations allowed.
  
  options(7)  = 10000;    % Maximum number of equality constraint
                          % evaluations allowed. Constrained problems only.  
  
  options(8)  = 10;       % Level of Display.
  
  options(9)  = 0;        % Iterate until convergence.
end
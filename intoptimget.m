function [] = intoptimget(options)

% Display parameters setted by an option vector. If no input arguments are 
% supplied display the default option vector.
% 
% Example of Usage
% >> opt = intoptimset();
% >> intoptimget(opt);
% Function Termination Tolerance:            1.000000e-08 
% Box Termination Tolerance:                 1.000000e-08 
% Equality Constraints Tolerance:            1.000000e-08 
% Initial guess for global minimum:          Inf 
% Number of Iterations Allowed:              1000 
% Number of Function evaluations allowed:    10000 
% Number of Constraints evaluations allowed: 10000 
% Level of Display:                          10 
% Iterate Until Convergence:                 0 
%
% This is an INTLAB file. It requires to have INTLAB installed under
% MATLAB to function properly.
%
% Arguments with * are needed. 
%
% Argument    i/o   description
%   options   i    A vector with actual settings.
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

  if nargin == 0
   options = intoptimset();
  end

  fprintf('Function Termination Tolerance:            %e \n', options(1));
  fprintf('Box Termination Tolerance:                 %e \n', options(2));
  fprintf('Equality Constraints Tolerance:            %e \n', options(3));
  fprintf('Initial guess for global minimum:          %e \n', options(4));
  fprintf('Number of Iterations Allowed:              %d \n', options(5));
  fprintf('Number of Function evaluations allowed:    %d \n', options(6));
  fprintf('Number of Constraints evaluations allowed: %d \n', options(7));
  fprintf('Level of Display:                          %d \n', options(8));
  fprintf('Iterate Until Convergence:                 %d \n', options(9));
end
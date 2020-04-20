function y = intdiv(a, b)

% Extended interval division. This implementation follows the technical 
% report
% Inclusion Isotone Interval Arithmetic. A toolbox update. Ratz D., 1996.
% Available at www.rz.uni-karlsruhe.de/~iam/html/reports/rep9605.ps.gz
% 
% Example of Usage
% >> a = infsup(2, 3);
% >> b = infsup(-4, 5);
% >> c = intdiv(a, b)
% [- Inf,-0.5000] [0.4000, Inf] 
%
% This is an INTLAB file. It requires to have INTLAB installed under
% MATLAB to function properly.
%
% Arguments with * are needed.
%
% Argument    i/o   description
%   a           i    numerator.
%   b           i    denominator.
%   y           o    a / b. An interval, an interval vetor or a NaN.
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

  if ~in(0,b)   
   y = a / b;   
   return;
  end

  if in(0, a) && in(0, b)
   y = infsup(-inf, inf);   
   return;
  end
  
  if a.sup < 0 && b.sup == 0   
   y = infsup(a.sup / b.inf, inf);   
   return;
  end
  
  if a.sup < 0 && b.inf == 0   
   y = infsup(-inf, a.sup / b.sup);
   return;
  end
  
  if a.inf > 0 && b.sup == 0   
   y = infsup(-inf, a.inf / b.inf);
   return;    
  end
  
  if a.inf > 0 && b.inf == 0
   y = infsup(a.inf / b.sup, inf);
   return;      
  end
  
  if a.sup < 0 && (b.inf < 0 && 0 < b.sup)
   y1 = infsup(-inf, a.sup / b.sup);
   y2 = infsup(a.sup / b.inf, inf); 
   y = [y1 y2];
   return;
  end

  if a.inf > 0 && (b.inf < 0 && 0 < b.sup)
   y1 = infsup(-inf, a.inf / b.inf);
   y2 = infsup(a.inf / b.sup, inf);
   y = [y1 y2];
   return;
  end 
  
  if ~in(0,a) && b == 0
   y = NaN;
   return
  end
end
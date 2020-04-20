function [x, ival, exitflag, output, L] = intsolve(f, x0, options)

% Find all zeros of a nonlinear system of equations inside an
% initial box x0.
%
% Example of Usage: 
% >> f = @(x) [x(1)^2+ x(2)^2 - 1, x(1)^2 - x(2)]
% >> x0 = [infsup(-5,5), infsup(-5,5)];
% >> x  = intsolve(f, x0)
%     ite 	 evf 	 evj 	 evinv 	 nsol 	 nlist 	 completed 
%
%      0 	 0 	     0 	     0 	     0 	     1	     0.000
%
%      10 	 87 	 10 	 6 	     1  	 3	     0.248
%
%      20 	 112 	 20 	 11 	 1 	     3	     0.436
%
%      30 	 199 	 30 	 17 	 2 	     3	     0.595
%
%      40 	 224 	 40 	 22 	 2 	     3	     0.901
% intval x = 
%      [0.7861, 0.7862] [0.6180, 0.6181] 
%      [-0.7862, -0.7861] [0.6180, 0.6181] 
%
% This example shows that our algorithm will enclosure all solutions for the
% problem min f(x) inside the box x0. Our algorithm can give some boxes as
% answer that do not contain any global minima.
%
% This is an INTLAB file. It requires to have INTLAB installed under
% MATLAB to function properly.
%
% This method uses automatic differentiation to evaluate the jacobian 
% matrix in a box X. The user just need to supply the objective function
% but he/she must to be sure that the function f is differentiable over
% the initial box x0.
%
% The algorithm presented here is a modification of the Hansen-Greenberg 
% algorithm proposed in 
% An Interval Newton Method. Hansen E., Greenberg R. I., 1983. Applied
% Mathematics and Computation. 12: 89-98.
%
% Users can set parameters such as box and function tolerance, maximum
% number of function and outer iterates and the level of display. See
% the intoptimset function.
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
%   x          o    Vector with interval Solutions.
%
%   ival       o    Vector with function evaluated at each solution.
%
%   exitflag   o    Integer identifying the reason the algorithm terminated
%                     0  Algorithm finished with success. 
%                     1  Algorithm finished with maximum iterates.
%                     2  Algorithm finished with maximum function
%                     evaluations
%
%   output     o    Information vector
%                     output(1) = number of outer iterations.
%                     output(2) = number of function evaluations.
%                     output(3) = number of jacobian evaluations.
%                     output(4) = number of matrix invertions.
%                     output(5) = number of non null gauss seidel steps
%                     realized. See the paper above.
%                     output(6) = number of null gauss seidel steps
%                     realized. See the paper above.
%                     output(7) = number of solutions.
%                     output(8) = percentage of domain reached by the
%                     method.
%
%   L          o    A list with unreached boxes. It must be empty if
%                   exitflag equals 0.
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
  
  epsF     = opt(1);
  epsX     = opt(2);
  maxIte   = opt(5);
  maxFunEv = opt(6);
  display  = opt(8);
  
  if opt(9) == 1
   maxIte   = inf;
   maxFunEv = inf;
   L        = intval(zeros(10^6, length(x0)));
  else
   L    = intval(zeros(2 * maxIte, length(x0)));
  end

  x          = intval(zeros(10^6, length(x0)));
  ival       = intval(zeros(10^6, length(x0)));
  exitflag   = 0;
  nlist      = 1;
  L(nlist,:) = x0;  
  total      = prod(diam(x0));
  ite        = 0;
  evf        = 0;
  evj        = 0;
  evi        = 0;
  evnngs     = 0;
  evngs      = 0;
  nsol       = 0;
  pct        = 0;
  
  while nlist > 0
    if display > 0
     if mod(ite, 20 * display) == 0
      clc;
      fprintf('\n ite \t evf \t evj \t evinv \t nsol \t nlist \t completed \n');
     end
     
     if mod(ite, display) == 0
      fprintf('\n %d \t %d \t %d \t %d \t %d \t %d\t %.3f\n', ite, evf, evj,...
                                                      evi, nsol, nlist, pct);
     end
    end

    if ite >= maxIte
     exitflag = 1;
     break;
    end

    ite = ite + 1;   
     
    X          = L(nlist,:);
    L(nlist,:) = intval(NaN);
    nlist      = nlist - 1;    
    vx         = prod(diam(X));
    
    if any(isNaN(X))
     continue;
    end
    
    eF = feval(f, gradientinit(X));
    evf = evf + 1;
    evj = evj + 1;
    
    if evf > maxFunEv
     exitflag = 2;
     break;
    end
    
    if ~all(in(0,eF.x))   
     pct = pct + vx / total;
     continue;
    end
    
    if max(diam(X)) < epsX && max(abss(eF.x)) < epsF
     pct          = pct + vx / total;
     nsol         = nsol + 1;
     x(nsol,:)    = X;
     ival(nsol,:) = eF.x;     
     continue;
    end
    
    A = mid(eF(:).dx);
    
    if issparse(A)
     s = svds(A);
    else
     s = svd(A);
    end
    
    if min(s) < 10^-6     
     [wid, ind]    = max(diam(X));
     nlist         = nlist + 1;
     L(nlist,:)    = X;
     L(nlist, ind) = infsup(X(ind).inf, mid(X(ind)));
     nlist         = nlist + 1;
     L(nlist,:)    = X;
     L(nlist,ind)  = infsup(mid(X(ind)), X(ind).sup);
     continue;
    end
    
    B       = inv(A);
    evi     = evi + 1;
    M       = B * eF(:).dx;
    XNew    = X;
    mx      = mid(XNew);
    fx      = feval(f, intval(mx));
    evf     = evf + 1;   
    
    if evf > maxFunEv
     exitflag = 2;
     break;
    end
    
    [r,c] = size(fx);
    
    if c > r
     fx = fx';
    end

    b       = -1 * B * fx;    
    satisfy = 1;
    findsol = 0;
    
    while satisfy
      [XNew, info]  = nonnullgs (X, mx, b, M);      
      evnngs        = evnngs + 1;

      if info == -1
       pct     = pct + vx / total;
       satisfy = 0;
       continue;
      end

      eF = feval(f, XNew);
      evf = evf + 1;

      if evf > maxFunEv
       exitflag   = 2;
       satisfy = 0;
       continue;
      end

      if ~all(in(0,eF))
       pct     = pct + vx / total;
       satisfy = 0;
       findsol = 1;
       continue;
      end
      
      if max(diam(XNew)) < epsX && max(abss(eF)) < epsF
       pct          = pct + vx / total;
       nsol         = nsol + 1;
       x(nsol,:)    = XNew;
       ival(nsol,:) = eF;
       
       findsol = 1;
       satisfy = 0;
       continue;
      end
      
      if ~all(XNew == X) 
       X  = XNew;
       mx = mid(XNew);
       fx      = feval(f, intval(mx));
       evf     = evf + 1;
    
       if evf > maxFunEv
        exitflag   = 2;
        satisfy = 0;
        continue;
       end
    
       [r,c] = size(fx);
    
       if c > r
        fx = fx';
       end
       
       b       = -1 * B * fx;
       satisfy = 1;
      else
       satisfy = 0;
      end
    end
          
    if exitflag == 2
     break;
    end
    
    if info == -1 || findsol == 1
     continue;
    end

    mx  = mid(XNew);
    fx  = feval(f, intval(mx));
    evf = evf + 1;
    
    if evf > maxFunEv
     exitflag = 2;
     break;
    end
    
    [r,c] = size(fx);
    
    if c > r
     fx = fx';
    end
    
    b                 = -1 * B * fx;
    [XNew, gap, info] = nullgs(XNew, mx, b, M);
    evngs             = evngs + 1;
    
    if info == -1
     pct = pct + vx / total;
     continue;
    end
    
    eF = feval(f, XNew);
    
    if evf > maxFunEv
     exitflag = 2;
     break;
    end
    
    if ~all(in(0,eF))
     pct     = pct + vx / total;
     continue;
    end
    
    if max(diam(XNew)) < epsX && max(abss(eF)) < epsF
     pct          = pct + vx / total;   
     nsol         = nsol + 1;   
     x(nsol,:)    = XNew;
     ival(nsol,:) = eF;    
     continue;
    end
    
    if info > 0     
     nlist         = nlist + 1;
     L(nlist,:)    = XNew;
     L(nlist,info) = gap(1);
     vg1           = prod(diam(L(nlist,:)));
     nlist         = nlist + 1;
     L(nlist,:)    = XNew;
     L(nlist,info) = gap(2);
     vg2           = prod(diam(L(nlist,:)));     
     vf            = vx - (vg1 + vg2);
     pct           = pct + vf / total;
     continue;
    else        
     [mdiam, ind] = max(diam(XNew));
     nlist        = nlist + 1;     
     L(nlist,:)   = XNew;     
     L(nlist,ind) = infsup(XNew(ind).inf, mid(XNew(ind)));
     nlist        = nlist + 1;
     L(nlist,:)   = XNew;     
     L(nlist,ind) = infsup(mid(XNew(ind)), XNew(ind).sup);
     vnx          = prod(diam(XNew));   
     pct          = pct + (vx - vnx) / total;   
    end
  end

  if nlist > 0
   L = L(1:nlist,:);
  else
   L = [];
  end
  
  if nsol > 0
   x    = x(1:nsol,:);
   ival = ival(1:nsol,:);
  else
   x    = [];
   ival = [];
  end
  
  if nargout >= 4
   output = [ite evf evj evi evnngs evngs nsol pct];
  end
end
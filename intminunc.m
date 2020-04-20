function [x, ival, exitflag, output, L] = intminunc(f, x0, options)
  
% Find the global minima of a function f inside the box x0. The user must
% to assure that the solutions of f lies in the interior of x0. This
% function may fail to find a minimum in the border of x0.
% 
% Example of Usage:
% >> f  = @(x) (x(1)^2+ x(2)^2 - 1)^2 +  (x(1)^2 - x(2))^2
% >> x0 = [infsup(-10,10), infsup(-10,10)];
% >> x  = intminunc(f, x0)
%  ite 	 evf  evg 	 evi 	 evh 	 min 	     nsol 	 nlist 	 completed 
%
%  0 	 1 	  0 	 0 	     0 	     Inf 	     0 	     1 	     0.000
%
%  10 	 71   50 	 10 	 10 	 1.000e+00 	 0 	     5 	     0.625
%
% 20 	 141  100 	 20 	 20 	 1.000e+00 	 0 	     5 	     0.938
%
% 30 	 211  150 	 30 	 30 	 1.028e-01 	 0 	     6 	     0.986
%
% 40 	 287  206 	 40 	 40 	 1.028e-01 	 0 	     13 	 0.991
% intval x = 
%  [-0.7862,-0.7861] [0.6180, 0.6181] 
%  [0.7861, 0.7862] [0.6180, 0.6181] 
%
% This example shows that our algorithm will enclosure all solutions for the
% problem min f(x) inside the box x0. Our algorithm can give some boxes as
% answer that do not contain any global minima.
%
% This is an INTLAB file. It requires to have INTLAB installed under
% MATLAB to function properly.
%
% This method uses automatic differentiation to evaluate gradients and 
% hessians on boxes. The user need to supply only the objective function
% but he/she must to be sure that the function f is differentiable over
% the initial box x0.
%
% Users can set parameters such as box and function tolerance, maximum
% number of function, outer iterates, the initial estimative of the global 
% minimum and the level of display. See the intoptimset function.
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
%                     output(3) = number of gradient evaluations.
%                     output(4) = number of matrix invertions.
%                     output(5) = number of hessian evaluations.
%                     output(6) = global minimum estimate.
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
  minglob  = opt(4);
  maxIte   = opt(5);
  maxFunEv = opt(6);
  display  = opt(8);
  
  if opt(9) == 1
   maxIte   = inf;
   maxFunEv = inf;
   L        = intval(zeros(10^6, length(x0)));
   minboxes = intval(zeros(10^6, 1));
  else
   L        = intval(zeros(5 * maxIte, length(x0)));
   minboxes = intval(zeros(5 * maxIte, 1));
  end
  
  x                 = intval(zeros(10^6, length(x0)));
  ival              = intval(zeros(10^6, 1));
  exitflag          = 0;
  nlist             = 1;
  L(nlist,:)        = x0;
  minboxes(nlist,:) = feval(f, x0);
  total             = prod(diam(x0));
  ite               = 0;
  evf               = 1;
  evg               = 0;
  evi               = 0;
  evh               = 0;
  nsol              = 0;
  pct               = 0;
  
  while nlist > 0
    if display > 0
     if mod(ite, 20 * display) == 0
      clc;
      fprintf('\n ite \t evf \t evg \t evi \t evh \t min \t\t nsol \t nlist \t completed \n');
     end
     
     if mod(ite, display) == 0
      fprintf('\n %d \t %d \t %d \t %d \t %d \t %.3e \t %d \t %d \t %.3f\n', ite,...
                                   evf, evg, evi, evh, minglob, nsol, nlist, pct);
     end
    end
       
    if ite >= maxIte
     exitflag = 1;
     break;
    end
    
    ite = ite + 1;    
    i   = 1;
    
    while i <= nlist
      if minboxes(i).inf > minglob
       for j = i: nlist
         pct         = pct + prod(diam(L(j,:))) / total;
         minboxes(j) = intval(NaN);
         L(j,:)      = intval(NaN);
         nlist       = nlist - 1;
       end
       break;
      else
       i = i + 1;
      end
    end
    
    if nlist == 0
     continue;
    end
    
    X  = L(1,:);
    eF = minboxes(1);
    
    if nlist > 1
     L(1:nlist-1,:)      = L(2:nlist,:);
     minboxes(1:nlist-1) = minboxes(2:nlist);
    end
    
    L(nlist,:)      = intval(NaN);
    minboxes(nlist) = intval(NaN);
    nlist           = nlist - 1;
    vx              = prod(diam(X));
    eH              = feval(f, hessianinit(X));
    evg             = evg + 1;
    evh             = evh + 1;
    
    if ~all(in(0,eH(:).dx))        
     pct = pct + vx / total;
     continue;
    end
    
    if any(diag(eH(:).hx.sup) < 0)
     pct = pct + vx / total;
     continue;
    end
    
    mx  = mid(X);
    fmx = feval(f, gradientinit(intval(mx)));
    evf = evf + 1;
    evg = evg + 1;
    
    if evf > maxFunEv
     exitflag = 2;
     break;
    end
    
    if fmx.x.sup < minglob
     minglob = fmx.x.sup;
     i       = 1;
     
     while i <= nsol
       if ival(i).inf > minglob
        for j = i: nsol
          ival(j) = intval(NaN);
          x(j,:)  = intval(NaN);
          nsol    = nsol - 1;
        end
        break;
       else
        i = i + 1;
       end
     end
     
     i = 1;
     
     while i <= nlist
       if minboxes(i).inf > minglob
        for j = i: nlist
          pct         = pct + prod(diam(L(j,:))) / total;
          minboxes(j) = intval(NaN);
          L(j,:)      = intval(NaN);
          nlist       = nlist - 1;
        end
        break;
       else
        i = i + 1;
       end
     end
    end

    if max(diam(X)) < epsX && max(abss(eH(:).dx)) < epsF
     pct          = pct + vx / total;
     nsol         = nsol + 1;
     x(nsol,:)    = X;
     ival(nsol)   = eF;
     [x,ival]     = sortlists(x, ival, nsol);
     continue;
    end

    A = mid(eH(:).hx);
       
    if issparse(A)
     s = svds(A);
    else
     s = svd(A);
    end 
      
    if min(s) < 10^-6
     [wid, ind]      = max(diam(X));
     nlist           = nlist + 1;
     L(nlist,:)      = X;
     L(nlist,ind)    = infsup(X(ind).inf, mid(X(ind)));
     minboxes(nlist) = f(L(nlist,:));
     [L, minboxes]   = sortlists(L, minboxes, nlist);
     nlist           = nlist + 1;
     L(nlist,:)      = X;
     L(nlist,ind)    = infsup(mid(X(ind)), X(ind).sup);
     minboxes(nlist) = f(L(nlist,:));
     [L, minboxes]   = sortlists(L, minboxes, nlist);
     evf             = evf + 2;
     
     if evf > maxFunEv
      exitflag = 2;
      break;
     end
     continue;
    end
    
    B       = inv(A);
    evi     = evi + 1;
    M       = B * eH(:).hx;
    gfx     = fmx(:).dx;
    [r,c]   = size(gfx);
    
    if c > r
     gfx = gfx';
    end
    
    b       = -1 * B * gfx;
    satisfy = 1;
    findsol = 0;
    
    while satisfy 
      [XNew, info]  = nonnullgs (X, mx, b, M);
      
      if info == -1
       pct     = pct + vx / total;
       satisfy = 0;
       continue;
      end

      eG  = feval(f, gradientinit(XNew));
      evf = evf + 1;
      evg = evg + 1;

      if evf > maxFunEv
       exitflag = 2;
       satisfy  = 0;
       continue;
      end
      
      if ~all(in(0,eG(:).dx))
       pct     = pct + vx / total;
       satisfy = 0;
       findsol = 1;
       continue;
      end
      
      if eG.x.inf > minglob
       pct     = pct + vx / total;
       satisfy = 0;
       findsol = 1;
       continue;
      end
      
      mx  = mid(XNew);
      fmx = feval(f, gradientinit(intval(mx)));
      evf = evf + 1;
      evg = evg + 1;
    
      if evf > maxFunEv
       exitflag = 2;
       satisfy  = 0;
       continue;
      end
      
      if fmx.x.sup < minglob
       minglob = fmx.x.sup;
       i       = 1;   
       
       while i <= nsol
         if ival(i).inf > minglob
          for j = i: nsol
            ival(j) = intval(NaN); 
            x(j,:)  = intval(NaN);
            nsol    = nsol - 1;
          end
          break;
         else
          i = i + 1;
         end
       end
     
       i = 1;
     
       while i <= nlist
         if minboxes(i).inf > minglob
          for j = i: nlist   
            pct         = pct + prod(diam(L(j,:))) / total;
            minboxes(j) = intval(NaN);
            L(j,:)      = intval(NaN);
            nlist       = nlist - 1;
          end
          break;
         else
          i = i + 1;
         end
       end
      end
      
      if max(diam(XNew)) < epsX && max(abss(eG(:).dx)) < epsF
       pct          = pct + vx / total;
       nsol         = nsol + 1;   
       x(nsol,:)    = XNew;
       ival(nsol)   = eG.x;
       [x, ival]    = sortlists(x, ival, nsol);
       findsol      = 1;
       satisfy      = 0;
       continue;
      end
      
      if ~all(XNew == X) 
       X     = XNew;   
       gfx   = fmx(:).dx;
       [r,c] = size(gfx);
    
       if c > r
        gfx = gfx';
       end
       
       b = -1 * B * gfx;       
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

    gfx   = fmx(:).dx;
    [r,c] = size(gfx);
    
    if c > r
     gfx = gfx';
    end
       
    b                 = -1 * B * gfx;
    [XNew, gap, info] = nullgs(XNew, mx, b, M);

    if info == -1
     pct = pct + vx / total;
     continue;
    end
    
    eG  = feval(f, gradientinit(XNew));
    evf = evf + 1;
    evg = evg + 1;

    if evf > maxFunEv
     exitflag = 2;
     break;
    end
   
    if ~all(in(0,eG(:).dx))
     pct = pct + vx / total;
     continue;
    end

    if eG.x.inf > minglob
     pct = pct + vx / total;
     continue;
    end
      
    mx  = mid(XNew);
    fmx = feval(f, intval(mx));
    evf = evf + 1;
  
    if evf > maxFunEv
     exitflag = 2;
     break;
    end
        
    if fmx.sup < minglob
     minglob = fmx.sup;
     i       = 1;
     
     while i <= nsol
       if ival(i).inf > minglob
        for j = i: nsol
          ival(j) = intval(NaN); 
          x(j,:)  = intval(NaN);
          nsol    = nsol - 1;
        end
        break;
       else
        i = i + 1;
       end
     end
     
     i = 1;
     
     while i <= nlist
       if minboxes(i).inf > minglob
        for j = i: nlist   
          pct         = pct + prod(diam(L(j,:))) / total;
          minboxes(j) = intval(NaN);
          L(j,:)      = intval(NaN);
          nlist       = nlist - 1;
        end
        break;
       else
        i = i + 1;
       end
     end
    end

    if max(diam(XNew)) < epsX && max(abss(eG(:).dx)) < epsF        
     pct        = pct + vx / total;
     nsol       = nsol + 1;
     x(nsol,:)  = XNew;
     ival(nsol) = eG.x;
     [x, ival]  = sortlists(x, ival, nsol);
     continue;
    end
    
    if info > 0    
     nlist           = nlist + 1;
     L(nlist,:)      = XNew;
     L(nlist,info)   = gap(1);
     minboxes(nlist) = f(L(nlist,:));
     vg1             = prod(diam(L(nlist,:)));
     [L, minboxes]   = sortlists(L, minboxes, nlist);     
     nlist           = nlist + 1;
     L(nlist,:)      = XNew;
     L(nlist,info)   = gap(2);
     minboxes(nlist) = f(L(nlist,:));
     vg2             = prod(diam(L(nlist,:)));
     [L, minboxes]   = sortlists(L, minboxes, nlist);     
     vf              = vx - (vg1 + vg2);
     pct             = pct + vf / total;
     evf             = evf + 2;
     
     if evf > maxFunEv
      exitflag = 2;
      break;
     end     
     continue;
     
    else     
     [mdiam, ind]    = max(diam(XNew));     
     nlist           = nlist + 1;
     L(nlist,:)      = XNew;
     L(nlist,ind)    = infsup(XNew(ind).inf, mid(XNew(ind)));     
     minboxes(nlist) = f(L(nlist,:));     
     [L, minboxes]   = sortlists(L, minboxes, nlist);     
     nlist           = nlist + 1;
     L(nlist,:)      = XNew;
     L(nlist,ind)    = infsup(mid(XNew(ind)), XNew(ind).sup);
     minboxes(nlist) = f(L(nlist,:));
     [L, minboxes]   = sortlists(L, minboxes, nlist);     
     vnx             = prod(diam(XNew));
     pct             = pct + (vx - vnx) / total;
     evf             = evf + 2;
     
     if evf > maxFunEv
      exitflag = 2;
      break;
     end
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
   output = [ite evf evg evi evh minglob nsol pct];
  end
end

function [A, B] = sortlists(L1, L2, n)

  A = L1;
  B = L2;
  
  for i = 1: n - 1
    if L2(n).inf < L2(i).inf
     tmp        = B(i:n-1);
     B(i)       = B(n);
     B(i+1:n)   = tmp;
     tmp        = A(i:n-1,:);
     A(i,:)     = A(n,:);
     A(i+1:n,:) = tmp;
     break
    end
  end
end
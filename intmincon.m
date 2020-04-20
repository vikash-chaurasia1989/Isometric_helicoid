function [x, ival, exitflag, output, L] = intmincon(f, x0, nonlcon, options)

% Find the global minima of a function f inside a box x0 subject to
% equality constraints nonlcon. The user must to assure that the solutions 
% of f lies in the interior of x0. This function may fail to find a minimum
% in the border of x0.
% 
% Example of Usage:
% >> f  = @(x) x(1) + x(2);
% >> c  = @(x) x(1)^2 + x(2)^2 - 1;
% >> x0 = [infsup(-1,1), infsup(-1,1)];
% >> x  = intmincon(f, x0, c)
% ite 	 evf 	 evg 	 evh 	 evc 	 evcg 	 evch 	 min 	 nsol 	 nlist 	 completed 
%
% 0 	 1 	     0 	     0 	     0 	     0 	     0 	     Inf 	 0 	     1 	     0.000 
%
% 10 	 59 	 49 	 5 	     49 	 49 	 5 	     Inf 	 0 	     7 	     0.051 
%
% 20 	 132 	 120 	 11 	 122 	 122 	 11 	 Inf 	 1 	     5 	     0.188 
%
% 30 	 160 	 144 	 16 	 148 	 148 	 16 	 Inf 	 1 	     5 	     0.375 
%
% 40 	 182 	 166 	 22 	 172 	 172 	 22 	 Inf 	 1 	     3 	     0.625 
%
% 50 	 210 	 190 	 27 	 198 	 198 	 27 	 Inf 	 1 	     3 	     0.812 
%
% 60 	 266 	 239 	 35 	 249 	 249 	 35 	 Inf 	 1 	     3 	     0.945 
% intval x = 
%   -0.7071   -0.7071
%    0.7071    0.7071
%
% This example shows that our algorithm will enclosure all solutions for the
% problem min f(x) s.t c(x) = 0 inside the box x0. Our algorithm can give
% some boxes as answer that do not contain any global minima. The second
% box in the example above show this.

% This is an INTLAB file. It requires to have INTLAB installed under
% MATLAB to function properly.
%
% This method uses automatic differentiation to evaluate gradients and 
% hessians on boxes. The user need to supply only the objective function
% and a equality constraints function but he/she must to be sure that the 
% function f is differentiable over the initial box x0.
%
% Users can set parameters such as box, objective function and constraint 
% function tolerance. Maximum number of objective function, constraint
% functions and outer iterates. The initial estimative of the global mini
% mum and the level of display. See the intoptimset function.
%
% Arguments with * are needed.
%
% Argument    i/o   description
%   f*         i    Objective Function.
%
%   X0*        i    Initial Interval(s). A row vector containing the
%                   initial box
%
%   nonlcon*   i    Equality constraints function. 
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
%                        evaluations.
%                     3  Algorithm finished with maximum constraint 
%                        evaluations.
%
%   output     o    Information vector
%                     output(1) = number of outer iterations.
%                     output(2) = number of function evaluations.
%                     output(3) = number of gradient evaluations.
%                     output(4) = number of hessian evaluations.
%                     output(5) = number of constraint evaluations.
%                     output(6) = number of constraints' gradient evaluations.
%                     output(7) = number of hessian's constraints evaluations.
%                     output(8) = global minimum estimate.
%                     output(9) = number of solutions.
%                     output(10)= percentage of domain reached by the
%                     method.
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
      
  if nargin < 4
   opt = intoptimset();
   
   if nargin < 3
    error('Missing Arguments. You must to specify an objective function, an initial box and a constraint function.');
   end
   
  else
   opt = options;
  end
  
  epsF     = opt(1);
  epsX     = opt(2);
  epsEC    = opt(3);
  minglob  = opt(4);
  maxIte   = opt(5);
  maxFunEv = opt(6);
  maxConEv = opt(7);
  display  = opt(8);
  
  if opt(9) == 1
   maxIte   = inf;
   maxFunEv = inf;
   maxConEv = inf;
   L        = intval(zeros(10^6, length(x0)));
   minboxes = intval(zeros(10^6, 1));
  else
   L        = intval(zeros(20 * maxIte, length(x0)));
   minboxes = intval(zeros(20    * maxIte, 1));
  end
  
  n                 = length(x0);
  x                 = intval(zeros(10^6, length(x0)));
  ival              = intval(zeros(10^6, 1));
  exitflag          = 0;
  nlist             = 1;
  L(nlist,:)        = x0;
  minboxes(nlist,:) = feval(f, x0);
  total             = prod(diam(x0));
  ite               = 0;
  evf               = 1;
  evc               = 0;
  evg               = 0;
  evh               = 0;
  evcg              = 0;
  evch              = 0;
  nsol              = 0;
  pct               = 0;
    
  while nlist > 0
    if display > 0
     if mod(ite, 20 * display) == 0
      clc;
      fprintf('\n ite \t evf \t evg \t evh \t evc \t evcg \t evch \t min \t nsol \t nlist \t completed \n');
     end
     
     if mod(ite, display) == 0
      fprintf('\n %d \t %d \t %d \t %d \t %d \t %d \t %d \t %.5e \t %d \t %d \t %.3f \n', ite,...
                                   evf, evg, evh, evc, evcg, evch, minglob, nsol, nlist, pct);
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
           
    if nlist > 1
     L(1:nlist-1,:)      = L(2:nlist,:);
     minboxes(1:nlist-1) = minboxes(2:nlist);
    end
    
    L(nlist,:)      = intval(NaN);
    minboxes(nlist) = intval(NaN);
    nlist           = nlist - 1;
    vx              = prod(diam(X));
    eC              = feval(nonlcon, gradientinit(X));
    evc             = evc + 1;
    evcg            = evcg + 1;
    
    if evc > maxConEv
     exitflag = 3;
     break;
    end
    
    if ~all(in(0, eC.x))
     pct = pct + vx / total;
     continue;
    end
    
    eF    = feval(f, gradientinit(X));
    evg   = evg + 1;
    A     = eC(:).dx';
    [r,c] = size(eF(:).dx);
    
    if c > r
     b = eF(:).dx';
    else
     b = eF(:).dx;
    end
    
    lambda = verifylss(A,b);
    
    if any(isNaN(lambda))
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
    
    m            = length(lambda);
    [flag, nmin] = verify_existence(f, nonlcon, X, A', x0);
    
    if flag == 1
     if nmin < minglob
      minglob = nmin;
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
    end
        
    XNew      = [X,lambda'];
    [eL, eDL] = evallag(f,nonlcon, XNew(1:n), XNew(n+1:n+m));
    evf       = evf + 1;
    
    if evf > maxFunEv
     exitflag = 2;
     break;
    end
    
    evg = evg + 1;
    evh = evh + 1;
    evc = evc + 1;
    
    if evc > maxConEv
     exitflag = 3;
     break;
    end
    
    evcg = evcg + 1;
    evch = evch + 1;
    A    = mid(eDL);
       
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
    
    B   = inv(A);
    M   = B * eDL;
    mx  = mid(XNew);
    gfx = evallag(f, nonlcon, intval(mx(1:n)), intval(mx(n+1:n+m)));
    evf = evf + 1;
    
    if evf > maxFunEv
     exitflag = 2;
     break;
    end
     
    evg = evg + 1;
    evc = evc + 1;
    
    if evc > maxConEv
     exitflag = 3;
     break;
    end
    
    evcg    = evcg + 1;
    b       = -1 * B * gfx;
    satisfy = 1;
    findsol = 0;
        
    while satisfy
      [XN, info]  = nonnullgs (XNew, mx, b, M);
     
      if info == -1
       pct     = pct + vx / total;
       satisfy = 0;
       continue;
      end

      mx  = mid(XN);
      eF  = feval(f, XN(1:n));
      eL  = evallag(f,nonlcon, XN(1:n), XN(n+1:n+m));
      evf = evf + 2;
      evg = evg + 1;

      if evf > maxFunEv
       exitflag = 2;
       satisfy  = 0;
       continue;
      end
      
      evc  = evc + 1;
      evcg = evcg + 1;

      if evc > maxConEv
       exitflag = 3;
       satisfy  = 0;
       continue;
      end
      
      if eF.inf > minglob
       pct     = pct + vx / total;
       satisfy = 0;
       findsol = 1;
       continue;
      end
      
      if ~all(in(0,eL))
       pct     = pct + vx / total;
       satisfy = 0;
       findsol = 1;
       continue;
      end
      
      if max(diam(XN(1:n))) < epsX && max(abss(eL(1:n))) < epsF && max(abss(eL(n+1:n+m))) < epsEC
       pct          = pct + vx / total;   
       nsol         = nsol + 1;   
       x(nsol,:)    = XN(1:n);
       ival(nsol)   = eF;
       [x, ival]    = sortlists(x, ival, nsol);
       findsol      = 1;
       satisfy      = 0;
       continue;
      end
            
      if ~all(XN == XNew) 
       XNew    = XN;
       gfx     = evallag(f, nonlcon, intval(mx(1:n)), intval(mx(n+1:n+m)));
       evf     = evf + 1;
       evg     = evg + 1;
    
       if evf > maxFunEv
        exitflag = 2;
        satisfy  = 0;
        continue;
       end
    
       evc  = evc + 1;
       evcg = evcg + 1;
    
       if evc > maxConEv
        exitflag = 3;
        satisfy  = 0;
        continue;
       end
    
       b = -1 * B * gfx;
      else
       satisfy = 0;
      end
    end

    if exitflag ~= 0
     nlist      = nlist + 1;
     L(nlist,:) = XNew(1:n);   
     break;
    end
    
    if info == -1 || findsol == 1
     continue;
    end

    gfx = evallag(f, nonlcon, intval(mx(1:n)), intval(mx(n+1:n+m)));
    evf = evf + 1;
    evg = evg + 1;
    
    if evf > maxFunEv
     exitflag   = 2;
     break;
    end
    
    evc  = evc + 1;
    evcg = evcg + 1;
    
    if evc > maxConEv
     exitflag   = 3;
     break;
    end
    
    b               = -1 * B * gfx;
    [XN, gap, info] = nullgs(XN, mx, b, M);

    if info == -1
     pct = pct + vx / total;   
     continue;
    end
    
    eF  = feval(f, XN(1:n));
    eL  = evallag(f,nonlcon,XN(1:n), XN(n+1:n+m));
    evf = evf + 2;
    evg = evg + 1;

    if evf > maxFunEv
     exitflag   = 2;
     break;
    end
      
    evc  = evc + 1;
    evcg = evcg + 1;
      
    if evc > maxConEv
     exitflag   = 3;
     break;
    end
      
    if eF.inf > minglob
     pct = pct + vx / total;   
     continue;
    end
      
    if ~all(in(0,eL))
     pct = pct + vx / total;   
     continue;
    end

    if max(diam(XN(1:n))) < epsX && max(abss(eL(1:n))) < epsF && max(abss(eL(n+1:n+m))) < epsEC
     nsol       = nsol + 1;
     x(nsol,:)  = XNew(1:n);
     ival(nsol) = eF;
     [x, ival]  = sortlists(x, ival, nsol);
     continue;
    end
    
    if info > 0 & info <=n
     nlist           = nlist + 1;
     L(nlist,:)      = XN(1:n);
     L(nlist,info)   = gap(1);
     minboxes(nlist) = f(L(nlist,:));
     vg1             = prod(diam(L(nlist,:)));
     [L, minboxes]   = sortlists(L, minboxes, nlist);
     nlist           = nlist + 1;
     L(nlist,:)      = XN(1:n);
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
     [mdiam, ind]    = max(diam(XN(1:n)));     
     nlist           = nlist + 1;
     L(nlist,:)      = XN(1:n);
     L(nlist,ind)    = infsup(XN(ind).inf, mid(XN(ind)));     
     minboxes(nlist) = f(L(nlist,:));     
     [L, minboxes]   = sortlists(L, minboxes, nlist);     
     nlist           = nlist + 1;
     L(nlist,:)      = XN(1:n);
     L(nlist,ind)    = infsup(mid(XN(ind)), XN(ind).sup);
     minboxes(nlist) = f(L(nlist,:));
     [L, minboxes]   = sortlists(L, minboxes, nlist);
     vnx             = prod(diam(XN(1:n)));
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
   output = [ite evf evg evh evc evcg evch minglob nsol pct];
  end
end

function [flag, nmin] = verify_existence(f, nonlcon, X, J, x0)
    
    [m,n] = size(J);
    flag  = 0;
    nmin  = inf;
    A     = J;
    XNew  = X;
    P     = 1:n;

    for k = 1:m
      %Changing Columns and Variables  
      [maxdist, ind] =  max(max(mig(A(k:m,k:n)),[],1));

      if maxdist <= 0
       return;
      end
     
      ind = ind + k - 1;
    
      if ind ~= k
       tmp       = A(:,k);
       A(:,k)    = A(:,ind);
       A(:,ind)  = tmp;
       tmp       = XNew(k);
       XNew(k)   = XNew(ind);
       XNew(ind) = tmp;
       tmp       = P(k);
       P(k)      = P(ind);
       P(ind)    = tmp;
      end
      
      %Changing Rows
      [maxdist, ind] = max(mig(A(k:m, k)));
      ind = ind + k - 1;
    
      if ind ~= k
       tmp      = A(k,:);
       A(k,:)   = A(ind,:);
       A(ind,:) = tmp;
      end

      %Givens
      for i = k+1: m
        Ratio  = A(i,k) / A(k,k);
        A(k,k) = sqrt(A(k,k)^2 + A(i,k)^2);
        A(i,k) = 0;
       
        for j = k+1: n
          u      = A(k,j);
          A(k,j) = (A(k,j) + A(i,j) * Ratio) / sqrt(1 + Ratio^2);
          A(i,j) = (-u * Ratio + A(i,j)) / sqrt(1 + Ratio^2);
        end
      end
    end
    
    y  = mid(XNew(1:m));
    fy = feval(nonlcon, intval(mid(X)));
    
    [r,c] = size(fy);
    
    if c > r
     fy = fy';
    end
        
    delta = verifylss(A(:,P(1:m)), fy);
    
    if any(isNaN(delta))    
     Y  = inv(mid(A(:,P(1:m))));
     K  = y' - Y * fy + (eye(m) - Y * A(:,P(1:m)))*(XNew(1:m) - y)';
    
     if all(in(K', XNew(1:m)))
      flag = 1;
      eF   = feval(f, XNew(P));
      nmin = eF.sup;  
     end
     return;
    end

    XNew(1:m) = intersect(XNew(1:m) + 2 * delta', x0(P(1:m)));
    y         = mid(XNew(1:m));
    fy        = feval(nonlcon, intval(mid(X)));
    
    [r,c] = size(fy);
    
    if c > r
     fy = fy';
    end
    
    X(P(1:m)) = XNew(1:m);
    eC        = feval(nonlcon, gradientinit(X));
    A         = eC(:).dx;
    Y         = inv(mid(A(:,P(1:m))));
    K         = y' - Y * fy + (eye(m) - Y * A(:,P(1:m)))*(XNew(1:m) - y)';
    
    if all(in(K', XNew(1:m)))
     flag = 1;
     eF   = feval(f, XNew(P));
     nmin = eF.sup;
    end
    return;
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

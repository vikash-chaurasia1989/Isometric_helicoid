function x1 = my_rk4(x0)

global ht bx0 by0 bz0 N N1 bx by bz d p1 branch path bxm1 bym1 bzm1 tau bx1 by1 bz1 bx2 by2 bz2  d21 dsum
smop = {'moving','lowess','lowess','sgolay','rlowess','rloess'};
opid  = 3;


k1 = fun_descent7(0,x0);

k2 = fun_descent7(0,x0+ht*k1/2);

k3 = fun_descent7(0,x0+ht*k2/2);

k4 = fun_descent7(0,x0+ht*k3);


x1 = x0 + ht/6*(k1 + 2*k2 + 2*k3 + k4);


%====== RK4 completed here ===========
%
%
%

%=== x is not on the constrained manifold ===
% we solve the minimization problem min |x2-x1|^2 with G(x2)=0, where G is
% the constraint matrix

bx0(2:N,1)    = x1(1:3:3*N1-2,1);
by0(2:N,1)    = x1(2:3:3*N1-1,1);
bz0(2:N,1)    = x1(3:3:3*N1,1);

bx0(1,1) = 1;
by0(1,1) = 0;
bz0(1,1) = 0;

bx0(N+1,1) = -1;
by0(N+1,1) =  0;
bz0(N+1,1) =  0;

% % 
% bx0 = smooth(bx0,smop{opid});
% by0 = smooth(by0,smop{opid});
% bz0 = smooth(bz0,smop{opid});
% % % 
%
bx0(1,1) = 1;
by0(1,1) = 0;
bz0(1,1) = 0;

bx0(N+1,1) = -1;
by0(N+1,1) =  0;
bz0(N+1,1) =  0;

%
%=== calculate the Frechet distance between bx0 and the solution at the
%previous step
str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*16)) '_step_' num2str(p1-1)   '.txt'] ;

%== previous steps ==
x0 =  load([path str0]);

bxm1(2:N,1)     = x0(1:3:3*N1-2,1)    ;
bym1(2:N,1)     = x0(2:3:3*N1-1,1)    ;
bzm1(2:N,1)     = x0(3:3:3*N1,1)      ;

bxm1(1,1)   =  1;
bym1(1,1)   =  0;
bzm1(1,1)   =  0;

bxm1(N+1,1) = -1;
bym1(N+1,1) =  0;
bzm1(N+1,1) =  0;

 


d =  1*sqrt(sum((bxm1-bx0).^2 + (bym1-by0).^2 + (bzm1-bz0).^2)) ;

dsum = dsum + d;   % 
%dsum = sqrt(sum((bx2-bx0).^2 + (by2-by0).^2 + (bz2-bz0).^2)) ; % distance between the solution from steepest descent and equilibrium solution
%==========================================================================
%

%=== Now searching for a solution at distance d from the previous step on
%the constrained manifold with least energy. Least energy ensures
%smoothness of the solution

%== initial guess

bx = bx0./sqrt(bx0.^2+by0.^2+bz0.^2);
by = by0./sqrt(bx0.^2+by0.^2+bz0.^2);
bz = bz0./sqrt(bx0.^2+by0.^2+bz0.^2);

% 
% bx = bxm1;
% by = bym1;
% bz = bzm1;

%   bx = smooth(bx,smop{opid});
%   by = smooth(by,smop{opid});
%   bz = smooth(bz,smop{opid});
 
%=== Using homotopy for initial guess =====


%=== Homotopy ==

bx = (1-dsum/d21)*bx2 + dsum/d21*bx1;
by = (1-dsum/d21)*by2 + dsum/d21*by1;
bz = (1-dsum/d21)*bz2 + dsum/d21*bz1;
 
bx = bx./sqrt(bx.^2+by.^2+bz.^2);
by = by./sqrt(bx.^2+by.^2+bz.^2);
bz = bz./sqrt(bx.^2+by.^2+bz.^2);

bx(1,1) = 1;
by(1,1) = 0;
bz(1,1) = 0;

bx(N+1,1) = -1;
by(N+1,1) =  0;
bz(N+1,1) =  0;


%=============================

% 
% plot3(bx0,by0,bz0)
% hold on
% plot3(bx,by,bz)
% hold on
% plot3(bxm1,bym1,bzm1)

% initial guess for x2 =====

var_in = zeros(5*N1+4,1);

var_in(1:3:3*N1-2,1) = bx(2:N,1)  ;
var_in(2:3:3*N1-1,1) = by(2:N,1)  ;
var_in(3:3:3*N1,1)   = bz(2:N,1)  ;
% 
% %=== lagrange multipliers for the initial guess using pseudoinverse of the
% %constraint gradient
% 
% 
% [F,J] = fun_pullback2(var_in);
% %
% % %== Constraint gradient
% %
% A = J(3*N1+1:5*N1+4,1:3*N1);
% 
% 
% var_in(3*N1+1:5*N1+4,1) = pinv(A')*F(1:3*N1,1);

%var_in(3*N1+1:5*N1+4,1) = 0.1*ones(2*N1+4,1);

%== solving the minimization problem ===
%     %options.Algorithm   = 'levenberg-marquardt'  ;    %'trust-region-reflective' ;

options             = optimset('Display','iter', 'Algorithm','levenberg-marquardt','Jacobian','on', 'TOlFun',10^(-9),'TOlX',10^(-9),'MaxFunEvals',3000  ) ;
options.MaxIter     = 2000  ;
% [x,fval,exitflag,output,qd1] =  fsolve(@fun_pullback2   , var_in ,options)    ;


var_in2 = var_in(1:3*N1,1);
[x,fval,exitflag,output,qd1] =  fsolve(@fun_pullback4  , var_in2 ,options)    ;

 
 %=== fun_pullback4 has work the best so far. 
 

%
% %
%  %=== the output from the above step may not be smooth. We smothen it and
%  %reuse it as guess to do the minimization again
    bx(2:N,1)     = x(1:3:3*N1-2,1)    ;
    by(2:N,1)     = x(2:3:3*N1-1,1)    ;
    bz(2:N,1)     = x(3:3:3*N1,1)      ;

    bx(1,1) = 1;
    by(1,1) = 0;
    bz(1,1) = 0;

    bx(N+1,1) = -1;
    by(N+1,1) =  0;
    bz(N+1,1) =  0;
    
    
    
%     
%     bx = smooth(bx,smop{opid});
%     by = smooth(by,smop{opid});
%     bz = smooth(bz,smop{opid});
%     
    bx(1,1) = 1;
    by(1,1) = 0;
    bz(1,1) = 0;
    bx(N+1,1) = -1;
    by(N+1,1) =  0;
    bz(N+1,1) =  0;
    

   
var_in2(1:3:3*N1-2,1) = bx(2:N,1)  ;
var_in2(2:3:3*N1-1,1) = by(2:N,1)  ;
var_in2(3:3:3*N1,1)   = bz(2:N,1)  ;


 %  [x,fval,exitflag,output,qd1] =  fsolve(@fun_pullback4  , var_in2 ,options)    ;
 
%  %  %=== the output from the above step may not be smooth. We smothen it and
% %  %reuse it as guess to do the minimization again
%     bx(2:N,1)     = x(1:3:3*N1-2,1)    ;
%     by(2:N,1)     = x(2:3:3*N1-1,1)    ;
%     bz(2:N,1)     = x(3:3:3*N1,1)      ;
% 
%     bx(1,1) = 1;
%     by(1,1) = 0;
%     bz(1,1) = 0;
% 
%     bx(N+1,1) = -1;
%     by(N+1,1) =  0;
%     bz(N+1,1) =  0;
%     
%     
%     
%     
%     bx = smooth(bx,smop{opid});
%     by = smooth(by,smop{opid});
%     bz = smooth(bz,smop{opid});
% %     
%     bx(1,1) = 1;
%     by(1,1) = 0;
%     bz(1,1) = 0;
%     bx(N+1,1) = -1;
%     by(N+1,1) =  0;
%     bz(N+1,1) =  0;
%     
% 
%    
% var_in2(1:3:3*N1-2,1) = bx(2:N,1)  ;
% var_in2(2:3:3*N1-1,1) = by(2:N,1)  ;
% var_in2(3:3:3*N1,1)   = bz(2:N,1)  ;
% 
% 
%    [x,fval,exitflag,output,qd1] =  fsolve(@fun_pullback5  , var_in2 ,options)    ;

% %     bx(2:N+1,1) = .5*(bx(1:N,1) + bx(2:N+1,1));
% %     by(2:N+1,1) = .5*(by(1:N,1) + by(2:N+1,1));
% %     bz(2:N+1,1) = .5*(bz(1:N,1) + bz(2:N+1,1));
% %
%
%     bx = smooth(bx);
%     by = smooth(by);
%     bz = smooth(bz);
%
%     bx(1,1) = 1;
%     by(1,1) = 0;
%     bz(1,1) = 0;
%     bx(N+1,1) = -1;
%     by(N+1,1) =  0;
%     bz(N+1,1) =  0;
%
%
%     x(1:3:3*N1-2,1) =  bx(2:N,1) ;
%     x(2:3:3*N1-1,1) =  by(2:N,1) ;
%     x(3:3:3*N1,1)   =  bz(2:N,1) ;
%
%
%  [x,fval,exitflag,output,qd1] =  fsolve(@fun_pullback2 , x ,options)    ;

%=== replacing the binormal vector b stored in x1 with the binormal on the
%constrained manifold ====




x1(1:3*N1,1) = x(1:3*N1,1);













end
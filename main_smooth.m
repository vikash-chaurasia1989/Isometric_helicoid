clear all 
clc 

%=== post processing the evolution data and searching for smooth solution 

global ht bx0 by0 bz0 N N1 bx by bz d p1 branch path bxm1 bym1 bzm1 tau bx1 by1 bz1 bx2 by2 bz2  d21 dsum 
smop = {'moving','lowess','lowess','sgolay','rlowess','rloess'};
opid  = 3;

parameters5();

for p1 = 405:505


%previous step
%previous step
str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*16)) '_step_' num2str(p1-1)   '.txt'] ;

%== previous steps ==
x0 =  load([path str0]);

bxm1(2:N,1)     = x0(1:3:3*N1-2,1)    ;
bym1(2:N,1)     = x0(2:3:3*N1-1,1)    ;
bzm1(2:N,1)     = x0(3:3:3*N1,1)      ;

bxm1(1,1) = 1;
bym1(1,1) = 0;
bzm1(1,1) = 0;

bxm1(N+1,1) = -1;
bym1(N+1,1) =  0;
bzm1(N+1,1) =  0;

 
%==== current step which serves as initial guess for finding smooth
%solution 


str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*16)) '_step_' num2str(p1)   '.txt'] ;

%== previous steps ==
x1 =  load([path str0]);

bx(2:N,1)     = x1(1:3:3*N1-2,1)    ;
by(2:N,1)     = x1(2:3:3*N1-1,1)    ;
bz(2:N,1)     = x1(3:3:3*N1,1)      ;

bx(1,1) = 1;
by(1,1) = 0;
bz(1,1) = 0;

bx(N+1,1) = -1;
by(N+1,1) =  0;
bz(N+1,1) =  0;

d =  sqrt(sum((bxm1-bx).^2 + (bym1-by).^2 + (bzm1-bz).^2)) ;
% initial guess for x2 =====

var_in = zeros(5*N1+5,1);

var_in(1:3:3*N1-2,1) = bx(2:N,1)  ;
var_in(2:3:3*N1-1,1) = by(2:N,1)  ;
var_in(3:3:3*N1,1)   = bz(2:N,1)  ;

%=== lagrange multipliers for the initial guess using pseudoinverse of the
%constraint gradient


[F,J] = fun_pullback6(var_in);
%
% %== Constraint gradient
%
A = J(3*N1+1:5*N1+5,1:3*N1);
var_in(3*N1+1:5*N1+5,1) = pinv(A')*F(1:3*N1,1);


%== solving the minimization problem ===
%     %options.Algorithm   = 'levenberg-marquardt'  ;    %'trust-region-reflective' ;

options             = optimset('Display','iter', 'Algorithm','levenberg-marquardt','Jacobian','on', 'TOlFun',10^(-9),'TOlX',10^(-9),'MaxFunEvals',3000  ) ;
options.MaxIter     = 2000  ;
[x,fval,exitflag,output,qd1] =  fsolve(@fun_pullback6  , var_in ,options)    ;

 
 

x1(1:3*N1,1) = x(1:3*N1,1);

    str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '_step_' num2str(p1)   '.txt'];
        
      % str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '_step_' '_ht_' num2str(ht) '_nstep_' num2str(nstep)   '.txt'];

        
        fileID = fopen([path str0],'w');
                  fprintf(fileID,'%30.16E   \r\n',x1 );
                 fclose(fileID);
 
 
end


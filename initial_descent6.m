function x =   initial_descent6()


global   N1  N     id   h   lm    p1 tau c ht N2 k bt

global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p bx0 by0 bz0 v0 branch x0 l path 
format longE
 

N1 = N-1;
h = 1/N;

 
str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*16)) '_step_' num2str(p1-1)   '.txt'] ;

%== previous steps ==
x0 =  load([path str0]);
 
bx0(2:N,1)     = x0(1:3:3*N1-2,1)    ;
by0(2:N,1)     = x0(2:3:3*N1-1,1)    ;
bz0(2:N,1)     = x0(3:3:3*N1,1)      ;

%==========================================================================

% Generating initial guess using negative eigen value direction ======

bx(1,1) = 1;
by(1,1) = 0;
bz(1,1) = 0;

bx(N+1,1) = -1;
by(N+1,1) =  0;
bz(N+1,1) =  0;
%-- hessian
id = [1 0 0;0 1 0;0 0 1];

%== Gradient of energy without the constraint term 

[fval0,jac] = fun_jacobian6([x0(1:3*N1,1)' zeros(1,2*N1+4)]');

gradG = jac(3*N1+1:5*N1+4,1:3*N1);
jac1  = jac(1:3*N1,1:3*N1);
 
%== lagrange multipliers corresponding to bx0, by0, and bz0 === 

lag = pinv(gradG')*fval0(1:3*N1,1); %(gradG*gradG')\(gradG*fval0(1:3*N1,1));  % lagrange multiplier corresponding to x0

% 
% 
[fval,jac] = fun_jacobian6([x0(1:3*N1,1)' lag']');   % Hessian evaluated at previous step
% 
  num_eq = 3*N1   ;    % number of equilibrium equations
  num_C  = 2*N1+4   ;   % number of constraints
% %
  A = (jac(num_eq+1:end,1:num_eq))';   % num_C x num_eq constraint gradient matrix
% 
% %--- Qr decomposition of the constraint matrix to obtain null space
% 
  [Q R] = qr(A)                ;
  Z     = Q(:,num_C+1:num_eq)  ;
% 
% %   %---- Projected hessian ----
% %
  Ac = Z'*jac(1:num_eq,1:num_eq)*Z     ;
%
%   %--- eigsen values and eigsen vectors ----
%

[V,D] = eigs(Ac,6,'smallestabs') ;
%
[B, I ] = sort(diag(D),'ascend');
 
[V1,D] = eigs(jac,6,'smallestabs') ;

[B, I ] = sort(diag(D),'ascend');
 

x(1:3*N1,1) = x0(1:3*N1,1)  + .1*V1(1:3*N1,I(1))/sqrt(sum(V1(1:3*N1,I(1))));



% %=== Evaluating hessian at the initial guess x ==
% [fval0,jac] = fun_jacobian6([x(1:3*N1,1)' zeros(1,2*N1+4)]');  % used only for function evaluation 
% 
% gradG = jac(3*N1+1:5*N1+4,1:3*N1);  % gradient of the constraint matrix at the initial guess 
% jac1  = jac(1:3*N1,1:3*N1);
% 
%  
% lag = (gradG*gradG')\(gradG*fval0(1:3*N1,1));  % lagrange multiplier corresponding to x 


%[fval,jac] = fun_jacobian6([x(1:3*N1,1)' lag']');   % Hessian evaluated at the initial guess 

%[V1,D1] = eig(jac1);

for i = 1:k 
    
    indv             = (3*N1)*(i-1)+1:(3*N1)*i           ;
    x(3*N1+indv,1) = V1(1:3*N1,I(i))./sqrt(sum(V1(1:3*N1,I(i)).^2));   
             v0(:,i) = x0(3*N1+indv,1)                 ; 
    
end

%== exponetial decay in the parameter l... E(l) = l^2/2

l0=.01;

l =1;% l0*exp(-p1*ht);



%=== eliminating the imaginary part =====

x = real(x);
x0 = real(x0);
v0 = real(v0);

bx0(1,1) = 1;
by0(1,1) = 0;
bz0(1,1) = 0;

bx0(N+1,1) = -1;
by0(N+1,1) =  0;
bz0(N+1,1) =  0;



%------------ Array initialization ----------------

bxp = zeros(N+1,1);
byp = zeros(N+1,1);
bzp = zeros(N+1,1);

bx2p = zeros(N+1,1);
by2p = zeros(N+1,1);
bz2p = zeros(N+1,1);

bx4p = zeros(N-1,1);
by4p = zeros(N-1,1);
bz4p = zeros(N-1,1);

 end
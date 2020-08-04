function x =   initial_descent5()


global   N1  N     id   h   lm    p1 tau c ht N2 k

global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p bx0 by0 bz0 v0 branch x0
format longE

 

N1 = N-1;
h = 1/N;

path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_descent7/'    ;

%path = '/Users/rtodres/Documents/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_descent5/'   ;

str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '_step_' num2str(p1-1)   '.txt'] ;

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

[fval,jac] = fun_jacobian6(x0(1:5*N1+4,1));   % Hessian evaluated at previous step

num_eq = 3*N1   ;    % number of equilibrium equations
num_C  = 2*N1+4   ;   % number of constraints
%
A = (jac(num_eq+1:end,1:num_eq))';   % num_C x num_eq constraint gradient matrix

%--- Qr decomposition of the constraint matrix to obtain null space

[Q R] = qr(A)                ;
Z     = Q(:,num_C+1:num_eq)  ;

%   %---- Projected hessian ----
%
Ac = Z'*jac(1:num_eq,1:num_eq)*Z     ;
%
%   %--- eigsen values and eigsen vectors ----
%
[V,D] = eigs(Ac,5,'smallestabs') ;
%
[B, I ] = sort(diag(real(D)),'ascend');
 
[V1,D1] = eig(jac);

x(1:5*N1+4,1) = x0(1:5*N1+4,1)  - ht*V1(I(1));


[fval,jac] = fun_jacobian6(x(1:5*N1+4,1));   % Hessian evaluated at the initial guess 

[V1,D1] = eig(jac);

for i = 1:k 
    
    indv             = (5*N1+4)*(i-1)+1:(5*N1+4)*i           ;
    x(5*N1+4+indv,1) = V1(:,I(i))./sqrt(sum(V1(:,I(i)).^2));   
             v0(:,i) = x0(5*N1+4+indv,1)                 ; 
    
end
    

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
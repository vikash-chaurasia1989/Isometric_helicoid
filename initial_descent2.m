function x =   initial_descent2()


global   N1  N     id   h   lm    p1 tau c ht N2

global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p bx0 by0 bz0  
format longE

N1 = N-1;
h = 1/N;

path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_descent2/'    ;
% %if(p1==2)
%     
%     %== previous step ==
%     
%     str0 = 'branch_213_N105_tau_160000000000_h_200_step_1.txt';
% %else
%     str0 = ['branch_' num2str(213) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '_h_' num2str(round(10^5*ht)) '_step_' num2str(p1-1) '.txt'];
% %end
str0 = ['branch_' num2str(21) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '_step_' num2str(N2) '_' num2str(p1-1)   '.txt'];


%== previous steps ==
x0 = load([path str0]);

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

[fval,jac] = fun_jacobian6(x0);   % Hessian evaluated at previous step

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

B

%-- perturbed shape with respect to negative eigen value ---

i11 = 1;

x1 = real(V(:,I(i11)));

[V1,D1] = eigs(jac);


% x(1:3*N1,1)  = x0(1:3*N1,1) + ht*Z*x1;
% x(3*N1+1:5*N1+4,1) = x0(3*N1+1:5*N1+4,1);

x = x0  - ht*V1(I(1));

% 
if(p1>2) 
    str0 = ['branch_' num2str(21) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '_step_' num2str(N2) '_' num2str(p1-2)   '.txt'];
   xm1 = load([path str0]);
   x = 2*x0-xm1;
end

%x(5*N1+4,1)=0;
     
    

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
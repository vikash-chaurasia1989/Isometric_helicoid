%function [] = pre_descent4()

global N fac tau branch bx by bz id k


N1 = N-1;


strpath = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch';
%strpath = '/Users/rtodres/Documents/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch';
path = [strpath num2str(branch) '/'];

strtau = ['tau_branch_' num2str(branch) '.txt'];

tau1 = load([path strtau]);

[val,ind] = min(abs(tau1-tau));

tau2 = tau1(ind);


str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau2)) '.txt'];



f1 = load([path str0]);
f1(3*N1+1:5*N1+4,1) = fac/0.001*f1(3*N1+1:5*N1+4,1);  % rescaled the lagrange multipliers for given sigma


 
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

[fval,jac] = fun_jacobian6(f1);   % Hessian evaluated at previous step

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
 
k= 0;
for i = 1:length(B)
    if(B(i) <0)
        k = k+1;
    end
end

for i = 1:k 
    
    indv = (N1-4)*(i-1)+1:(N1-4)*i;
    
    f1(5*N1+4+indv,1) = V(:,I(i))./sqrt(sum(V(:,I(i)).^2));   
end

 %=== saving it into data_descent folder ====
path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_descent4/'    ;
str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '_step_0'   '.txt'];

fileID = fopen([path str0],'w');
fprintf(fileID,'%30.16E   \r\n',f1 );
fclose(fileID);

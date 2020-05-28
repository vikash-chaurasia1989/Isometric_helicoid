function [] = pre_descent2(N2)

global N fac tau d21 d13 d23 


N1 = N-1;


%======================== Branch 1 ========================================
branch = 1;
path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];

strtau = ['tau_branch_' num2str(branch) '.txt'];

tau1 = load([path strtau]);

[val,ind] = min(abs(tau1-tau));

tau2 = tau1(ind);


str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau2)) '.txt'];



f1 = load([path str0]);
f1(3*N1+1:5*N1+4,1) = fac/0.001*f1(3*N1+1:5*N1+4,1);  % rescaled the lagrange multipliers for given sigma

%=== saving it into data_descent folder ====
path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_descent2/'    ;
str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '_step_' num2str(N2) '_1'   '.txt'];

fileID = fopen([path str0],'w');
fprintf(fileID,'%30.16E   \r\n',f1 );
fclose(fileID);


%==========================================================================
%==========================================================================


%======================== Branch 1 ========================================
branch = 2;
path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];

strtau = ['tau_branch_' num2str(branch) '.txt'];

tau1 = load([path strtau]);

[val,ind] = min(abs(tau1-tau));

tau2 = tau1(ind);


%str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau2)) '.txt'];

str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau2)) '.txt'];

f2 = load([path str0]);
f2(3*N1+1:5*N1+4,1) = fac/0.001*f2(3*N1+1:5*N1+4,1);  % rescaled the lagrange multipliers for given

%=== saving it into data_descent folder ====
path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_descent2/'    ;
str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '_step_' num2str(N2) '_1'   '.txt'];

str0 = ['branch_' num2str(21) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '_step_' num2str(N2) '_1'   '.txt'];


fileID = fopen([path str0],'w');
fprintf(fileID,'%30.16E   \r\n',f2 );
fclose(fileID);


%==========================================================================
%==========================================================================

%======================== Branch 1 ========================================
branch =3;
path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];

strtau = ['tau_branch_' num2str(branch) '.txt'];

tau1 = load([path strtau]);

[val,ind] = min(abs(tau1-tau));

tau2 = tau1(ind);


str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau2)) '.txt'];


f3 = load([path str0]);
f3(3*N1+1:5*N1+4,1) = fac/0.001*f3(3*N1+1:5*N1+4,1);  % rescaled the lagrange multipliers for given

%=== saving it into data_descent folder ====
path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_descent2/'    ;
str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '_step_' num2str(N2) '_1'   '.txt'];

fileID = fopen([path str0],'w');
fprintf(fileID,'%30.16E   \r\n',f3 );
fclose(fileID);


%==========================================================================
%==========================================================================

%==== Calculating Frechet distance between the three curves

%=== Loading binormal curve ====

bx1(2:N,1)     = f1(1:3:3*N1-2,1)    ;
by1(2:N,1)     = f1(2:3:3*N1-1,1)    ;
bz1(2:N,1)     = f1(3:3:3*N1,1)      ;

bx1(1,1) = 1;
by1(1,1) = 0;
bz1(1,1) = 0;

bx1(N+1,1) = -1;
by1(N+1,1) =  0;
bz1(N+1,1) =  0;


%=====================================
bx2(2:N,1)     = f2(1:3:3*N1-2,1)    ;
by2(2:N,1)     = f2(2:3:3*N1-1,1)    ;
bz2(2:N,1)     = f2(3:3:3*N1,1)      ;

bx2(1,1) = 1;
by2(1,1) = 0;
bz2(1,1) = 0;

bx2(N+1,1) = -1;
by2(N+1,1) =  0;
bz2(N+1,1) =  0;

%====================================

bx3(2:N,1)     = f3(1:3:3*N1-2,1)    ;
by3(2:N,1)     = f3(2:3:3*N1-1,1)    ;
bz3(2:N,1)     = f3(3:3:3*N1,1)      ;

bx3(1,1) = 1;
by3(1,1) = 0;
bz3(1,1) = 0;

bx3(N+1,1) = -1;
by3(N+1,1) =  0;
bz3(N+1,1) =  0;

%==================================== 

d21 = sum(sqrt((bx2-bx1).^2 + (by2-by1).^2 + (bz2-bz1).^2)) ;

d13 = sum(sqrt((bx3-bx1).^2 + (by3-by1).^2 + (bz3-bz1).^2)) ;

d23 = sum(sqrt((bx2-bx3).^2 + (by2-by3).^2 + (bz2-bz3).^2)) ;

end

function f = pre_descent()

global N fac tau


N1 = N-1;


%======================== Branch 1 ========================================
branch = 1;
path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];

strtau = ['tau_branch_' num2str(branch) '.txt'];

tau1 = load([path strtau]);

[val,ind] = min(abs(tau1-tau));

tau2 = tau1(ind);


str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau2)) '.txt'];


f = load([path str0]);
f(3*N1+1:5*N1+4,1) = fac/0.001*f(3*N1+1:5*N1+4,1);  % rescaled the lagrange multipliers for given

%=== saving it into data_descent folder ====
path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_descent/'    ;
str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) 'step0'   '.txt'];

fileID = fopen([path str0],'w');
fprintf(fileID,'%30.16E   \r\n',f );
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


str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau2)) '.txt'];


f = load([path str0]);
f(3*N1+1:5*N1+4,1) = fac/0.001*f(3*N1+1:5*N1+4,1);  % rescaled the lagrange multipliers for given

%=== saving it into data_descent folder ====
path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_descent/'    ;
str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) 'step0'   '.txt'];

fileID = fopen([path str0],'w');
fprintf(fileID,'%30.16E   \r\n',f );
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


f = load([path str0]);
f(3*N1+1:5*N1+4,1) = fac/0.001*f(3*N1+1:5*N1+4,1);  % rescaled the lagrange multipliers for given

%=== saving it into data_descent folder ====
path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_descent/'    ;
str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) 'step0'   '.txt'];

fileID = fopen([path str0],'w');
fprintf(fileID,'%30.16E   \r\n',f );
fclose(fileID);


%==========================================================================
%==========================================================================
 


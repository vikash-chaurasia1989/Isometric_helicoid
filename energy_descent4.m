clear all
clc

global N N1 h tau sig fac branch ht bt gamma t1 

parameters5(); 


%== prefactor for aspect ratio
 

branch = 2;

path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_descent4/'    ;
% path = '/Users/rtodres/Documents/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_descent5/';

        
for p1= 1:971%length(t1)

    str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*16)) '_step_' num2str(p1-1)   '.txt'];


    x = load([path str0]);

    bx(2:N,1)     = x(1:3:3*N1-2,1)    ;
    by(2:N,1)     = x(2:3:3*N1-1,1)    ;
    bz(2:N,1)     = x(3:3:3*N1,1)      ;

    bx(1,1) = 1;
    by(1,1) = 0;
    bz(1,1) = 0;

    bx(N+1,1) = -1;
    by(N+1,1) =  0;
    bz(N+1,1) =  0;

    E(p1) = energy_b(bx,by,bz);
end

branch = 2;
strpath = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch';
%strpath = '/Users/rtodres/Documents/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch';
path = [strpath num2str(branch) '/'];

strtau = ['tau_branch_' num2str(branch) '.txt'];

tau1 = load([path strtau]);

[val,ind] = min(abs(tau1-tau));

tau2 = tau1(ind);


str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau2)) '.txt'];


x = load([path str0]);

bx(2:N,1)     = x(1:3:3*N1-2,1)    ;
by(2:N,1)     = x(2:3:3*N1-1,1)    ;
bz(2:N,1)     = x(3:3:3*N1,1)      ;

bx(1,1) = 1;
by(1,1) = 0;
bz(1,1) = 0;

bx(N+1,1) = -1;
by(N+1,1) =  0;
bz(N+1,1) =  0;

%E(p1+1) = energy_b(bx,by,bz);



plot(E,'-o')
hold on
plot(E(1),'-or')
hold on
%plot(p1+1,E(end),'-or')

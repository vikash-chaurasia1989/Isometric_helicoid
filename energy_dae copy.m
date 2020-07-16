clear all
clc

global N tau sig N1 ht nnt tspan tau

parameters();


 


%== prefactor for aspect ratio
sig = 0.01;

branch = 2;

 path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_dae/'    ;
    str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '_step_' '_ht_' num2str(ht) '_tspan_' num2str(length(tspan))   '.txt'];
    
    temp = load([path str0]);
    
    sz = size(temp);
for p1= 1:sz(2)
    
     
    x = temp(:,p1);

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

 
%E(p1+1) = energy_b(bx,by,bz);



plot(E,'-o')
hold on
plot(E(1),'-or')
hold on
%plot(p1+1,E(end),'-or')

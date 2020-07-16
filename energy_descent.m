clear all
clc

global N tau sig N1


N   = 105 ;
N1=N-1;
tau = 16  ;
sig = 0.01;

%ht = 0.005;
%c =  0.2;

ht = 0.04;
 
t1 = 0:ht:1;

N2 = 1000;



for p1= 1:length(t1)
    
   % str0 = ['branch_' num2str(213) '_N' num2str(N) '_tau_' num2str(round(10^10*16)) 'step' num2str(p1) '.txt'];
    
%    
%    if(p1==1)
%        path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_descent/';
%        
%        str0 = ['branch_' num2str(2) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) 'step0'   '.txt'];
%    elseif(p1<52)
%        path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_descent/'    ;
%        
%        str0 = ['branch_' num2str(213) '_N' num2str(N) '_tau_' num2str(round(10^10*16)) 'step' num2str(p1-1) '_c_' num2str(round(1000*c)) '_ht_' num2str(round(1000*ht)) '.txt'];
%        
%    else
%        path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_descent/';
%        
%        str0 = ['branch_' num2str(1) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) 'step0'   '.txt'];
%    end
%     
path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_descent2/'    ;

if(p1==1)
    
    %== previous step ==
    
    str0 = 'branch_213_N105_tau_160000000000_h_200_step_1.txt';

else
     str0 = ['branch_' num2str(213) '_N' num2str(N) '_tau_' num2str(round(10^10*round(tau))) '_h_' num2str(round(10^5*ht)) '_step_' num2str(p1) '.txt'];
end 

    str0 = ['branch_' num2str(21) '_N' num2str(N) '_tau_' num2str(round(10^10*16)) '_step_' num2str(N2) '_' num2str(p1)   '.txt'];

    
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

path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_descent/';

path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch2/';


str0 = ['branch_' num2str(2) '_N' num2str(N) '_tau_' num2str(round(10^10*round(tau))) 'step0'   '.txt'];

str0 = 'branch_2_N105_tau_160000000000.txt';

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

E(p1+1) = energy_b(bx,by,bz);



plot(E,'-o')
hold on
plot(E(1),'-or')
hold on
plot(p1+1,E(end),'-or')
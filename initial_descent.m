function x =   initial_descent()


global   N1  N     id   h   lm    p1 tau c ht

global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p bx0 by0 bz0 steps
format longE

N1 = N-1;
h = 1/N;


path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_descent/';
 
p1==1;
if(p1==1)
    
    str0 = ['branch_' num2str(2) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) 'step0'   '.txt'];
    
    path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_evolution2/';

    str0 =   'b_branch_213_N_105_t_5_tau_160000000000.txt';%'b_branch_213_N_105_t_10_tau_160000000000.txt';
    %== previous steps ==
    x0 = load([path str0]);
%     
%     bx0(2:N,1)     = x0(1:3:3*N1-2,1)    ;
%     by0(2:N,1)     = x0(2:3:3*N1-1,1)    ;
%     bz0(2:N,1)     = x0(3:3:3*N1,1)      ;
 
   bx0 = x0(:,1);
   by0 = x0(:,2);
   bz0 = x0(:,3);
 
else
   % str0 = ['branch_' num2str(213) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) 'step' num2str(p1-1) '.txt'];
    str0 = ['branch_' num2str(213) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) 'step' num2str(p1-1) '_c_' num2str(round(1000*c)) '_ht_' num2str(round(1000*ht)) '.txt'];

     %== previous steps ==
    x0 = load([path str0]);
    bx0(2:N,1)     = x0(1:3:3*N1-2,1)    ;
    by0(2:N,1)     = x0(2:3:3*N1-1,1)    ;
    bz0(2:N,1)     = x0(3:3:3*N1,1)      ;
end

path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_descent/';

      str0 = ['branch_' num2str(2) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) 'step0'   '.txt'];
    
       x0 = load([path str0]);
 
 
    %=== now generating initial guess using isotopy between two branches ===

    str1 = ['branch_' num2str(1) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) 'step0'   '.txt'];
 
    x1 = load([path str1]);
    
    if(p1<length(steps))
        t = (steps(p1+1)-steps(1))/(steps(end)-steps(1));
    else
        t = 1;
    end
    
    x = (1-t)*x0 + t*x1;
    
 
    if(p1>2)
     str0 = ['branch_' num2str(213) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) 'step' num2str(p1-2) '_c_' num2str(round(1000*c)) '_ht_' num2str(round(1000*ht)) '.txt'];
      xm1 = load([path str0]);
      
    x = 2*x0-xm1;
    end
%   str0 = ['branch_' num2str(213) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) 'step' num2str(p1) '_c_' num2str(round(1000*c)) '_ht_' num2str(round(1000*ht)) '.txt'];
%       x  = load([path str0]);
    
%   str0 = ['branch_' num2str(2) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) 'step0'   '.txt'];
%   x0 = load([path str0]);
%   str0 = ['branch_' num2str(213) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) 'step' num2str(p1+1) '_c_' num2str(round(1000*c)) '_ht_' num2str(round(1000*ht)) '.txt'];
%   x1 = load([path str0]); 
%     
%   x = .5*(x0+x1);
      
bx(1,1) = 1;
by(1,1) = 0;
bz(1,1) = 0;

bx(N+1,1) = -1;
by(N+1,1) =  0;
bz(N+1,1) =  0;


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



%-- hessian
id = [1 0 0;0 1 0;0 0 1];


 
end
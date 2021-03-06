function f =   initial_evolution()


global   N1  N    id   h    lm  tau   p1 fac p q

global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p nfold b0x b0y b0z b2x b2y b2z sig 
format longE



path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_evolution/';

 
nfold = 3;
 
%N = 72;
if(p1==1)
    tau2 =  tau ;
    
    if(nfold==3)
        str0 = ['3fold_N' num2str(N) '_tau_' num2str(10^10*tau2) '.txt'];
    elseif(nfold==5)
        str0 =['5pi_knot_N' num2str(N) '_tau_' num2str(10^10*tau2) '.txt'];
        %
    else
        str0 = ['7pi_knot_N' num2str(N) '_tau_' num2str(10^10*tau2) '.txt'];
    end
    
    
    strb =  [path 'b_' str0] ;
    
    
    %==== load b ===========
    temp = load(strb);
    N   = length(temp)-1;
    N1 = N-1;
    
    h = 1/N;
  
    
    %== initial guess for b
    bx = temp(1:N+1,1);
    by = temp(1:N+1,2);
    bz = temp(1:N+1,3);
    
    %=== base curve from which distance is to be maintained ===
    
    b0x = bx;
    b0y = by;
    b0z = bz;   
elseif(p1==2)
    
     tau2 =  tau ;
    
    if(nfold==3)
        str0 = ['3fold_N' num2str(N) '_tau_' num2str(10^10*tau2) '.txt'];
    elseif(nfold==5)
        str0 =['5pi_knot_N' num2str(N) '_tau_' num2str(10^10*tau2) '.txt'];
        %
    else
        str0 = ['7pi_knot_N' num2str(N) '_tau_' num2str(10^10*tau2) '.txt'];
    end
    
    
     
    %==== load b ===========
    temp = load([path 'b_' str0]);
    N   = length(temp)-1;
    N1 = N-1;
    
    h = 1/N;
  
    
    %== initial guess for b
    bxm1 = temp(1:N+1,1);
    bym1 = temp(1:N+1,2);
    bzm1 = temp(1:N+1,3);
    
    %=== base curve from which distance is to be maintained ===
    
    str1 = ['branch_' num2str(p) num2str(q) '_N_' num2str(N) '_step_' num2str(p1-1) '_tau_' num2str(10^10*tau) '.txt'];
    temp = load([path 'b_' str1]);
    
      %== initial guess for b
    b0x = temp(1:N+1,1);
    b0y = temp(1:N+1,2);
    b0z = temp(1:N+1,3);
    
    
    bx = 2*b0x-bxm1;
    by = 2*b0y-bym1;
    bz = 2*b0z-bzm1;
    
    

    %== initial guess is linear extrapolation from last two steps 
    
else 
    str0 = ['branch_' num2str(p) num2str(q) '_N_' num2str(N) '_step_' num2str(p1-2) '_tau_' num2str(10^10*tau) '.txt'];
     
    %==== load b ===========
    temp = load([path 'b_' str0]);
    N   = length(temp)-1;
    N1 = N-1;
    
    h = 1/N;
  
    
    %== initial guess for b
    bxm1 = temp(1:N+1,1);
    bym1 = temp(1:N+1,2);
    bzm1 = temp(1:N+1,3);
    
    %=== base curve from which distance is to be maintained ===
    
    str1 = ['branch_' num2str(p) num2str(q) '_N_' num2str(N) '_step_' num2str(p1-1) '_tau_' num2str(10^10*tau) '.txt'];
    temp = load([path 'b_' str1]);
    
      %== initial guess for b
    b0x = temp(1:N+1,1);
    b0y = temp(1:N+1,2);
    b0z = temp(1:N+1,3);
    
    
    bx = 2*b0x-bxm1;
    by = 2*b0y-bym1;
    bz = 2*b0z-bzm1;
    
    
    
    
    
end

  str0 =['5pi_knot_N' num2str(N) '_tau_' num2str(10^10*tau) '.txt'];
 %==== load b ===========
    temp = load([path 'b_' str0]);
   
    
    %== initial guess for b
    b2x = temp(1:N+1,1);
    b2y = temp(1:N+1,2);
    b2z = temp(1:N+1,3);
        
%   
%     bx = 2*b0x-bxm1;
%     by = 2*b0y-bym1;
%     bz = 2*b0z-bzm1;
%           
norm = sqrt(bx.^2+by.^2+bz.^2);

bx = bx./norm;
by = by./norm;
bz = bz./norm;

strlm = [path 'lm_' str0];
strrho = [path 'rho_' str0];
struvw = [path 'uvw_' str0];



b(1:3:3*N1-2,1) = bx(2:N,1)  ;
b(2:3:3*N1-1,1) = by(2:N,1)  ;
b(3:3:3*N1,1)   = bz(2:N,1)  ;
%====================


 

temp =     load(struvw);
%
u  =    temp(1,1);
v  =    temp(2,1);
w =     temp(3,1);

%===============================

%=== load rho ==========
rho =      load(strrho);


% %=== load lm ====
lm =    load(strlm);
del = 0;
%f =  [b(1:3*N1,1)'  rho' lm' u v    w  del ]'  ;         % for fun_jacobian6
f = b(1:3*N1,1);

%--- boundary points --

bx(1,1) = 1;
by(1,1) = 0;
bz(1,1) = 0;

bx(N+1,1) = -1;
by(N+1,1) =  0;
bz(N+1,1) =  0;


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
function f =   initial_guess_evolution()


global   N1  N  qd   id   h   rho lm  tau tau1 p1 fac

global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p nfold sig branch
format longE

p1 = 101;
%N = 72;
if(p1>0)
    tau2 =  tau1(p1);
    
    if(nfold==3)
        str0 = ['3fold_N' num2str(N) '_tau_' num2str(10^10*tau2) '.txt'];
        %  str0 = ['3fold_N' num2str(N) '_tau_' num2str(10^10*tau2) '_sig_' num2str(sig) '.txt'];
        
        path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_3pi_3/';
    elseif(nfold==5)
        str0 =['5pi_knot_N' num2str(N) '_tau_' num2str(10^10*tau2) '.txt'];
        path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_5pi_3/';
        %
    else
        str0 = ['7pi_knot_N' num2str(N) '_tau_' num2str(10^10*tau2) '.txt'];
        path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_7pi_3/';
    end
    
    
end


strb =  [path 'b_' str0] ;
strlm = [path 'lm_' str0];
strrho = [path 'rho_' str0];
struvw = [path 'uvw_' str0];


%==== load b ===========
temp = load(strb);
N   = length(temp)-1;
N1 = N-1;

h = 1/N;
% N = N/2;
bx1 = temp(1:N+1,1);
by1 = temp(1:N+1,2);
bz1 = temp(1:N+1,3);





temp =     load(struvw);
%
u  =    temp(1,1);
v  =    temp(2,1);
w =     temp(3,1);

%===============================

%=== load rho ==========
rho1 =      load(strrho);


% %=== load lm ====
lm1 =    load(strlm);
 

%================  Interpolation of the initial guess ==============

%======== b =========
% %
N2 = 105;        % new number of points -1


 s1 = (linspace(0,1,N+1))';
 s2 = (linspace(0,1,N2+1))';
 
 bx = interp1q(s1,bx1,s2)  ;
 by = interp1q(s1,by1,s2)  ;
 bz = interp1q(s1,bz1,s2)  ;

 s1 = (linspace(0,1,N))';
 s2 = (linspace(0,1,N2))';
 
 
lm  = (interp1q(s1,lm1,s2));
 
 s1 = (linspace(0,1,N-1))';
 s2 = (linspace(0,1,N2-1))';
 
 
 rho = (interp1q(s1,rho1,s2));
 
 

N = N2;
N1 = N-1;

h = 1/N;




b(1:3:3*N1-2) = bx(2:N,1)  ;
b(2:3:3*N1-1) = by(2:N,1)  ;
b(3:3:3*N1)   = bz(2:N,1)  ;
%====================



f =  [b(1:3*N1)  rho' lm' u v    1*w  ]'  ;         % for fun_jacobian6



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


figure(1)

plotbrowser on
title('Initial guess')
% plot3(bx,by,bz,'LineWidth',4)
plot(lm,'-o')
set(gca,'DefaultAxesFontSize',4,'FontSize',25,'LineWidth',2)
box on
grid on

end
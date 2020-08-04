clear all
clc

global    lm  N1 p1 count tau  rho f_b  h tx ty tz  u v w tau1 p1 fac

global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p  N   len rx ry rz Lk err   nfold branch sig w1 w2 w3 tau2

format longE
tic


col1 = [0,     0.976, 0.968]  ;      % -- color for mode n = 2
col2 = [0.831, 0.545, 0.247]  ;      % -- color for mode n = 3
col3 = [0.341, 0.505, 0.819]  ;      % -- color for mode n = 4
col4 = [0.705, 0.701, 0.070]  ;      % -- color for mode n = 5
col5 = [0.301, 0.811, 0.498]  ;      % -- color for mode n = 6
col6 = [1, 0.929, 0.278];
col = {col1,col2,col3,col4,col5,col6} ;
colw = [41,44,52]/255;




branch = 3;

path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];
%path = ['/Users/rtodres/Documents/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];


strtau = ['tau_branch_' num2str(branch) '.txt'];

tau1 = load([path strtau]);

p2 = 22 ;

N = 105;
N1 = N-1;
h = 1/N;

%=== most symmetric 5pi knot ==
%str0 = 'branch_2_N120_tau_119842500480_symmetry.txt';%;  % tau = 11.9842500480

for p1  = p2:p2%2:length(tau1)%53:length(tau1)%53:length(tau1)%308:332%1:1%length(tau1)%length(tau1)%length(tau1)%%262:length(tau1)%1:length(tau1)%3:length(tau1)%48:length(tau1)
    
    
    tau =  16% tau1(p1);%8.093946633549898e+00
    
    str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '.txt'];
    
    x = load([path str0]);
    
    %----- reading input and constructing bx by bz ----
    bx(1,1) = 1;
    by(1,1) = 0;
    bz(1,1) = 0;
    
    bx(N+1,1) = -1;
    by(N+1,1) =  0;
    bz(N+1,1) =  0;
    bx(2:N,1)     = x(1:3:3*N1-2,1)    ;
    by(2:N,1)     = x(2:3:3*N1-1,1)    ;
    bz(2:N,1)     = x(3:3:3*N1,1)      ;
    
    
    rho           = x(3*N1+1:4*N1,1)   ;
    lm            = x(4*N1+1:5*N1+1,1) ;
    
    u             = x(5*N1+2 ,1)       ;
    v             = x(5*N1+3,1)        ;
    w             = x(5*N1+4,1)        ;
    
    bext(1:3,1) = -1*[bx(N,1);by(N,1);bz(N,1)];
    bext(4:6,1) = -1*[bx(2,1);by(2,1);bz(2,1)];
    
    %---- derivative calculation for b2 --- bN  ---
    
    ig = 2:N;
    
    
    
    
    %=== O(h)====
    bxp(1,1) = (bx(1,1)-bext(1,1))/(h);
    byp(1,1) = (by(1,1)-bext(2,1))/(h);
    bzp(1,1) = (bz(1,1)-bext(3,1))/(h);
    
    bxp(ig,1) = (bx(ig,1)-bx(ig-1,1))/(h);
    byp(ig,1) = (by(ig,1)-by(ig-1,1))/(h);
    bzp(ig,1) = (bz(ig,1)-bz(ig-1,1))/(h);
    
    
    bxp(N+1,1) = (bx(N+1,1)-bx(N,1))/(h);
    byp(N+1,1) = (by(N+1,1)-by(N,1))/(h);
    bzp(N+1,1) = (bz(N+1,1)-bz(N,1))/(h);
    
    
    
    %--  b''  (O(h^2))----
    
    bx2p(1,1) = (bx(2,1) + bext(1,1) -2*bx(1,1))/h^2 ;
    by2p(1,1) = (by(2,1) + bext(2,1) -2*by(1,1))/h^2 ;
    bz2p(1,1) = (bz(2,1) + bext(3,1) -2*bz(1,1))/h^2 ;
    
    bx2p(ig,1) = (bx(ig+1,1)  + bx(ig-1,1) - 2*bx(ig,1))/h^2 ;
    by2p(ig,1) = (by(ig+1,1)  + by(ig-1,1) - 2*by(ig,1))/h^2 ;
    bz2p(ig,1) = (bz(ig+1,1)  + bz(ig-1,1) - 2*bz(ig,1))/h^2 ;
    
    bx2p(N+1,1) = (bext(4,1) + bx(N,1)   - 2*bx(N+1,1))/h^2 ;
    by2p(N+1,1) = (bext(5,1) + by(N,1)   - 2*by(N+1,1))/h^2 ;
    bz2p(N+1,1) = (bext(6,1) + bz(N,1)   - 2*bz(N+1,1))/h^2 ;
    
    %
    i = 1:N+1   ;
    
    %
    tx(i,1) = by(i,1).*bzp(i,1) - bz(i).*byp(i,1);
    ty(i,1) = bz(i,1).*bxp(i,1) - bx(i).*bzp(i,1);
    tz(i,1) = bx(i,1).*byp(i,1) - by(i).*bxp(i,1);
    tx = tx/tau;
    ty = ty/tau;
    tz = tz/tau;
    
    %--------------------
    
    % %--- position vector using integration of tangent
    % % initialization
    %h=1;
    rx(1) = 0; ry(1) = 0; rz(1) = 0;
    
    
    
    for i = 1:N
        
        rx(i+1) = rx(i) + h*tx(i+1);
        ry(i+1) = ry(i) + h*ty(i+1);
        rz(i+1) = rz(i) + h*tz(i+1);
        
    end
    
    
    
    err1  = u*tx + v*ty + w*tz;
    
    kappa = sqrt((bx2p.*bx2p + by2p.*by2p + bz2p.*bz2p -tau^4)/tau^2);
    s0 = linspace(0,1,length(kappa));
    s1 = linspace(0,1,18*20);
    
    kappa2 = interp1q(s0',kappa,s1');
    

    
end

plot(kappa2)

% % golden ratio = 1.618033988749895e+00
% 
% 
% %=== comparing numerical and elliptic curvature for tau = 8.093946633549898e+00
% path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];
% strkp = 'curvature_numerical';% ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '.txt'];
% 
% kappa_num = load([path strkp]);
% s_num = linspace(0,1,length(kappa_num));
% 
% %=== Elliptic curvature ===
% A = 4.037960852507271e+01;
% B = 1.679359081062158e+04;
% al1 = 2*A + 2*sqrt(A^2+B);
% 
% al2 = 0;
% al3 = -2*A +2*sqrt(A^2+B);
% 
% p = sqrt((al3-al2)/(al3+al1));
% q = sqrt(1-al2/al3);
% r = 1/2*sqrt(al3+al1);
% 
% 
% [K,E] = ellipke(p);
% N = 1440;
% 
% s = linspace(0,2*nfold*K/r,N+1);
% tau =  8.093946628372475e+00;%tau1*2*nfold*K/r;
% 
% h = 1/N*2*nfold*K/r;
% 
% [sn,cn,dn] = ellipj(r*s,p);
% 
% 
% kappa =  sqrt(al3*(1-q^2*sn.^2));
% nfold = 3;
% for m = 1:2:nfold-1
%     
%     kappa((2*m-1)*N/(2*nfold)+1:((2*m-1)*N/(2*nfold) + N/(nfold))) = -kappa((2*m-1)*N/(2*nfold)+1:((2*m-1)*N/(2*nfold) + N/(nfold)));
%     
% end
% m = nfold;
% kappa((2*m-1)*N/(2*nfold)+1:N+1) = -kappa((2*m-1)*N/(2*nfold)+1:N+1);
% kappa = kappa';
% 
% 
% kappa_num = interp1q(s_num',kappa_num,s');
% 
% plot(s,(kappa-kappa_num)/max(kappa))
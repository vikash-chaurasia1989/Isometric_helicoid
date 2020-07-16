clear all
clc

global N tau sig N1 ht tau h tspan 

parameters();


%== prefactor for aspect ratio
sig = 0.01;

branch = 2;

path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_dae/'    ;
str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '_step_' '_ht_' num2str(ht) '_tspan_' num2str(length(tspan))   '.txt'];

temp = load([path str0]);

sz = size(temp);
bx(1,1) = 1;
by(1,1) = 0;
bz(1,1) = 0;

bx(N+1,1) = -1;
by(N+1,1) =  0;
bz(N+1,1) =  0;
for p1= 1:sz(2)
    
    x = temp(:,p1);
    
    %----- reading input and constructing bx by bz ----
    
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
    
    om =  [ 0  -w   v  ;
            w   0  -u  ;
           -v   u   0] ;   %-- \gamma_x tensor
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
    
    
    %   lmp = (lm(ig+1,1)-lm(ig-1))/(2*h);
    lmp = (lm(ig,1)-lm(ig-1,1))/(h);
    
    
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
    
    %---------------------------------------
    
    f_c = [0;0;0];
    
    
    p = 2;
    
    f_b(p-1,1) = 1/2*(bx(p,1)^2 + by(p,1)^2 + bz(p,1)^2 -1);
    
    f_b(N1+p-1,1) = 1/2*(bxp(p,1)^2 + byp(p,1)^2 + bzp(p,1)^2 - tau^2);
    
    %--- closure -----
    %
    f_c  = f_c       + [by(p-1,1)*bz(p,1) - bz(p-1,1)*by(p,1)    ;
        bz(p-1,1)*bx(p,1) - bx(p-1,1)*bz(p,1)    ;
        bx(p-1,1)*by(p,1) - by(p-1,1)*bx(p,1) ] ;
    %
    
    %----- Euler--Lagrange and Jacobian at the boundary points
    
    p = 3;
    
    
    f_b(p-1,1) = 1/2*(bx(p,1)^2 + by(p,1)^2 + bz(p,1)^2 -1);
    
    f_b(N1+p-1,1) = 1/2*(bxp(p,1)^2 + byp(p,1)^2 + bzp(p,1)^2 - tau^2);
    
    %--- closure -----
    
    f_c  = f_c       + [by(p-1,1)*bz(p,1) - bz(p-1,1)*by(p,1) ;
        bz(p-1,1)*bx(p,1) - bx(p-1,1)*bz(p,1) ;
        bx(p-1,1)*by(p,1) - by(p-1,1)*bx(p,1) ] ;
    
    
    %======================================================================
    
    
    
    for p = 4:N-2
        
        
        f_b(p-1,1) = 1/2*(bx(p,1)^2 + by(p,1)^2 + bz(p,1)^2 -1);
        
        f_b(N1+p-1,1) = 1/2*(bxp(p,1)^2 + byp(p,1)^2 + bzp(p,1)^2 - tau^2);
        
        %--- closure -----
        
        f_c  = f_c       + [by(p-1,1)*bz(p,1) - bz(p-1,1)*by(p,1) ;
            bz(p-1,1)*bx(p,1) - bx(p-1,1)*bz(p,1) ;
            bx(p-1,1)*by(p,1) - by(p-1,1)*bx(p,1) ] ;
        
        
        
        %======================================================================
    end
    
    %----- Euler--Lagrange and Jacobian at the boundary points
    
    p = N-1;
    
    f_b(p-1,1) = 1/2*(bx(p,1)^2 + by(p,1)^2 + bz(p,1)^2 -1);
    
    f_b(N1+p-1,1) = 1/2*(bxp(p,1)^2 + byp(p,1)^2 + bzp(p,1)^2 - tau^2);
    
    %--- closure -----
    
    f_c  = f_c       + [by(p-1,1)*bz(p,1) - bz(p-1,1)*by(p,1) ;
        bz(p-1,1)*bx(p,1) - bx(p-1,1)*bz(p,1) ;
        bx(p-1,1)*by(p,1) - by(p-1,1)*bx(p,1) ] ;
    
    
    %----- Euler--Lagrange and Jacobian at the boundary points
    
    p = N;
    
    
    f_b(p-1,1) = 1/2*(bx(p,1)^2 + by(p,1)^2 + bz(p,1)^2 -1)         ;
    
    f_b(N1+p-1,1) = 1/2*(bxp(p,1)^2 + byp(p,1)^2 + bzp(p,1)^2 - tau^2) ;
    
    %--- closure -----
    
    f_c  = f_c       + [  by(p-1,1)*bz(p,1) - bz(p-1,1)*by(p,1)   ;
        bz(p-1,1)*bx(p,1) - bx(p-1,1)*bz(p,1)   ;
        bx(p-1,1)*by(p,1) - by(p-1,1)*bx(p,1) ] ;
    
    
    f_b(2*N1+2:2*N1+4,1) =  0.*[0;-bz(N);by(N)] + f_c +     [by(N,1)*bz(N+1,1) - bz(N,1)*by(N+1,1)   ;
        bz(N,1)*bx(N+1,1) - bx(N,1)*bz(N+1,1)   ;
        bx(N,1)*by(N+1,1) - by(N,1)*bx(N+1,1) ] ;
    
    
    f_b(2*N1+1,1) = 1/2*(bxp(N+1,1)^2 + byp(N+1,1)^2 + bzp(N+1,1)^2 - tau^2)   ;
    
    
    
    err(p1) = sqrt(sum(f_b.^2)) ;
    
    figure(1)
    plotbrowser on 
    plot3(bx,by,bz)
    hold on 
    
    
    
end

figure(2)
plotbrowser on 
plot(err,'--o')

function [f_b,jac] =   fun_pullback4(x)

%========================================================================
%====== In this function file, Lagrangian is given by
 %=== This is continuous  version of pullback5
 % This is constrained optimization. Given an initial guess, this code
 % finds point on the constrained manifold. No minimization of energy or
 % distance. 
 % 
%
%========================================================================
%========================================================================


global  N  N1  lm  err1 l    count  u v w  bx by bz     A tau h
%     x = initial_guess()     ;
count = count+1                  ;

%--- variables -----
bx(2:N,1)     = x(1:3:3*N1-2,1)    ;
by(2:N,1)     = x(2:3:3*N1-1,1)    ;
bz(2:N,1)     = x(3:3:3*N1,1)      ;

bext(1:3,1) = -1*[bx(N,1);by(N,1);bz(N,1)];
bext(4:6,1) = -1*[bx(2,1);by(2,1);bz(2,1)];

ig = 2:N;


f_b = zeros(2*N1+4,1);
jac = zeros(2*N1+4,3*N1);

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


f_c = [0;0;0];



p = 2;

ind = 3*p-5:3*p-3;
%---- constraints ---------
%--  Unity constraint -----

f_b(p-1,1) = 1/2*(bx(p,1)^2+by(p,1)^2 + bz(p,1)^2-1);

%--- inextensibility -----
f_b(N1+p-1,1)=  1/2*(bxp(p,1)^2 + byp(p,1)^2 + bzp(p,1)^2 - tau^2); % sum((bi-bim1).*(bi-bim1)) -h^2;

%--- closure -----

f_c  = f_c       + [by(p-1,1)*bz(p,1) - bz(p-1,1)*by(p,1)   ;
                    bz(p-1,1)*bx(p,1) - bx(p-1,1)*bz(p,1)   ;
                    bx(p-1,1)*by(p,1) - by(p-1,1)*bx(p,1) ] ;
%======================================================================
%----- constraint gradient -----
%--- unit vector constraint --
jac(p-1,ind) = [bx(p) by(p) bz(p)];
%--- unit speed constraint --
jac(N1+p-1,ind)  =  1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)];

%---- closure ----
jac(2*N1+2:2*N1+4,ind) =   [   0                 -(bz(p-1)-bz(p+1))   (by(p-1)-by(p+1))  ;
                              (bz(p-1)-bz(p+1))       0             -(bx(p-1)-bx(p+1))  ;
                             -(by(p-1)-by(p+1))  (bx(p-1)-bx(p+1))           0        ] ;





for p = 3:N
    
    
    ind = 3*p-5:3*p-3;
    %---- constraints ---------
    %--  Unity constraint -----
    
    f_b(p-1,1) = 1/2*(bx(p,1)^2+by(p,1)^2 + bz(p,1)^2-1);
    
    %--- inextensibility -----
    f_b(N1+p-1,1)=  1/2*(bxp(p,1)^2 + byp(p,1)^2 + bzp(p,1)^2 - tau^2); % sum((bi-bim1).*(bi-bim1)) -h^2;
    
    %--- closure -----
    
    f_c  = f_c       + [by(p-1,1)*bz(p,1) - bz(p-1,1)*by(p,1) ;
                        bz(p-1,1)*bx(p,1) - bx(p-1,1)*bz(p,1) ;
                        bx(p-1,1)*by(p,1) - by(p-1,1)*bx(p,1) ] ;
    %======================================================================
    %----- constraint gradient -----
    %--- unit vector constraint --
    jac(p-1,ind) = [bx(p) by(p) bz(p)];
    %--- unit speed constraint --
    jac(N1+p-1,ind-3) = -1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)];
    
    jac(N1+p-1,ind)  =  1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)];
    
    %---- closure ----
    jac(2*N1+2:2*N1+4,ind) =   [ 0                 -(bz(p-1)-bz(p+1))   (by(p-1)-by(p+1))  ;
                                 (bz(p-1)-bz(p+1))       0             -(bx(p-1)-bx(p+1))  ;
                                -(by(p-1)-by(p+1))  (bx(p-1)-bx(p+1))           0        ] ;
    
    
 end

%--- end points ----
 
%=============== Function evaluation =========================

f_b(2*N1+1,1) = 1/2*(bxp(N+1,1)^2 + byp(N+1,1)^2 + bzp(N+1,1)^2 - tau^2)   ;


f_b(2*N1+2:2*N1+4,1) = f_c +     [by(N,1)*bz(N+1,1) - bz(N,1)*by(N+1,1)   ;
                                  bz(N,1)*bx(N+1,1) - bx(N,1)*bz(N+1,1)   ;
                                  bx(N,1)*by(N+1,1) - by(N,1)*bx(N+1,1) ] ;

jac(2*N1+1,ind)   = -1/h*[bxp(N+1,1) byp(N+1,1) bzp(N+1,1)];  %b'(N+1)


 


%
err1 = sqrt(sum(f_b.^2)) ;
end
function [F,jac] = fun_evolution(x)

%==== This is the version that gives same results as if discrete
%formulation is used.. We use this function to generate all the results for
%data_3pi_2

%=== Difference between fun_jacobian and fun_jacobian2 is that in fun_jacobian
%we used O(h^2) central difference scheme for x' while in fun_jacobian2, we
%used O(h) backward difference method

% My current understanding is that one should always use backward Euler for
% 1st order disretization.


%--- in this file, we solve for minimizing bending energy of the ruled surface --

%  b'''' + lm'*b' + lm*b'' - rho b + gm x b' = 0;
%
% --- constraints --
%   |b| = 1 and |b|' = \tau
%--- This is the final version of the analytical jacobian for finite
%element discretization of the continouous formulation

%=== Difference between fun_jacobian2 and this file is that in this file,
%we don't arrest the rotation of the binormal curve. Thus, in this file,
%we also solve for bz(2,1) and bz(N,1). We realized that arresting bz(2,1)
%and bz(N,1) to 0 in fun_jacobian2 led to overconstraining the system
%which resulted in incorrect stability results

%====  This is the version we will use for the paper
%===== Data generated using this code is saved in data_3pi_3, data_4pi_3, and data_5pi_3


global  N  N1  u v w      len id uni err qd bx2 E count bxN h tau rho lm bxg byg bzg sig E
global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p bext J2 fac b0x b0y b0z  d d1 b2x b2y b2z


% x = initial_guess_evolution();

count = count+1 ;
%----- reading input and constructing bx by bz ----

bx(2:N,1)     = x(1:3:3*N1-2,1)    ;
by(2:N,1)     = x(2:3:3*N1-1,1)    ;
bz(2:N,1)     = x(3:3:3*N1,1)      ;




bext(1:3,1) = -[bx(N,1);by(N,1);bz(N,1)];
bext(4:6,1) = -[bx(2,1);by(2,1);bz(2,1)];



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

%
%
%--  b''  (O(h^2))----

bx2p(1,1) = (bx(2,1) + bext(1,1) -2*bx(1,1))/h^2 ;
by2p(1,1) = (by(2,1) + bext(2,1) -2*by(1,1))/h^2 ;
bz2p(1,1) = (bz(2,1) + bext(3,1) -2*bz(1,1))/h^2 ;

bx2p(ig,1) = (bx(ig+1,1)  + bx(ig-1,1) - 2*bx(ig,1))/h^2 ;
by2p(ig,1) = (by(ig+1,1)  + by(ig-1,1)  - 2*by(ig,1))/h^2 ;
bz2p(ig,1) = (bz(ig+1,1)  + bz(ig-1,1)  - 2*bz(ig,1))/h^2 ;

bx2p(N+1,1) = (bext(4,1) + bx(N,1)   - 2*bx(N+1,1))/h^2 ;
by2p(N+1,1) = (bext(5,1) + by(N,1)   - 2*by(N+1,1))/h^2 ;
bz2p(N+1,1) = (bext(6,1) + bz(N,1)   - 2*bz(N+1,1))/h^2 ;
%

%    %=====================================================





bx4p(1,1) = (bx(4,1) - 4*bx(3,1) + 6*bx(2,1)-4*bx(1,1)+bext(1,1))/h^4;
by4p(1,1) = (by(4,1) - 4*by(3,1) + 6*by(2,1)-4*by(1,1)+bext(2,1))/h^4;
bz4p(1,1) = (bz(4,1) - 4*bz(3,1) + 6*bz(2,1)-4*bz(1,1)+bext(3,1))/h^4;

ig = 3:N-1;

bx4p(ig-1,1) = (bx(ig+2,1) - 4*bx(ig+1,1) + 6*bx(ig,1) - 4*bx(ig-1,1) + bx(ig-2,1))/(h^4);
by4p(ig-1,1) = (by(ig+2,1) - 4*by(ig+1,1) + 6*by(ig,1) - 4*by(ig-1,1) + by(ig-2,1))/(h^4);
bz4p(ig-1,1) = (bz(ig+2,1) - 4*bz(ig+1,1) + 6*bz(ig,1) - 4*bz(ig-1,1) + bz(ig-2,1))/(h^4);

bx4p(N-1,1)   = (bext(4,1)  - 4*bx(N+1,1)  + 6*bx(N,1)  - 4*bx(N-1,1) + bx(N-2,1) )/h^4;
by4p(N-1,1)   = (bext(5,1)  - 4*by(N+1,1)  + 6*by(N,1)  - 4*by(N-1,1) + by(N-2,1) )/h^4;
bz4p(N-1,1)   = (bext(6,1)  - 4*bz(N+1,1)  + 6*bz(N,1)  - 4*bz(N-1,1) + bz(N-2,1) )/h^4;


%---------------------------------------
 %===
n1 = tau/2/pi;
fac = (asinh(n1*pi*sig)/(4*pi^3*n1^3))  ;
f_c = [0;0;0];

 
p = 2;




f_b(p-1,1)    = 1/2*(bx(p,1)^2 + by(p,1)^2 + bz(p,1)^2 -1);
f_b(N1+p-1,1) = 1/2*(bxp(p,1)^2 + byp(p,1)^2 + bzp(p,1)^2 - tau^2);

%--- closure -----
%
f_c  = f_c       + [by(p-1,1)*bz(p,1) - bz(p-1,1)*by(p,1)    ;
                    bz(p-1,1)*bx(p,1) - bx(p-1,1)*bz(p,1)    ;
                    bx(p-1,1)*by(p,1) - by(p-1,1)*bx(p,1) ]  ;


%======================================================================
%====    jacobian ===
ind = 3*p-5:3*p-3;

jac(p-1,ind)           = [bx(p,1) by(p,1) bz(p,1)];
jac(N1+p-1,ind)        = 1/h*[bxp(p,1) byp(p,1) bzp(p,1)];
jac(2*N1+2:2*N1+4,ind) = [ 0                 -(bz(p-1)-bz(p+1))   (by(p-1)-by(p+1))  ;
                           (bz(p-1)-bz(p+1))       0             -(bx(p-1)-bx(p+1))  ;
                          -(by(p-1)-by(p+1))  (bx(p-1)-bx(p+1))           0        ] ; 
                        
jac(2*N1+5,ind)    = h*fac/h^2*([(bx2p(p-1,1)-2*bx2p(p,1)+ bx2p(p+1,1)) (by2p(p-1,1)-2*by2p(p,1)+ by2p(p+1,1)) (bz2p(p-1,1)-2*bz2p(p,1)+ bz2p(p+1,1))])       ;                  
                      

%----- Euler--Lagrange and Jacobian at the boundary points

p = 3;



f_b(p-1,1) = 1/2*(bx(p,1)^2 + by(p,1)^2 + bz(p,1)^2 -1);

f_b(N1+p-1,1) = 1/2*(bxp(p,1)^2 + byp(p,1)^2 + bzp(p,1)^2 - tau^2);

%--- closure -----

f_c  = f_c       + [by(p-1,1)*bz(p,1) - bz(p-1,1)*by(p,1) ;
    bz(p-1,1)*bx(p,1) - bx(p-1,1)*bz(p,1) ;
    bx(p-1,1)*by(p,1) - by(p-1,1)*bx(p,1) ] ;


%======================================================================
%====    jacobian ===
ind = 3*p-5:3*p-3;

jac(p-1,ind)           = [bx(p,1) by(p,1) bz(p,1)];

jac(N1+p-1,ind-3)      =-1/h*[bxp(p,1) byp(p,1) bzp(p,1)];
jac(N1+p-1,ind)        = 1/h*[bxp(p,1) byp(p,1) bzp(p,1)];

jac(2*N1+2:2*N1+4,ind) = [ 0                 -(bz(p-1)-bz(p+1))   (by(p-1)-by(p+1))  ;
                           (bz(p-1)-bz(p+1))       0             -(bx(p-1)-bx(p+1))  ;
                          -(by(p-1)-by(p+1))  (bx(p-1)-bx(p+1))           0        ] ; 
                        
jac(2*N1+5,ind)    = h*fac/h^2*([(bx2p(p-1,1)-2*bx2p(p,1)+ bx2p(p+1,1)) (by2p(p-1,1)-2*by2p(p,1)+ by2p(p+1,1)) (bz2p(p-1,1)-2*bz2p(p,1)+ bz2p(p+1,1))])       ;                  
                      

   

for p = 4:N-2
    
    
    f_b(p-1,1) = 1/2*(bx(p,1)^2 + by(p,1)^2 + bz(p,1)^2 -1);
    
    f_b(N1+p-1,1) = 1/2*(bxp(p,1)^2 + byp(p,1)^2 + bzp(p,1)^2 - tau^2);
    
    %--- closure -----
    
    f_c  = f_c       + [by(p-1,1)*bz(p,1) - bz(p-1,1)*by(p,1) ;
        bz(p-1,1)*bx(p,1) - bx(p-1,1)*bz(p,1) ;
        bx(p-1,1)*by(p,1) - by(p-1,1)*bx(p,1) ] ;
    
%======================================================================
%====    jacobian ===
ind = 3*p-5:3*p-3;

jac(p-1,ind)           = [bx(p,1) by(p,1) bz(p,1)];

jac(N1+p-1,ind-3)      =-1/h*[bxp(p,1) byp(p,1) bzp(p,1)];
jac(N1+p-1,ind)        = 1/h*[bxp(p,1) byp(p,1) bzp(p,1)];

jac(2*N1+2:2*N1+4,ind) = [ 0                 -(bz(p-1)-bz(p+1))   (by(p-1)-by(p+1))  ;
                           (bz(p-1)-bz(p+1))       0             -(bx(p-1)-bx(p+1))  ;
                          -(by(p-1)-by(p+1))  (bx(p-1)-bx(p+1))           0        ] ; 
                        
jac(2*N1+5,ind)    = h*fac/h^2*([(bx2p(p-1,1)-2*bx2p(p,1)+ bx2p(p+1,1)) (by2p(p-1,1)-2*by2p(p,1)+ by2p(p+1,1)) (bz2p(p-1,1)-2*bz2p(p,1)+ bz2p(p+1,1))])       ;                  
     
end

%----- Euler--Lagrange and Jacobian at the boundary points

p = N-1;


f_b(p-1,1) = 1/2*(bx(p,1)^2 + by(p,1)^2 + bz(p,1)^2 -1);

f_b(N1+p-1,1) = 1/2*(bxp(p,1)^2 + byp(p,1)^2 + bzp(p,1)^2 - tau^2);

%--- closure -----

f_c  = f_c       + [by(p-1,1)*bz(p,1) - bz(p-1,1)*by(p,1) ;
    bz(p-1,1)*bx(p,1) - bx(p-1,1)*bz(p,1) ;
    bx(p-1,1)*by(p,1) - by(p-1,1)*bx(p,1) ] ;

%======================================================================
%====    jacobian ===
ind = 3*p-5:3*p-3;

jac(p-1,ind)           = [bx(p,1) by(p,1) bz(p,1)];

jac(N1+p-1,ind-3)      =-1/h*[bxp(p,1) byp(p,1) bzp(p,1)];
jac(N1+p-1,ind)        = 1/h*[bxp(p,1) byp(p,1) bzp(p,1)];

jac(2*N1+2:2*N1+4,ind) = [ 0                 -(bz(p-1)-bz(p+1))   (by(p-1)-by(p+1))  ;
                           (bz(p-1)-bz(p+1))       0             -(bx(p-1)-bx(p+1))  ;
                          -(by(p-1)-by(p+1))  (bx(p-1)-bx(p+1))           0        ] ; 
                        
jac(2*N1+5,ind)    = h*fac/h^2*([(bx2p(p-1,1)-2*bx2p(p,1)+ bx2p(p+1,1)) (by2p(p-1,1)-2*by2p(p,1)+ by2p(p+1,1)) (bz2p(p-1,1)-2*bz2p(p,1)+ bz2p(p+1,1))])       ;                  


%----- Euler--Lagrange and Jacobian at the boundary points

p = N;


f_b(p-1,1) = 1/2*(bx(p,1)^2 + by(p,1)^2 + bz(p,1)^2 -1);

f_b(N1+p-1,1) = 1/2*(bxp(p,1)^2 + byp(p,1)^2 + bzp(p,1)^2 - tau^2);

%--- closure -----

f_c  = f_c       + [  by(p-1,1)*bz(p,1) - bz(p-1,1)*by(p,1)   ;
    bz(p-1,1)*bx(p,1) - bx(p-1,1)*bz(p,1)   ;
    bx(p-1,1)*by(p,1) - by(p-1,1)*bx(p,1) ] ;


f_b(2*N1+2:2*N1+4,1) =   f_c +     [by(N,1)*bz(N+1,1) - bz(N,1)*by(N+1,1)   ;
    bz(N,1)*bx(N+1,1) - bx(N,1)*bz(N+1,1)   ;
    bx(N,1)*by(N+1,1) - by(N,1)*bx(N+1,1) ] ;


%======================================================================
%====    jacobian ===
ind = 3*p-5:3*p-3;

jac(p-1,ind)           = [bx(p,1) by(p,1) bz(p,1)];

jac(N1+p-1,ind-3)      =-1/h*[bxp(p,1) byp(p,1) bzp(p,1)];
jac(N1+p-1,ind)        = 1/h*[bxp(p,1) byp(p,1) bzp(p,1)];

jac(2*N1+2:2*N1+4,ind) = [ 0                 -(bz(p-1)-bz(p+1))   (by(p-1)-by(p+1))  ;
                           (bz(p-1)-bz(p+1))       0             -(bx(p-1)-bx(p+1))  ;
                          -(by(p-1)-by(p+1))  (bx(p-1)-bx(p+1))           0        ] ; 
                        
jac(2*N1+5,ind)    = h*fac/h^2*([(bx2p(p-1,1)-2*bx2p(p,1)+ bx2p(p+1,1)) (by2p(p-1,1)-2*by2p(p,1)+ by2p(p+1,1)) (bz2p(p-1,1)-2*bz2p(p,1)+ bz2p(p+1,1))])       ;                  

 
%=== derivative at N+1
f_b(2*N1+1,1) = 1/2*(bxp(N+1,1)^2 + byp(N+1,1)^2 + bzp(N+1,1)^2 - tau^2)   ;

jac(2*N1+1,ind) = -1/h*[bxp(N+1,1) byp(N+1,1) bzp(N+1,1)];

%f_b(2*N1+5,1) = sum(sqrt((bx-b0x).^2 + (by-b0y).^2 + (bz-b0z).^2)) - d^2 + sum(sqrt((bx-b2x).^2 + (by-b2y).^2 + (bz-b2z).^2)) - (d1-d)^2;

 
bpp = bx2p.*bx2p + by2p.*by2p + bz2p.*bz2p;

f_b(2*N1+5,1) =   fac/2*(h*(sum(bpp(1:N))) - 16*pi^4*n1^4) - E;


F = f_b;

err = sqrt(sum(f_b.^2)) ;



end
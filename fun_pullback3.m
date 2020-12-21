      function [F,J] = fun_pullback3(x)
 
%==== This is the version that gives same results as if discrete
% We use fun_descent7 to determine the next point in gradient descent
% direction. Then we calculate the distance d between this point and the
% previous point. Then this fold finds minimium energy configuration which
% is at distance d from the previous step at distance d. 

% This code uses continuous spatial formulation. fun_pullback6 is the
% disretized version of this code. 
 
 
global  N  N1  u v w    f_b    len id uni err qd bx2 E count bxN h tau rho lm bxg byg bzg sig
global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p bext J2 fac bxm1 bym1 bzm1 d del

 
 % x = initial_guess_evolution();
      
   %  count = count+1 ;
  %----- reading input and constructing bx by bz ----
    
  bx(2:N,1)     = x(1:3:3*N1-2,1)    ;
  by(2:N,1)     = x(2:3:3*N1-1,1)    ;
  bz(2:N,1)     = x(3:3:3*N1,1)      ;
  
  
  rho           = x(3*N1+1:4*N1,1)   ;
  lm            = x(4*N1+1:5*N1+1,1) ;
  
  u             = x(5*N1+2 ,1)       ;
  v             = x(5*N1+3,1)        ;
  w             = x(5*N1+4,1)        ;
  
 del           = x(5*N1+5,1)        ;
  
  
  bext(1:3,1) = -[bx(N,1);by(N,1);bz(N,1)];
  bext(4:6,1) = -[bx(2,1);by(2,1);bz(2,1)];
  
  om =  [ 0  -w   v  ;
         w   0  -u  ;
        -v   u   0]  ;   %-- \gamma_x tensor
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
    
  %     
   %   lmp = (lm(ig+1,1)-lm(ig-1))/(2*h);
     lmp = (lm(ig,1)-lm(ig-1,1))/(h);
   
          
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
    
   f_c = [0;0;0];
    
  
  %----- Euler--Lagrange and Jacobian at the boundary points 
  fac = 0.0001;
  
  %== For studying the effect of aspect ratio on the lagrange multipliers
  %===
  
  %fac = asinh(tau*pi*sig)/(4*pi^3*tau^3);
   
  
  p = 2;
  
           
        gm = om*[bx(p+1)-bx(p-1);by(p+1)-by(p-1);bz(p+1)-bz(p-1)];
        gmx = gm(1);
        gmy = gm(2);
        gmz = gm(3);
        
        
        f_b(3*p-5,1) = fac*bx4p(p-1,1) + lmp(p-1)*bxp(p,1) + lm(p)*bx2p(p,1) - rho(p-1)*bx(p,1) + gmx + del*(bx(p,1)-bxm1(p,1));
        f_b(3*p-4,1) = fac*by4p(p-1,1) + lmp(p-1)*byp(p,1) + lm(p)*by2p(p,1) - rho(p-1)*by(p,1) + gmy + del*(by(p,1)-bym1(p,1));  
        f_b(3*p-3,1) = fac*bz4p(p-1,1) + lmp(p-1)*bzp(p,1) + lm(p)*bz2p(p,1) - rho(p-1)*bz(p,1) + gmz + del*(bz(p,1)-bzm1(p,1));
        
        
        f_b(3*N1+p-1,1) = 1/2*(bx(p,1)^2 + by(p,1)^2 + bz(p,1)^2 -1);
        
        f_b(4*N1+p-1,1) = 1/2*(bxp(p,1)^2 + byp(p,1)^2 + bzp(p,1)^2 - tau^2); 
        
       %--- closure -----
%        
        f_c  = f_c       + [by(p-1,1)*bz(p,1) - bz(p-1,1)*by(p,1)    ;
                            bz(p-1,1)*bx(p,1) - bx(p-1,1)*bz(p,1)    ;
                            bx(p-1,1)*by(p,1) - by(p-1,1)*bx(p,1) ] ; 
%   
  %==================== Jacobian entry ================================
  
  ind = 3*p-5:3*p-3;
  
  jac(ind,ind)     = id*(6*fac/h^4 + lmp(p-1,1)/h -2*lm(p)/h^2-rho(p-1) + del)    ;  
  jac(ind,ind+3)   = id*(-4*fac/h^4   +lm(p)/h^2)  +om                              ; 
  jac(ind,ind+6)   = id*1*fac/h^4                                                ;
  
  jac(ind,3*N1-2:3*N1) = -id*1*fac/h^4  ;
  

  
%   jac(ind,ind)     = id*(1*fac/h^4 + lmp(p-1,1)/h -2*lm(p)/h^2-rho(p-1))   ;  
%   jac(ind,ind+3)   = id*(-4*fac/h^4   +lm(p)/h^2)     + om                 ; 
%   jac(ind,ind+6)   = id*6*fac/h^4                                          ;
%   jac(ind,ind+9)   = id*(-4)*fac/h^4                                       ;
%   jac(ind,ind+12)  = id*1*fac/h^4                                          ;
%   
  %--- Constraint gradient ---
  
  %--- unit vector constraint --
  jac(3*N1+p-1,ind) = [bx(p) by(p) bz(p)];
  %--- unit speed constraint --
  jac(4*N1+p-1,ind)  =  1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)];
   
   %---- closure ----
  jac(5*N1+2:5*N1+4,ind) =   [ 0                 -(bz(p-1)-bz(p+1))   (by(p-1)-by(p+1))  ;
                               (bz(p-1)-bz(p+1))       0             -(bx(p-1)-bx(p+1))  ;
                              -(by(p-1)-by(p+1))  (bx(p-1)-bx(p+1))           0        ] ;  
 
  jac(5*N1+5,ind) = [(bx(p,1)-bxm1(p,1)) (by(p,1)-bym1(p,1)) (bz(p,1)-bzm1(p,1))]     ;                  
   %-- differentiation with rho
   
   jac(ind,3*N1+p-1) = -[bx(p) ;by(p) ;bz(p)];
  
   %-- differentiation with lambda 
   jac(ind,4*N1+p-1)   = -1/h*[bxp(p) ;byp(p) ;bzp(p)]; 
   jac(ind,4*N1+p)     =  1/h*[bxp(p) ;byp(p) ;bzp(p)] +  [bx2p(p) ;by2p(p) ;bz2p(p)]; 

  %-- differentiation with u v w 
   jac(ind,5*N1+2:5*N1+4)       = - [ 0                 -(bz(p-1)-bz(p+1))   (by(p-1)-by(p+1))  ;
                                    (bz(p-1)-bz(p+1))       0             -(bx(p-1)-bx(p+1))  ;
                                   -(by(p-1)-by(p+1))  (bx(p-1)-bx(p+1))           0        ] ;  
    
   jac(ind,5*N1+5) =   [(bx(p,1)-bxm1(p,1)) (by(p,1)-bym1(p,1)) (bz(p,1)-bzm1(p,1))]'     ;                                     
    %======================================================================
  
  %----- Euler--Lagrange and Jacobian at the boundary points 
  
  p = 3;
   
        
        gm = om*[bx(p+1)-bx(p-1);by(p+1)-by(p-1);bz(p+1)-bz(p-1)];
        gmx = gm(1);
        gmy = gm(2);
        gmz = gm(3);
% %         
        
        f_b(3*p-5,1) = fac*bx4p(p-1,1) + lmp(p-1)*bxp(p,1) + lm(p)*bx2p(p,1) - rho(p-1)*bx(p,1) + gmx + del*(bx(p,1)-bxm1(p,1));
        f_b(3*p-4,1) = fac*by4p(p-1,1) + lmp(p-1)*byp(p,1) + lm(p)*by2p(p,1) - rho(p-1)*by(p,1) + gmy + del*(by(p,1)-bym1(p,1));  
        f_b(3*p-3,1) = fac*bz4p(p-1,1) + lmp(p-1)*bzp(p,1) + lm(p)*bz2p(p,1) - rho(p-1)*bz(p,1) + gmz + del*(bz(p,1)-bzm1(p,1));
        
        
        
        f_b(3*N1+p-1,1) = 1/2*(bx(p,1)^2 + by(p,1)^2 + bz(p,1)^2 -1);
        
        f_b(4*N1+p-1,1) = 1/2*(bxp(p,1)^2 + byp(p,1)^2 + bzp(p,1)^2 - tau^2); 
        
       %--- closure -----
       
        f_c  = f_c       + [by(p-1,1)*bz(p,1) - bz(p-1,1)*by(p,1) ;
                            bz(p-1,1)*bx(p,1) - bx(p-1,1)*bz(p,1) ;
                            bx(p-1,1)*by(p,1) - by(p-1,1)*bx(p,1) ] ; 
                    
  %==================== Jacobian entry ================================
  
  ind = 3*p-5:3*p-3;
  
  jac(ind,ind-3)   = id*(-4*fac/h^4 - lmp(p-1,1)/h + lm(p)/h^2) - om        ;
  jac(ind,ind)     = id*(6*fac/h^4 + lmp(p-1,1)/h -2*lm(p)/h^2-rho(p-1) + del)        ;  
  jac(ind,ind+3)   = id*(-4*fac/h^4   +lm(p)/h^2)    +om  ; 
  jac(ind,ind+6)   = id*1*fac/h^4                                             ;
  
  
  
  %--- Constraint gradient ---
  
  %--- unit vector constraint --
  jac(3*N1+p-1,ind) = [bx(p) by(p) bz(p)];
  %--- unit speed constraint --
  jac(4*N1+p-1,ind-3) = -1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)];
  jac(4*N1+p-1,ind)   =  1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)];
   
   %---- closure ----
  jac(5*N1+2:5*N1+4,ind) =   [ 0                 -(bz(p-1)-bz(p+1))   (by(p-1)-by(p+1))  ;
                               (bz(p-1)-bz(p+1))       0             -(bx(p-1)-bx(p+1))  ;
                              -(by(p-1)-by(p+1))  (bx(p-1)-bx(p+1))           0        ] ;  
   
  jac(5*N1+5,ind) = [(bx(p,1)-bxm1(p,1)) (by(p,1)-bym1(p,1)) (bz(p,1)-bzm1(p,1))]     ;                  

   %-- differentiation with rho
   
   jac(ind,3*N1+p-1) = -[bx(p) ;by(p) ;bz(p)];
  
   %-- differentiation with lambda 
   jac(ind,4*N1+p-1)   = -1/h*[bxp(p) ;byp(p) ;bzp(p)]; 
   jac(ind,4*N1+p)     =  1/h*[bxp(p) ;byp(p) ;bzp(p)] +  [bx2p(p) ;by2p(p) ;bz2p(p)]; 

  %-- differentiation with u v w 
   jac(ind,5*N1+2:5*N1+4)       = - [ 0                 -(bz(p-1)-bz(p+1))   (by(p-1)-by(p+1))  ;
                                    (bz(p-1)-bz(p+1))       0             -(bx(p-1)-bx(p+1))  ;
                                      -(by(p-1)-by(p+1))  (bx(p-1)-bx(p+1))           0        ] ;  
  
  jac(ind,5*N1+5) =   [(bx(p,1)-bxm1(p,1)) (by(p,1)-bym1(p,1)) (bz(p,1)-bzm1(p,1))]'     ;                                     
                                       
    %======================================================================
  
  
  
    for p = 4:N-2
          
         
        gm = om*[bx(p+1)-bx(p-1);by(p+1)-by(p-1);bz(p+1)-bz(p-1)];
        gmx = gm(1);
        gmy = gm(2);
        gmz = gm(3);
        
        
        f_b(3*p-5,1) = fac*bx4p(p-1,1) + lmp(p-1)*bxp(p,1) + lm(p)*bx2p(p,1) - rho(p-1)*bx(p,1) + gmx + del*(bx(p,1)-bxm1(p,1));
        f_b(3*p-4,1) = fac*by4p(p-1,1) + lmp(p-1)*byp(p,1) + lm(p)*by2p(p,1) - rho(p-1)*by(p,1) + gmy + del*(by(p,1)-bym1(p,1));  
        f_b(3*p-3,1) = fac*bz4p(p-1,1) + lmp(p-1)*bzp(p,1) + lm(p)*bz2p(p,1) - rho(p-1)*bz(p,1) + gmz + del*(bz(p,1)-bzm1(p,1));
        
        
        
        f_b(3*N1+p-1,1) = 1/2*(bx(p,1)^2 + by(p,1)^2 + bz(p,1)^2 -1);
        
        f_b(4*N1+p-1,1) = 1/2*(bxp(p,1)^2 + byp(p,1)^2 + bzp(p,1)^2 - tau^2); 
        
       %--- closure -----
       
          f_c  = f_c       + [by(p-1,1)*bz(p,1) - bz(p-1,1)*by(p,1) ;
                              bz(p-1,1)*bx(p,1) - bx(p-1,1)*bz(p,1) ;
                              bx(p-1,1)*by(p,1) - by(p-1,1)*bx(p,1) ] ;  
       
                       
  %==================== Jacobian entry ================================
  
  ind = 3*p-5:3*p-3;
  
  jac(ind,ind-6)   = id*fac/h^4                                                  ;
  jac(ind,ind-3)   = id*(-4*fac/h^4 - lmp(p-1,1)/h + lm(p)/h^2) - om            ;
  jac(ind,ind)     = id*(6*fac/h^4 + lmp(p-1,1)/h -2*lm(p)/h^2-rho(p-1) + del)     ;  
  jac(ind,ind+3)   = id*(-4*fac/h^4   +lm(p)/h^2)    + om                            ; 
  jac(ind,ind+6)   = id*1*fac/h^4                                                ;
  
  
  
  %--- Constraint gradient ---
  
  %--- unit vector constraint --
  jac(3*N1+p-1,ind) = [bx(p) by(p) bz(p)];
  %--- unit speed constraint --
  jac(4*N1+p-1,ind-3) = -1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)];
  jac(4*N1+p-1,ind)   =  1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)];
   
   %---- closure ----
  jac(5*N1+2:5*N1+4,ind) =   [ 0                 -(bz(p-1)-bz(p+1))   (by(p-1)-by(p+1))  ;
                               (bz(p-1)-bz(p+1))       0             -(bx(p-1)-bx(p+1))  ;
                              -(by(p-1)-by(p+1))  (bx(p-1)-bx(p+1))           0        ] ;  
   
  jac(5*N1+5,ind) = [(bx(p,1)-bxm1(p,1)) (by(p,1)-bym1(p,1)) (bz(p,1)-bzm1(p,1))]     ;                  

   %-- differentiation with rho
   
   jac(ind,3*N1+p-1) = -[bx(p) ;by(p) ;bz(p)];
  
   %-- differentiation with lambda 
   jac(ind,4*N1+p-1)   = -1/h*[bxp(p) ;byp(p) ;bzp(p)]; 
   jac(ind,4*N1+p)     =  1/h*[bxp(p) ;byp(p) ;bzp(p)] +  [bx2p(p) ;by2p(p) ;bz2p(p)]; 

  %-- differentiation with u v w 
     jac(ind,5*N1+2:5*N1+4)       = - [ 0                 -(bz(p-1)-bz(p+1))   (by(p-1)-by(p+1))  ;
                                    (bz(p-1)-bz(p+1))       0             -(bx(p-1)-bx(p+1))  ;
                                      -(by(p-1)-by(p+1))  (bx(p-1)-bx(p+1))           0        ] ;  

      jac(ind,5*N1+5) =   [(bx(p,1)-bxm1(p,1)) (by(p,1)-bym1(p,1)) (bz(p,1)-bzm1(p,1))]'     ;                                     
                                       
    %======================================================================
    end
     
 %----- Euler--Lagrange and Jacobian at the boundary points 
  
  p = N-1;
 
        
        gm = om*[bx(p+1)-bx(p-1);by(p+1)-by(p-1);bz(p+1)-bz(p-1)];
        gmx = gm(1);
        gmy = gm(2);
        gmz = gm(3);
        
        
        f_b(3*p-5,1) = fac*bx4p(p-1,1) + lmp(p-1)*bxp(p,1) + lm(p)*bx2p(p,1) - rho(p-1)*bx(p,1) + gmx + del*(bx(p,1)-bxm1(p,1));
        f_b(3*p-4,1) = fac*by4p(p-1,1) + lmp(p-1)*byp(p,1) + lm(p)*by2p(p,1) - rho(p-1)*by(p,1) + gmy + del*(by(p,1)-bym1(p,1));  
        f_b(3*p-3,1) = fac*bz4p(p-1,1) + lmp(p-1)*bzp(p,1) + lm(p)*bz2p(p,1) - rho(p-1)*bz(p,1) + gmz + del*(bz(p,1)-bzm1(p,1));
        
        
        
        f_b(3*N1+p-1,1) = 1/2*(bx(p,1)^2 + by(p,1)^2 + bz(p,1)^2 -1);
        
        f_b(4*N1+p-1,1) = 1/2*(bxp(p,1)^2 + byp(p,1)^2 + bzp(p,1)^2 - tau^2); 
        
       %--- closure -----
       
        f_c  = f_c       + [by(p-1,1)*bz(p,1) - bz(p-1,1)*by(p,1) ;
                            bz(p-1,1)*bx(p,1) - bx(p-1,1)*bz(p,1) ;
                            bx(p-1,1)*by(p,1) - by(p-1,1)*bx(p,1) ] ; 
                    
  %==================== Jacobian entry ================================
  
  ind = 3*p-5:3*p-3;
  
  jac(ind,ind-6)   = id*fac/h^4                                                  ;
  jac(ind,ind-3)   = id*(-4*fac/h^4 - lmp(p-1,1)/h + lm(p)/h^2) - om            ;
  jac(ind,ind)     = id*(6*fac/h^4 + lmp(p-1,1)/h -2*lm(p)/h^2-rho(p-1) + del)   ;  
  jac(ind,ind+3)   = id*(-4*fac/h^4   +lm(p)/h^2)             +om                   ; 
   
  
  
  %--- Constraint gradient ---
  
  %--- unit vector constraint --
  jac(3*N1+p-1,ind) = [bx(p) by(p) bz(p)];
  %--- unit speed constraint --
  jac(4*N1+p-1,ind-3) = -1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)];
  jac(4*N1+p-1,ind)   =  1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)];
   
   %---- closure ----
  jac(5*N1+2:5*N1+4,ind) =   [ 0                 -(bz(p-1)-bz(p+1))   (by(p-1)-by(p+1))  ;
                               (bz(p-1)-bz(p+1))       0             -(bx(p-1)-bx(p+1))  ;
                              -(by(p-1)-by(p+1))  (bx(p-1)-bx(p+1))           0        ] ;  
 
  jac(5*N1+5,ind) = [(bx(p,1)-bxm1(p,1)) (by(p,1)-bym1(p,1)) (bz(p,1)-bzm1(p,1))]     ;                  

   %-- differentiation with rho
   
   jac(ind,3*N1+p-1) = -[bx(p) ;by(p) ;bz(p)];
  
   %-- differentiation with lambda 
   jac(ind,4*N1+p-1)   = -1/h*[bxp(p) ;byp(p) ;bzp(p)]; 
   jac(ind,4*N1+p)     =  1/h*[bxp(p) ;byp(p) ;bzp(p)] +  [bx2p(p) ;by2p(p) ;bz2p(p)]; 

   %-- differentiation with u v w 
   jac(ind,5*N1+2:5*N1+4)       = - [ 0                 -(bz(p-1)-bz(p+1))   (by(p-1)-by(p+1))  ;
                                    (bz(p-1)-bz(p+1))       0             -(bx(p-1)-bx(p+1))  ;
                                      -(by(p-1)-by(p+1))  (bx(p-1)-bx(p+1))           0        ] ;  

   jac(ind,5*N1+5) =   [(bx(p,1)-bxm1(p,1)) (by(p,1)-bym1(p,1)) (bz(p,1)-bzm1(p,1))]'     ;                                     
                                       
    %======================================================================
   
 %----- Euler--Lagrange and Jacobian at the boundary points 
  
 p = N;
 
         gm = om*[bx(p+1)-bx(p-1);by(p+1)-by(p-1);bz(p+1)-bz(p-1)];
         gmx = gm(1);
         gmy = gm(2);
         gmz = gm(3);
 
 
        f_b(3*p-5,1) = fac*bx4p(p-1,1) + lmp(p-1)*bxp(p,1) + lm(p)*bx2p(p,1) - rho(p-1)*bx(p,1) + gmx + del*(bx(p,1)-bxm1(p,1));
        f_b(3*p-4,1) = fac*by4p(p-1,1) + lmp(p-1)*byp(p,1) + lm(p)*by2p(p,1) - rho(p-1)*by(p,1) + gmy + del*(by(p,1)-bym1(p,1));  
        f_b(3*p-3,1) = fac*bz4p(p-1,1) + lmp(p-1)*bzp(p,1) + lm(p)*bz2p(p,1) - rho(p-1)*bz(p,1) + gmz + del*(bz(p,1)-bzm1(p,1));
        
 
 
 f_b(3*N1+p-1,1) = 1/2*(bx(p,1)^2 + by(p,1)^2 + bz(p,1)^2 -1)         ;
 
 f_b(4*N1+p-1,1) = 1/2*(bxp(p,1)^2 + byp(p,1)^2 + bzp(p,1)^2 - tau^2) ;
 
 %--- closure -----
 
 f_c  = f_c       + [  by(p-1,1)*bz(p,1) - bz(p-1,1)*by(p,1)   ;
                       bz(p-1,1)*bx(p,1) - bx(p-1,1)*bz(p,1)   ;
                       bx(p-1,1)*by(p,1) - by(p-1,1)*bx(p,1) ] ;
 
 
 f_b(5*N1+2:5*N1+4,1) =  0.*[0;-bz(N);by(N)] + f_c +     [by(N,1)*bz(N+1,1) - bz(N,1)*by(N+1,1)   ;
                                                          bz(N,1)*bx(N+1,1) - bx(N,1)*bz(N+1,1)   ;
                                                          bx(N,1)*by(N+1,1) - by(N,1)*bx(N+1,1) ] ;
 
 
 f_b(5*N1+1,1) = 1/2*(bxp(N+1,1)^2 + byp(N+1,1)^2 + bzp(N+1,1)^2 - tau^2)   ;
 
 %==================== Jacobian entry ================================
 
 ind = 3*p-5:3*p-3;
 
 jac(ind,ind-6)   = id*fac/h^4                                                  ;
 jac(ind,ind-3)   = id*(-4*fac/h^4 - lmp(p-1,1)/h + lm(p)/h^2) - om           ;
 jac(ind,ind)     = id*(6*fac/h^4 + lmp(p-1,1)/h -2*lm(p)/h^2-rho(p-1) + del)    ;
 
 jac(ind,1:3)     = -id*fac/h^4;
    
  
  %jac(5*N1+2:5*N1+4,3*N1-2:3*N1) = .5*[0 0 0; 0 0 -1;0 1 0];
  %--- Constraint gradient ---
  
  %--- unit vector constraint --
  jac(3*N1+p-1,ind) = [bx(p) by(p) bz(p)];
  %--- unit speed constraint --
  jac(4*N1+p-1,ind-3) = -1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)];     % b'(N)
  jac(4*N1+p-1,ind)   =  1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)];
  
  jac(5*N1+1,ind)   = -1/h*[bxp(N+1,1) byp(N+1,1) bzp(N+1,1)];  %b'(N+1)
   
   %---- closure ----
  jac(5*N1+2:5*N1+4,ind) =    [ 0                 -(bz(p-1)-bz(p+1))   (by(p-1)-by(p+1))  ;
                               (bz(p-1)-bz(p+1))       0             -(bx(p-1)-bx(p+1))  ;
                              -(by(p-1)-by(p+1))  (bx(p-1)-bx(p+1))           0        ] ;  
  jac(5*N1+5,ind) = [(bx(p,1)-bxm1(p,1)) (by(p,1)-bym1(p,1)) (bz(p,1)-bzm1(p,1))]     ;                  

   %-- differentiation with rho
   
   jac(ind,3*N1+p-1) = -[bx(p) ;by(p) ;bz(p)];
  
   %-- differentiation with lambda 
   jac(ind,4*N1+p-1)   = -1/h*[bxp(p) ;byp(p) ;bzp(p)]; 
   jac(ind,4*N1+p)     =  1/h*[bxp(p) ;byp(p) ;bzp(p)] +  [bx2p(p) ;by2p(p) ;bz2p(p)]; 

   %-- differentiation with u v w 
   jac(ind,5*N1+2:5*N1+4)       = - [ 0                 -(bz(p-1)-bz(p+1))   (by(p-1)-by(p+1))  ;
                                    (bz(p-1)-bz(p+1))       0             -(bx(p-1)-bx(p+1))  ;
                                      -(by(p-1)-by(p+1))  (bx(p-1)-bx(p+1))           0        ] ;  

   jac(ind,5*N1+5) =   [(bx(p,1)-bxm1(p,1)) (by(p,1)-bym1(p,1)) (bz(p,1)-bzm1(p,1))]'     ;                                     
                               
    jac(:,5*N1+2:5*N1+4) = -jac(:,5*N1+2:5*N1+4);
                                          
    %======================================================================
  
    %=================  Deleting columns and rows with respect to entries for
  %bz(2) and bz(N)
   
 
% %   

f_b(5*N1+5,1) = 1/2*(sum((bx-bxm1).^2 + (by-bym1).^2 + (bz-bzm1).^2) - d^2);

  
 F = f_b;
  J = jac;
 err = sqrt(sum(f_b.^2)) ;  
 
      end
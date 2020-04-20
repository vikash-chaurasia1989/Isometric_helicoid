     function [F,J] = fun_jacobian(x)
 

 %--- in this file, we solve for minimizing bending energy of the ruled surface --
   
   %  b'''' + lm'*b' + lm*b'' - rho b + gm x b' = 0;
   %
   % --- constraints --
   %   |b| = 1 and |b|' = \tau
   %--- This is the final version of the analytical jacobian for finite
   %element discretization of the continouous formulation 
   
 
global  N  N1  u v w    f_b    len id uni err qd bx2 E count bxN h tau rho lm bxg byg bzg
global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p bext J2

 
 %  x = initial_guess3();
      
     count = count+1 ;
  %----- reading input and constructing bx by bz ----
  bx(2,1)       = x(1,1);
  by(2,1)       = x(2,1);
  
  bx(3:N,1)     = x(3:3:3*N1-3,1)    ;
  by(3:N,1)     = x(4:3:3*N1-2,1)    ;
  bz(3:N-1,1)   = x(5:3:3*N1-4,1)    ;
  
  
  rho           = x(3*N1-1:4*N1-2,1) ;
  lm            = x(4*N1-1:5*N1,1)   ;
  
  u             = x(5*N1+1,1)        ;
  v             = x(5*N1+2,1)        ;
  w             = x(5*N1+3,1)        ;
  
  bext          = x(5*N1+4:5*N1+9,1) ;
 
 
    
  om = [ 0  -w   v  ;
         w   0  -u  ;
        -v   u   0] ;   %-- \gamma_x tensor
 %---- derivative calculation for b2 --- bN  ---
 
   ig = 2:N;
   
   
   %--  b'   (O(h^2))----
   
   %-- Central difference
   %---O(h^2)---
   
     bxp(1,1) = (bx(2,1)-bext(1,1))/(2*h);
     byp(1,1) = (by(2,1)-bext(2,1))/(2*h);
     bzp(1,1) = (bz(2,1)-bext(3,1))/(2*h);
 
     bxp(ig,1) = (bx(ig+1,1)-bx(ig-1,1))/(2*h);
     byp(ig,1) = (by(ig+1,1)-by(ig-1,1))/(2*h);
     bzp(ig,1) = (bz(ig+1,1)-bz(ig-1,1))/(2*h);
     
 
    bxp(N+1,1) = (bext(4,1)-bx(N,1))/(2*h);
    byp(N+1,1) = (bext(5,1)-by(N,1))/(2*h);
    bzp(N+1,1) = (bext(6,1)-bz(N,1))/(2*h);
    
             %  lmp = (lm(ig+1,1)-lm(ig))/(1*h);
      
              
  %=== O(h)====
%      bxp(1,1) = (bx(1,1)-bext(1,1))/(h);
%      byp(1,1) = (by(1,1)-bext(2,1))/(h);
%      bzp(1,1) = (bz(1,1)-bext(3,1))/(h);
%  
%      bxp(ig,1) = (bx(ig,1)-bx(ig-1,1))/(h);
%      byp(ig,1) = (by(ig,1)-by(ig-1,1))/(h);
%      bzp(ig,1) = (bz(ig,1)-bz(ig-1,1))/(h);
%      
%  
%     bxp(N+1,1) = (bx(N+1,1)-bx(N,1))/(h);
%     byp(N+1,1) = (by(N+1,1)-by(N,1))/(h);
%     bzp(N+1,1) = (bz(N+1,1)-bz(N,1))/(h);
% %     
      lmp = (lm(ig+1,1)-lm(ig-1))/(2*h);
   %    lmp = (lm(ig,1)-lm(ig-1,1))/(h);
      
  
          
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
  
  p = 2;
  
        gmx = v*bzp(p,1) - w*byp(p,1);
        gmy = w*bxp(p,1) - u*bzp(p,1);
        gmz = u*byp(p,1) - v*bxp(p,1);
        
        f_b(3*p-5,1) = .001*bx4p(p-1,1) + lmp(p-1)*bxp(p,1) + lm(p)*bx2p(p,1) - rho(p-1)*bx(p,1) + gmx;
        f_b(3*p-4,1) = .001*by4p(p-1,1) + lmp(p-1)*byp(p,1) + lm(p)*by2p(p,1) - rho(p-1)*by(p,1) + gmy;
        f_b(3*p-3,1) = .001*bz4p(p-1,1) + lmp(p-1)*bzp(p,1) + lm(p)*bz2p(p,1) - rho(p-1)*bz(p,1) + gmz;
        
        
        f_b(3*N1+p-1,1) = 1/2*(bx(p,1)^2 + by(p,1)^2 + bz(p,1)^2 -1);
        
        f_b(4*N1+p,1) = 1/2*(bxp(p,1)^2 + byp(p,1)^2 + bzp(p,1)^2 - tau^2); 
        
       %--- closure -----
       
%         f_c  = f_c       + [by(p-1,1)*bz(p,1) - bz(p-1,1)*by(p,1)    ;
%                                  bz(p-1,1)*bx(p,1) - bx(p-1,1)*bz(p,1)    ;
%                                  bx(p-1,1)*by(p,1) - by(p-1,1)*bx(p,1) ] ; 
      
        f_c  = f_c     + [by(p,1)*bzp(p,1) - bz(p,1)*byp(p,1) ;
                              bz(p,1)*bxp(p,1) - bx(p,1)*bzp(p,1) ;
                              bx(p,1)*byp(p,1) - by(p,1)*bxp(p,1) ] ;
  %==================== Jacobian entry ================================
  
  ind = 3*p-5:3*p-3;
  
  
  jac(ind,ind)     = id*(6*.001/h^4-2*lm(p)/h^2-rho(p-1))                    ;  
  jac(ind,ind+3) = id*(-4*.001/h^4 +lmp(p-1)/(2*h) +lm(p)/h^2) + om/(2*h)    ; 
  jac(ind,ind+6) = id*1*.001/h^4                                             ;
   
  jac(ind,5*N1+6:5*N1+8) = id*1*.001/h^4                                     ;
  
  %--- Constraint gradient ---
  
  %--- unit vector constraint --
  jac(3*N1+p-1,ind) = [bx(p) by(p) bz(p)];
  %--- unit speed constraint --
   jac(4*N1+p-1,5*N1+6:5*N1+8) = -1/(2*h)*[bxp(p-1,1) byp(p-1,1) bzp(p-1,1)];
   jac(4*N1+p-1,ind)                     =  1/(2*h)*[bxp(p-1,1) byp(p-1,1) bzp(p-1,1)];
   
   %-- differentiation with rho
   
   jac(ind,3*N1+p-1) = -[bx(p) ;by(p) ;bz(p)];
  
   %-- differentiation with lambda 
     jac(ind,4*N1+p-1:4*N1+p+1)   = [ -1/(2*h)*bxp(p)   bx2p(p) 1/(2*h)*bxp(p)  ; 
                                                         -1/(2*h)*byp(p)    by2p(p) 1/(2*h)*byp(p)  ;
                                                         -1/(2*h)*bzp(p)   bz2p(p)  1/(2*h)*bzp(p) ] ;
  %-- differentiation with u v w 
     jac(ind,5*N1+3:5*N1+5)       = -[   0          -bzp(p)    byp(p)  ;
                                                      bzp(p)       0        -bxp(p) ;
                                                    -byp(p)      bxp(p)      0    ] ;
    
  %---- closure ----
%        
%     
%   jac(5*N1+3:5*N1+5,ind) =   [              0                  -(bz(p-1)-bz(p+1))     (by(p-1)-by(p+1))    ;
%                                                (bz(p-1)-bz(p+1))       0                           -(bx(p-1)-bx(p+1))   ;
%                                              -(by(p-1)-by(p+1))       (bx(p-1)-bx(p+1))      0                        ] ;  
                                         
   jac(5*N1+3:5*N1+5,ind) =   [              0                  -(bz(p-1)/2-bz(p+1))     (by(p-1)/2-by(p+1))    ;
                                               (bz(p-1)/2-bz(p+1))       0                           -(bx(p-1)/2-bx(p+1))   ;
                                             -(by(p-1)/2-by(p+1))       (bx(p-1)/2-bx(p+1))      0                        ] ;                                          
  %======================================================================
  
  %----- Euler--Lagrange and Jacobian at the boundary points 
  
  p = 3;
  
        gmx = v*bzp(p,1) - w*byp(p,1);
        gmy = w*bxp(p,1) - u*bzp(p,1); 
        gmz = u*byp(p,1) - v*bxp(p,1);
        
        f_b(3*p-5,1) = .001*bx4p(p-1,1) + lmp(p-1)*bxp(p,1) + lm(p)*bx2p(p,1) - rho(p-1)*bx(p,1) + gmx;
        f_b(3*p-4,1) = .001*by4p(p-1,1) + lmp(p-1)*byp(p,1) + lm(p)*by2p(p,1) - rho(p-1)*by(p,1) + gmy;
        f_b(3*p-3,1) = .001*bz4p(p-1,1) + lmp(p-1)*bzp(p,1) + lm(p)*bz2p(p,1) - rho(p-1)*bz(p,1) + gmz;
        
        
        f_b(3*N1+p-1,1) = 1/2*(bx(p,1)^2 + by(p,1)^2 + bz(p,1)^2 -1);
        
        f_b(4*N1+p,1) = 1/2*(bxp(p,1)^2 + byp(p,1)^2 + bzp(p,1)^2 - tau^2); 
        
       %--- closure -----
       
%         f_c  = f_c       + [by(p-1,1)*bz(p,1) - bz(p-1,1)*by(p,1) ;
%                             bz(p-1,1)*bx(p,1) - bx(p-1,1)*bz(p,1) ;
%                             bx(p-1,1)*by(p,1) - by(p-1,1)*bx(p,1) ] ; 
           f_c  = f_c     + [by(p,1)*bzp(p,1) - bz(p,1)*byp(p,1) ;
                              bz(p,1)*bxp(p,1) - bx(p,1)*bzp(p,1) ;
                              bx(p,1)*byp(p,1) - by(p,1)*bxp(p,1) ] ;                    
  %==================== Jacobian entry ================================
  
  ind = 3*p-5:3*p-3;
  
  
  jac(ind,ind-3)   = id*(-4*.001/h^4 - lmp(p-1)/(2*h) + lm(p)/h^2)  - om/(2*h) ;
  jac(ind,ind)     = id*( 6*.001/h^4 - 2*lm(p)/h^2-rho(p-1))                 ; 
  jac(ind,ind+3)   = id*(-4*.001/h^4 + lmp(p-1)/(2*h) + lm(p)/h^2) +om/(2*h) ;  
  jac(ind,ind+6)   = id*1*.001/h^4  ;  
  
   
   
  %======================================================================
   %--- Constraint gradient ---
  
  %--- unit vector constraint --
  jac(3*N1+p-1,ind) = [bx(p) by(p) bz(p)];
  %--- unit speed constraint --
  jac(4*N1+p-1,ind)           =  1/(2*h)*[bxp(p-1,1) byp(p-1,1) bzp(p-1,1)];
  
     %-- differentiation with rho
   
   jac(ind,3*N1+p-1) = -[bx(p) ;by(p) ;bz(p)];
  
   %-- differentiation with lambda 
   %-- differentiation with lambda 
     jac(ind,4*N1+p-1:4*N1+p+1)   = [ -1/(2*h)*bxp(p)   bx2p(p) 1/(2*h)*bxp(p)  ; 
                                                         -1/(2*h)*byp(p)    by2p(p) 1/(2*h)*byp(p)  ;
                                                         -1/(2*h)*bzp(p)   bz2p(p)  1/(2*h)*bzp(p) ] ;
    %-- differentiation with u v w 
     jac(ind,5*N1+3:5*N1+5)       = -[   0          -bzp(p)    byp(p)  ;
                                                      bzp(p)       0        -bxp(p) ;
                                                    -byp(p)      bxp(p)      0    ] ;
  %---- closure ----
       
    
  jac(5*N1+3:5*N1+5,ind) =  [              0                  -(bz(p-1)-bz(p+1))     (by(p-1)-by(p+1))    ;
                                               (bz(p-1)-bz(p+1))       0                           -(bx(p-1)-bx(p+1))   ;
                                             -(by(p-1)-by(p+1))       (bx(p-1)-bx(p+1))      0                        ] ;  
  %======================================================================
  
  
  
    for p = 4:N-2
          
        
        
        gmx = v*bzp(p,1) - w*byp(p,1);
        gmy = w*bxp(p,1) - u*bzp(p,1);
        gmz = u*byp(p,1) - v*bxp(p,1);
        
        f_b(3*p-5,1) = .001*bx4p(p-1,1) + lmp(p-1)*bxp(p,1) + lm(p)*bx2p(p,1) - rho(p-1)*bx(p,1) + gmx;
        f_b(3*p-4,1) = .001*by4p(p-1,1) + lmp(p-1)*byp(p,1) + lm(p)*by2p(p,1) - rho(p-1)*by(p,1) + gmy;
        f_b(3*p-3,1) = .001*bz4p(p-1,1) + lmp(p-1)*bzp(p,1) + lm(p)*bz2p(p,1) - rho(p-1)*bz(p,1) + gmz;
        
        
        f_b(3*N1+p-1,1) = 1/2*(bx(p,1)^2 + by(p,1)^2 + bz(p,1)^2 -1);
        
        f_b(4*N1+p,1) = 1/2*(bxp(p,1)^2 + byp(p,1)^2 + bzp(p,1)^2 - tau^2); 
        
       %--- closure -----
       
%         f_c  = f_c       + [by(p-1,1)*bz(p,1) - bz(p-1,1)*by(p,1) ;
%                             bz(p-1,1)*bx(p,1) - bx(p-1,1)*bz(p,1) ;
%                             bx(p-1,1)*by(p,1) - by(p-1,1)*bx(p,1) ] ;  
       
        f_c  = f_c     + [by(p,1)*bzp(p,1) - bz(p,1)*byp(p,1) ;
                              bz(p,1)*bxp(p,1) - bx(p,1)*bzp(p,1) ;
                              bx(p,1)*byp(p,1) - by(p,1)*bxp(p,1) ] ;                           
  %==================== Jacobian entry ================================
  
  ind = 3*p-5:3*p-3;
  
  jac(ind,ind-6)   = id*1*.001/h^4;
  jac(ind,ind-3)   = id*(-4*.001/h^4 - lmp(p-1)/(2*h) + lm(p)/h^2) -om/(2*h) ;
  jac(ind,ind)       = id*( 6*.001/h^4 - 2*lm(p)/h^2-rho(p-1))                 ; 
  jac(ind,ind+3)   = id*(-4*.001/h^4 + lmp(p-1)/(2*h) + lm(p)/h^2) +om/(2*h) ;                                            ;
  jac(ind,ind+6)   = id*1*.001/h^4;
   
  %======================================================================
   %--- Constraint gradient ---
  
  %--- unit vector constraint --
  jac(3*N1+p-1,ind) = [bx(p) by(p) bz(p)];
  %--- unit speed constraint --
  
  jac(4*N1+p-1,ind-6)  =  -1/(2*h)*[bxp(p-1,1) byp(p-1,1) bzp(p-1,1)]; 
  jac(4*N1+p-1,ind)    =   1/(2*h)*[bxp(p-1,1) byp(p-1,1) bzp(p-1,1)];
  
    %-- differentiation with rho
   
   jac(ind,3*N1+p-1) = -[bx(p) ;by(p) ;bz(p)];
  
    %-- differentiation with lambda 
     jac(ind,4*N1+p-1:4*N1+p+1)   = [ -1/(2*h)*bxp(p)   bx2p(p) 1/(2*h)*bxp(p)  ; 
                                                         -1/(2*h)*byp(p)    by2p(p) 1/(2*h)*byp(p)  ;
                                                         -1/(2*h)*bzp(p)   bz2p(p)  1/(2*h)*bzp(p) ] ;
  %-- differentiation with u v w 
     jac(ind,5*N1+3:5*N1+5)       = -[   0          -bzp(p)    byp(p)  ;
                                                      bzp(p)       0        -bxp(p) ;
                                                    -byp(p)      bxp(p)      0    ] ;
  %---- closure ----
       
    
  jac(5*N1+3:5*N1+5,ind) =  [  0                  -(bz(p-1)-bz(p+1))  (by(p-1)-by(p+1))   ;
                               (bz(p-1)-bz(p+1))   0                 -(bx(p-1)-bx(p+1))   ;
                              -(by(p-1)-by(p+1))   (bx(p-1)-bx(p+1))  0                 ] ;  
  %======================================================================
                      
                        
     %======================================================================
    end
     
 %----- Euler--Lagrange and Jacobian at the boundary points 
  
  p = N-1;
  
        gmx = v*bzp(p,1) - w*byp(p,1);
        gmy = w*bxp(p,1) - u*bzp(p,1);
        gmz = u*byp(p,1) - v*bxp(p,1);
        
        f_b(3*p-5,1) = .001*bx4p(p-1,1) + lmp(p-1)*bxp(p,1) + lm(p)*bx2p(p,1) - rho(p-1)*bx(p,1) + gmx;
        f_b(3*p-4,1) = .001*by4p(p-1,1) + lmp(p-1)*byp(p,1) + lm(p)*by2p(p,1) - rho(p-1)*by(p,1) + gmy;
        f_b(3*p-3,1) = .001*bz4p(p-1,1) + lmp(p-1)*bzp(p,1) + lm(p)*bz2p(p,1) - rho(p-1)*bz(p,1) + gmz;
        
        
        f_b(3*N1+p-1,1) = 1/2*(bx(p,1)^2 + by(p,1)^2 + bz(p,1)^2 -1);
        
        f_b(4*N1+p,1) = 1/2*(bxp(p,1)^2 + byp(p,1)^2 + bzp(p,1)^2 - tau^2); 
        
       %--- closure -----
       
%         f_c  = f_c       + [by(p-1,1)*bz(p,1) - bz(p-1,1)*by(p,1) ;
%                             bz(p-1,1)*bx(p,1) - bx(p-1,1)*bz(p,1) ;
%                             bx(p-1,1)*by(p,1) - by(p-1,1)*bx(p,1) ] ; 
        f_c  = f_c     + [by(p,1)*bzp(p,1) - bz(p,1)*byp(p,1) ;
                              bz(p,1)*bxp(p,1) - bx(p,1)*bzp(p,1) ;
                              bx(p,1)*byp(p,1) - by(p,1)*bxp(p,1) ] ;                       
  %==================== Jacobian entry ================================
  
  ind = 3*p-5:3*p-3;
  
  jac(ind,ind-6)   = id*1*.001/h^4 ;
  jac(ind,ind-3)   = id*(-4*.001/h^4 - lmp(p-1)/(2*h) + lm(p)/h^2) -om/(2*h) ;
  jac(ind,ind)     = id*( 6*.001/h^4 - 2*lm(p)/h^2-rho(p-1))                 ; 
  jac(ind,ind+3)   = id*(-4*.001/h^4 + lmp(p-1)/(2*h) + lm(p)/h^2) +om/(2*h) ;                                            ;
   
   %--- Constraint gradient ---
  
  %--- unit vector constraint --
  jac(3*N1+p-1,ind) = [bx(p) by(p) bz(p)];
  
  %--- unit speed constraint --
  jac(4*N1+p-1,ind-6)    = -1/(2*h)*[bxp(p-1,1) byp(p-1,1) bzp(p-1,1)];
  jac(4*N1+p-1,ind)      =  1/(2*h)*[bxp(p-1,1) byp(p-1,1) bzp(p-1,1)];
  
      %-- differentiation with rho
   
     jac(ind,3*N1+p-1) = -[bx(p) ;by(p) ;bz(p)];
     %-- differentiation with lambda 
     jac(ind,4*N1+p-1:4*N1+p+1)   = [ -1/(2*h)*bxp(p)   bx2p(p) 1/(2*h)*bxp(p)  ; 
                                                         -1/(2*h)*byp(p)    by2p(p) 1/(2*h)*byp(p)  ;
                                                         -1/(2*h)*bzp(p)   bz2p(p)  1/(2*h)*bzp(p) ] ;

       %-- differentiation with u v w 
     jac(ind,5*N1+3:5*N1+5)       = -[   0          -bzp(p)    byp(p)  ;
                                                      bzp(p)       0        -bxp(p) ;
                                                    -byp(p)      bxp(p)      0    ] ;                                        
  %---- closure ----
       
    
  jac(5*N1+3:5*N1+5,ind) =  [  0                  -(bz(p-1)-bz(p+1))  (by(p-1)-by(p+1))   ;
                               (bz(p-1)-bz(p+1))   0                 -(bx(p-1)-bx(p+1))   ;
                              -(by(p-1)-by(p+1))   (bx(p-1)-bx(p+1))  0                 ] ;  
                          
    %--- CLOSURE CONSTRAINT WITH RESPECT TO bext(1:3) and bext(4:6);
    
%   jac(5*N1+3:5*N1+5,5*N1+6:5*N1+8) =  -0.5*[  0     -bz(1)  by(1)   ;
%                                              bz(1)   0     -bx(1)   ;
%                                             -by(1)   bx(1)   0      ] ;  
  jac(5*N1+3:5*N1+5,5*N1+9:5*N1+11) =   0.5*[  0     -bz(N+1)   by(N+1)   ;
                                                                         bz(N+1)   0     -bx(N+1)   ;
                                                                       -by(N+1)  bx(N+1)   0      ] ;  
  %======================================================================
   
  %======================================================================
   
 %----- Euler--Lagrange and Jacobian at the boundary points 
  
  p = N;
  
        gmx = v*bzp(p,1) - w*byp(p,1);
        gmy = w*bxp(p,1) - u*bzp(p,1);
        gmz = u*byp(p,1) - v*bxp(p,1);
        
        f_b(3*p-5,1) = .001*bx4p(p-1,1) + lmp(p-1)*bxp(p,1) + lm(p)*bx2p(p,1) - rho(p-1)*bx(p,1) + gmx;
        f_b(3*p-4,1) = .001*by4p(p-1,1) + lmp(p-1)*byp(p,1) + lm(p)*by2p(p,1) - rho(p-1)*by(p,1) + gmy;
        f_b(3*p-3,1) = .001*bz4p(p-1,1) + lmp(p-1)*bzp(p,1) + lm(p)*bz2p(p,1) - rho(p-1)*bz(p,1) + gmz;
        
        
        f_b(3*N1+p-1,1) = 1/2*(bx(p,1)^2 + by(p,1)^2 + bz(p,1)^2 -1);
        
        f_b(4*N1+p,1) = 1/2*(bxp(p,1)^2 + byp(p,1)^2 + bzp(p,1)^2 - tau^2); 
        
       %--- closure -----
       
%         f_c  = f_c       + [by(p-1,1)*bz(p,1) - bz(p-1,1)*by(p,1) ;
%                             bz(p-1,1)*bx(p,1) - bx(p-1,1)*bz(p,1) ;
%                             bx(p-1,1)*by(p,1) - by(p-1,1)*bx(p,1) ] ; 
       f_c  = f_c     + [by(p,1)*bzp(p,1) - bz(p,1)*byp(p,1) ;
                              bz(p,1)*bxp(p,1) - bx(p,1)*bzp(p,1) ;
                              bx(p,1)*byp(p,1) - by(p,1)*bxp(p,1) ] ;                        
  %==================== Jacobian entry ================================
  
  ind = 3*p-5:3*p-3;
  
  jac(ind,ind-6)   = id*1*.001/h^4 ;
  jac(ind,ind-3)   = id*(-4*.001/h^4 - lmp(p-1)/(2*h) + lm(p)/h^2) -om/(2*h) ;
  jac(ind,ind)     = id*( 6*.001/h^4 - 2*lm(p)/h^2-rho(p-1))                 ; 
  
  jac(ind,5*N1+9:5*N1+11)   = id*1*.001/h^4  ;       
   
   %-- differentiation with rho
   
   jac(ind,3*N1+p-1) = -[bx(p) ;by(p) ;bz(p)];
     %-- differentiation with lambda 
   jac(ind,4*N1+p-1:4*N1+p+1)   = [ -1/(2*h)*bxp(p)   bx2p(p) 1/(2*h)*bxp(p)  ; 
                                                         -1/(2*h)*byp(p)    by2p(p) 1/(2*h)*byp(p)  ;
                                                         -1/(2*h)*bzp(p)   bz2p(p)  1/(2*h)*bzp(p) ] ;
   %---- closure ----
       
    
  jac(5*N1+3:5*N1+5,ind) =  [  0                  -(bz(p-1)-bz(p+1))  (by(p-1)-by(p+1))   ;
                               (bz(p-1)-bz(p+1))   0                 -(bx(p-1)-bx(p+1))   ;
                              -(by(p-1)-by(p+1))   (bx(p-1)-bx(p+1))  0                 ] ;  

  %======================================================================  
    %--- Constraint gradient ---
  
  %--- unit vector constraint --
  jac(3*N1+p-1,ind) = [bx(p) by(p) bz(p)];
  %--- unit speed constraint --
  jac(4*N1+p-1,ind-6)    = -1/(2*h)*[bxp(p-1,1) byp(p-1,1) bzp(p-1,1)];
  jac(4*N1+p-1,ind)      =  1/(2*h)*[bxp(p-1,1) byp(p-1,1) bzp(p-1,1)];
  
  
   %-- differentiation with lambda 
     jac(ind,4*N1+p-1:4*N1+p+1)   = [ -1/(2*h)*bxp(p)   bx2p(p) 1/(2*h)*bxp(p)  ; 
                                                         -1/(2*h)*byp(p)    by2p(p) 1/(2*h)*byp(p)  ;
                                                         -1/(2*h)*bzp(p)   bz2p(p)  1/(2*h)*bzp(p) ] ;
                                                     
   %-- differentiation with u v w 
     jac(ind,5*N1+3:5*N1+5)       = -[   0          -bzp(p)    byp(p)  ;
                                                      bzp(p)       0        -bxp(p) ;
                                                    -byp(p)      bxp(p)      0    ] ;
  %-- unit speed constraint at i = N
  
  p = N+1;
  
  ind = 3*p-5:3*p-3;  % = 3*N1+1:3*N1+3
  
  jac(4*N1+p-1,ind-6)    = -1/(2*h)*[bxp(p-1,1) byp(p-1,1) bzp(p-1,1)];   % differentiation with b(N-1)
  
  
  %--unit speed at i = N+1 
  p = N+2 ;
  ind = 3*p-5:3*p-3;     % = 3*N1+4:3*N1+6
  % differentiation with bext(4:6)
  jac(4*N1+p-1,5*N1+9:5*N1+11)  =  1/(2*h)*[bxp(p-1,1) byp(p-1,1) bzp(p-1,1)];  % b'(N+1)   
  % differentiation with b(N) 
  jac(4*N1+p-1,ind-6)           = -1/(2*h)*[bxp(p-1,1) byp(p-1,1) bzp(p-1,1)];

                            
    
  %======================================================================
   
    
    
%     len(p) = sqrt(sum((bi-bip1).^2)) ; 
     %--- end points ----
   
   p = N1+2;
   
     
  f_b(5*N1+3:5*N1+5,1) =  h*f_c +   h*[by(p,1)*bzp(p,1)  - bz(p,1)*byp(p,1)   ;
                                                            bz(p,1)*bxp(p,1) - bx(p,1)*bzp(p,1)   ;
                                                            bx(p,1)*byp(p,1) - by(p,1)*bxp(p,1) ] ;  

                                                  
 %=============== Function evaluation =========================
        
   %--- |b'(0)| = tau
  f_b(4*N1+1,1) = 1/2*(bxp(1,1)^2 + byp(1,1)^2 + bzp(1,1)^2 - tau^2)   ;
  
  % |b'(1)| = tau
  f_b(5*N1+2,1) = 1/2*(bxp(N+1,1)^2 + byp(N+1,1)^2 + bzp(N+1,1)^2 - tau^2)   ;
  
  % ---- b'(1) + b'(N+1) = 0 ----
  
  f_b(5*N1+6,1) = bxp(1,1) + bxp(N+1,1);% bxp1^2  + byp1^2  + bzp1^2  - tau^2 ;
  f_b(5*N1+7,1) = byp(1,1) + byp(N+1,1);%bxpN1^2 + bypN1^2 + bzpN1^2 - tau^2 ;
  f_b(5*N1+8,1) = bzp(1,1) + bzp(N+1,1);
  
  %============= Jacobian evaluation 
  
  jac(5*N1+6:5*N1+8,1:3)           =  1/(2*h)*id;
  jac(5*N1+6:5*N1+8,5*N1+6:5*N1+8) = -1/(2*h)*id;
  
  jac(5*N1+6:5*N1+8,3*N1-2:3*N1)   = -1/(2*h)*id;
  jac(5*N1+6:5*N1+8,5*N1+9:5*N1+11) = 1/(2*h)*id;
 
                       
  %----- b''(1) + b''(N+1) = 0 ----
 
  
  f_b(5*N1+9,1)   = bx2p(1,1) + bx2p(N+1,1);
  f_b(5*N1+10,1)  = by2p(1,1) + by2p(N+1,1);
  f_b(5*N1+11,1)  = bz2p(1,1) + bz2p(N+1,1);
  
   %============= Jacobian evaluation 
  
  jac(5*N1+9:5*N1+11,1:3)      = 1/h^2*id;
  jac(5*N1+9:5*N1+11,5*N1+6:5*N1+8)  = 1/h^2*id;
  
  jac(5*N1+9:5*N1+11,3*N1-2:3*N1)    = 1/h^2*id;
  jac(5*N1+9:5*N1+11,5*N1+9:5*N1+11) = 1/h^2*id;
  
 % jac(:,3*N1+1:5*N1+11) = jac(3*N1+1:5*N1+11,:)';   
  
  %=================  Deleting columns and rows with respect to entries for
  %bz(2) and bz(N)
  J2 = jac;
  
  f_b(3,:)      = [];
  f_b(3*N1-1,:) = [];
  
  jac(3,:)      = [];
  jac(:,3)      = [];
  
  jac(3*N1-1,:) = [];
  jac(:,3*N1-1) = [];
  
  
 F = f_b;
  J = jac;
 err = sqrt(sum(f_b.^2)) ;  
 
     end
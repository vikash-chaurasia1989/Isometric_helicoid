function [F,J] = fun_pullback2(x)
 
%====  
%=== This code gives configuration on the constrained manifold nearest to
%given point 
 
 
global  N  N1  u v w       len id uni err qd bx2 E count bxN h tau rho lm bxg byg bzg sig bx0 by0 bz0
global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p bext J2 fac

 
  %x = initial_guess_symmetry();%initial_guess_evolution();
      
     count = count+1 ;
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
    
  
  %----- Euler--Lagrange and Jacobian at the boundary points 
 % fac = 0.001;
  
  %== For studying the effect of aspect ratio on the lagrange multipliers
  %===
  
  %fac = asinh(tau*pi*sig)/(4*pi^3*tau^3);
   
  
  p = 2;
  
           
        gm = om*[bx(p+1)-bx(p-1);by(p+1)-by(p-1);bz(p+1)-bz(p-1)];
        gmx = gm(1);
        gmy = gm(2);
        gmz = gm(3);
        
        
        F(3*p-5,1) = bx(p,1)-bx0(p,1) + lmp(p-1)*bxp(p,1) + lm(p)*bx2p(p,1) - rho(p-1)*bx(p,1) + gmx;
        F(3*p-4,1) = by(p,1)-by0(p,1) + lmp(p-1)*byp(p,1) + lm(p)*by2p(p,1) - rho(p-1)*by(p,1) + gmy;
        F(3*p-3,1) = bz(p,1)-bz0(p,1) + lmp(p-1)*bzp(p,1) + lm(p)*bz2p(p,1) - rho(p-1)*bz(p,1) + gmz;
        
        
        F(3*N1+p-1,1) = 1/2*(bx(p,1)^2 + by(p,1)^2 + bz(p,1)^2 -1);
        
        F(4*N1+p-1,1) = 1/2*(bxp(p,1)^2 + byp(p,1)^2 + bzp(p,1)^2 - tau^2); 
        
       %--- closure -----
%        
        f_c  = f_c       + [by(p-1,1)*bz(p,1) - bz(p-1,1)*by(p,1)    ;
                            bz(p-1,1)*bx(p,1) - bx(p-1,1)*bz(p,1)    ;
                            bx(p-1,1)*by(p,1) - by(p-1,1)*bx(p,1) ] ; 
%   
  %==================== Jacobian entry ================================
  
  ind = 3*p-5:3*p-3;
  
  J(ind,ind)     = id*(1 + lmp(p-1,1)/h -2*lm(p)/h^2-rho(p-1))    ;  
  J(ind,ind+3)   = id*(lm(p)/h^2)  +om                              ; 
    
 %   
  %--- Constraint gradient ---
  
  %--- unit vector constraint --
  J(3*N1+p-1,ind) = [bx(p) by(p) bz(p)];
  %--- unit speed constraint --
  J(4*N1+p-1,ind)  =  1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)];
   
   %---- closure ----
  J(5*N1+2:5*N1+4,ind) =   [ 0                 -(bz(p-1)-bz(p+1))   (by(p-1)-by(p+1))  ;
                               (bz(p-1)-bz(p+1))       0             -(bx(p-1)-bx(p+1))  ;
                              -(by(p-1)-by(p+1))  (bx(p-1)-bx(p+1))           0        ] ;  
 
   %-- differentiation with rho
   
   J(ind,3*N1+p-1) = -[bx(p) ;by(p) ;bz(p)];
  
   %-- differentiation with lambda 
   J(ind,4*N1+p-1)   = -1/h*[bxp(p) ;byp(p) ;bzp(p)]; 
   J(ind,4*N1+p)     =  1/h*[bxp(p) ;byp(p) ;bzp(p)] +  [bx2p(p) ;by2p(p) ;bz2p(p)]; 

  %-- differentiation with u v w 
   J(ind,5*N1+2:5*N1+4)       = - [ 0                 -(bz(p-1)-bz(p+1))   (by(p-1)-by(p+1))  ;
                                    (bz(p-1)-bz(p+1))       0             -(bx(p-1)-bx(p+1))  ;
                                   -(by(p-1)-by(p+1))  (bx(p-1)-bx(p+1))           0        ] ;  
    
                                          
    %======================================================================
  
  %----- Euler--Lagrange and Jacobian at the boundary points 
  
  p = 3;
  
          
        gm = om*[bx(p+1)-bx(p-1);by(p+1)-by(p-1);bz(p+1)-bz(p-1)];
        gmx = gm(1);
        gmy = gm(2);
        gmz = gm(3);
% %         
        
        F(3*p-5,1) = bx(p,1)-bx0(p,1) + lmp(p-1)*bxp(p,1) + lm(p)*bx2p(p,1) - rho(p-1)*bx(p,1) + gmx;
        F(3*p-4,1) = by(p,1)-by0(p,1) + lmp(p-1)*byp(p,1) + lm(p)*by2p(p,1) - rho(p-1)*by(p,1) + gmy;
        F(3*p-3,1) = bz(p,1)-bz0(p,1) + lmp(p-1)*bzp(p,1) + lm(p)*bz2p(p,1) - rho(p-1)*bz(p,1) + gmz;
        
        
        F(3*N1+p-1,1) = 1/2*(bx(p,1)^2 + by(p,1)^2 + bz(p,1)^2 -1);
        
        F(4*N1+p-1,1) = 1/2*(bxp(p,1)^2 + byp(p,1)^2 + bzp(p,1)^2 - tau^2); 
        
       %--- closure -----
       
        f_c  = f_c       + [by(p-1,1)*bz(p,1) - bz(p-1,1)*by(p,1) ;
                            bz(p-1,1)*bx(p,1) - bx(p-1,1)*bz(p,1) ;
                            bx(p-1,1)*by(p,1) - by(p-1,1)*bx(p,1) ] ; 
                    
  %==================== Jacobian entry ================================
  
  ind = 3*p-5:3*p-3;
  
  J(ind,ind-3)   = id*(  - lmp(p-1,1)/h + lm(p)/h^2) - om        ;
  J(ind,ind)     = id*(1 + lmp(p-1,1)/h -2*lm(p)/h^2-rho(p-1))        ;  
  J(ind,ind+3)   = id*(   +lm(p)/h^2)    +om  ; 
   
  
  
  %--- Constraint gradient ---
  
  %--- unit vector constraint --
  J(3*N1+p-1,ind) = [bx(p) by(p) bz(p)];
  %--- unit speed constraint --
  J(4*N1+p-1,ind-3) = -1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)];
  J(4*N1+p-1,ind)   =  1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)];
   
   %---- closure ----
  J(5*N1+2:5*N1+4,ind) =   [ 0                 -(bz(p-1)-bz(p+1))   (by(p-1)-by(p+1))  ;
                               (bz(p-1)-bz(p+1))       0             -(bx(p-1)-bx(p+1))  ;
                              -(by(p-1)-by(p+1))  (bx(p-1)-bx(p+1))           0        ] ;  
 
   %-- differentiation with rho
   
   J(ind,3*N1+p-1) = -[bx(p) ;by(p) ;bz(p)];
  
   %-- differentiation with lambda 
   J(ind,4*N1+p-1)   = -1/h*[bxp(p) ;byp(p) ;bzp(p)]; 
   J(ind,4*N1+p)     =  1/h*[bxp(p) ;byp(p) ;bzp(p)] +  [bx2p(p) ;by2p(p) ;bz2p(p)]; 

  %-- differentiation with u v w 
   J(ind,5*N1+2:5*N1+4)       = - [ 0                 -(bz(p-1)-bz(p+1))   (by(p-1)-by(p+1))  ;
                                    (bz(p-1)-bz(p+1))       0             -(bx(p-1)-bx(p+1))  ;
                                      -(by(p-1)-by(p+1))  (bx(p-1)-bx(p+1))           0        ] ;  
  
                                          
    %======================================================================
  
  
  
    for p = 4:N-2
          
        
    
        
        gm = om*[bx(p+1)-bx(p-1);by(p+1)-by(p-1);bz(p+1)-bz(p-1)];
        gmx = gm(1);
        gmy = gm(2);
        gmz = gm(3);
        
        
        F(3*p-5,1) = bx(p,1)-bx0(p,1) + lmp(p-1)*bxp(p,1) + lm(p)*bx2p(p,1) - rho(p-1)*bx(p,1) + gmx;
        F(3*p-4,1) = by(p,1)-by0(p,1) + lmp(p-1)*byp(p,1) + lm(p)*by2p(p,1) - rho(p-1)*by(p,1) + gmy;
        F(3*p-3,1) = bz(p,1)-bz0(p,1) + lmp(p-1)*bzp(p,1) + lm(p)*bz2p(p,1) - rho(p-1)*bz(p,1) + gmz;
        
        
        F(3*N1+p-1,1) = 1/2*(bx(p,1)^2 + by(p,1)^2 + bz(p,1)^2 -1);
        
        F(4*N1+p-1,1) = 1/2*(bxp(p,1)^2 + byp(p,1)^2 + bzp(p,1)^2 - tau^2); 
        
       %--- closure -----
       
          f_c  = f_c       + [by(p-1,1)*bz(p,1) - bz(p-1,1)*by(p,1) ;
                              bz(p-1,1)*bx(p,1) - bx(p-1,1)*bz(p,1) ;
                              bx(p-1,1)*by(p,1) - by(p-1,1)*bx(p,1) ] ;  
       
                       
  %==================== Jacobian entry ================================
  
  ind = 3*p-5:3*p-3;
  
  J(ind,ind-3)  = id*( - lmp(p-1,1)/h + lm(p)/h^2) - om            ;
  J(ind,ind)     = id*(1 + lmp(p-1,1)/h -2*lm(p)/h^2-rho(p-1))     ;  
  J(ind,ind+3)   = id*(   +lm(p)/h^2)    + om                      ; 
    
  
  %--- Constraint gradient ---
  
  %--- unit vector constraint --
  J(3*N1+p-1,ind) = [bx(p) by(p) bz(p)];
  %--- unit speed constraint --
  J(4*N1+p-1,ind-3) = -1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)];
  J(4*N1+p-1,ind)   =  1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)];
   
   %---- closure ----
  J(5*N1+2:5*N1+4,ind) =   [ 0                 -(bz(p-1)-bz(p+1))   (by(p-1)-by(p+1))  ;
                               (bz(p-1)-bz(p+1))       0             -(bx(p-1)-bx(p+1))  ;
                              -(by(p-1)-by(p+1))  (bx(p-1)-bx(p+1))           0        ] ;  
 
   %-- differentiation with rho
   
   J(ind,3*N1+p-1) = -[bx(p) ;by(p) ;bz(p)];
  
   %-- differentiation with lambda 
   J(ind,4*N1+p-1)   = -1/h*[bxp(p) ;byp(p) ;bzp(p)]; 
   J(ind,4*N1+p)     =  1/h*[bxp(p) ;byp(p) ;bzp(p)] +  [bx2p(p) ;by2p(p) ;bz2p(p)]; 

  %-- differentiation with u v w 
     J(ind,5*N1+2:5*N1+4)       = - [ 0                 -(bz(p-1)-bz(p+1))   (by(p-1)-by(p+1))  ;
                                    (bz(p-1)-bz(p+1))       0             -(bx(p-1)-bx(p+1))  ;
                                      -(by(p-1)-by(p+1))  (bx(p-1)-bx(p+1))           0        ] ;  

                                          
    %======================================================================
    end
     
 %----- Euler--Lagrange and Jacobian at the boundary points 
  
  p = N-1;
  
            
        gm = om*[bx(p+1)-bx(p-1);by(p+1)-by(p-1);bz(p+1)-bz(p-1)];
        gmx = gm(1);
        gmy = gm(2);
        gmz = gm(3);
        
        
        F(3*p-5,1) = bx(p,1)-bx0(p,1) + lmp(p-1)*bxp(p,1) + lm(p)*bx2p(p,1) - rho(p-1)*bx(p,1) + gmx;
        F(3*p-4,1) = by(p,1)-by0(p,1) + lmp(p-1)*byp(p,1) + lm(p)*by2p(p,1) - rho(p-1)*by(p,1) + gmy;
        F(3*p-3,1) = bz(p,1)-bz0(p,1) + lmp(p-1)*bzp(p,1) + lm(p)*bz2p(p,1) - rho(p-1)*bz(p,1) + gmz;
        
        
        F(3*N1+p-1,1) = 1/2*(bx(p,1)^2 + by(p,1)^2 + bz(p,1)^2 -1);
        
        F(4*N1+p-1,1) = 1/2*(bxp(p,1)^2 + byp(p,1)^2 + bzp(p,1)^2 - tau^2); 
        
       %--- closure -----
       
        f_c  = f_c       + [by(p-1,1)*bz(p,1) - bz(p-1,1)*by(p,1) ;
                            bz(p-1,1)*bx(p,1) - bx(p-1,1)*bz(p,1) ;
                            bx(p-1,1)*by(p,1) - by(p-1,1)*bx(p,1) ] ; 
                    
  %==================== Jacobian entry ================================
  
  ind = 3*p-5:3*p-3;
  
  J(ind,ind-3)   = id*(  - lmp(p-1,1)/h + lm(p)/h^2) - om            ;
  J(ind,ind)     = id*(1 + lmp(p-1,1)/h -2*lm(p)/h^2-rho(p-1))   ;  
  J(ind,ind+3)   = id*(0   +lm(p)/h^2)             +om                   ; 
   
  
  
  %--- Constraint gradient ---
  
  %--- unit vector constraint --
  J(3*N1+p-1,ind) = [bx(p) by(p) bz(p)];
  %--- unit speed constraint --
  J(4*N1+p-1,ind-3) = -1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)];
  J(4*N1+p-1,ind)   =  1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)];
   
   %---- closure ----
  J(5*N1+2:5*N1+4,ind) =   [ 0                 -(bz(p-1)-bz(p+1))   (by(p-1)-by(p+1))  ;
                               (bz(p-1)-bz(p+1))       0             -(bx(p-1)-bx(p+1))  ;
                              -(by(p-1)-by(p+1))  (bx(p-1)-bx(p+1))           0        ] ;  
 
   %-- differentiation with rho
   
   J(ind,3*N1+p-1) = -[bx(p) ;by(p) ;bz(p)];
  
   %-- differentiation with lambda 
   J(ind,4*N1+p-1)   = -1/h*[bxp(p) ;byp(p) ;bzp(p)]; 
   J(ind,4*N1+p)     =  1/h*[bxp(p) ;byp(p) ;bzp(p)] +  [bx2p(p) ;by2p(p) ;bz2p(p)]; 

   %-- differentiation with u v w 
   J(ind,5*N1+2:5*N1+4)       = - [ 0                 -(bz(p-1)-bz(p+1))   (by(p-1)-by(p+1))  ;
                                    (bz(p-1)-bz(p+1))       0             -(bx(p-1)-bx(p+1))  ;
                                      -(by(p-1)-by(p+1))  (bx(p-1)-bx(p+1))           0        ] ;  

                                          
    %======================================================================
   
 %----- Euler--Lagrange and Jacobian at the boundary points 
  
 p = N;
 
  
         gm = om*[bx(p+1)-bx(p-1);by(p+1)-by(p-1);bz(p+1)-bz(p-1)];
         gmx = gm(1);
         gmy = gm(2);
         gmz = gm(3);
 
 
 
 F(3*p-5,1) = bx(p,1)-bx0(p,1) + lmp(p-1)*bxp(p,1) + lm(p)*bx2p(p,1) - rho(p-1)*bx(p,1) + gmx;
 F(3*p-4,1) = by(p,1)-by0(p,1) + lmp(p-1)*byp(p,1) + lm(p)*by2p(p,1) - rho(p-1)*by(p,1) + gmy;
 F(3*p-3,1) = bz(p,1)-bz0(p,1) + lmp(p-1)*bzp(p,1) + lm(p)*bz2p(p,1) - rho(p-1)*bz(p,1) + gmz;
 
 
 F(3*N1+p-1,1) = 1/2*(bx(p,1)^2 + by(p,1)^2 + bz(p,1)^2 -1)         ;
 
 F(4*N1+p-1,1) = 1/2*(bxp(p,1)^2 + byp(p,1)^2 + bzp(p,1)^2 - tau^2) ;
 
 %--- closure -----
 
 f_c  = f_c       + [  by(p-1,1)*bz(p,1) - bz(p-1,1)*by(p,1)   ;
                       bz(p-1,1)*bx(p,1) - bx(p-1,1)*bz(p,1)   ;
                       bx(p-1,1)*by(p,1) - by(p-1,1)*bx(p,1) ] ;
 
 
 F(5*N1+2:5*N1+4,1) =    f_c +     [by(N,1)*bz(N+1,1) - bz(N,1)*by(N+1,1)   ;
                                                          bz(N,1)*bx(N+1,1) - bx(N,1)*bz(N+1,1)   ;
                                                          bx(N,1)*by(N+1,1) - by(N,1)*bx(N+1,1) ] ;
 
 
 F(5*N1+1,1) = 1/2*(bxp(N+1,1)^2 + byp(N+1,1)^2 + bzp(N+1,1)^2 - tau^2)   ;
 
 %==================== Jacobian entry ================================
 
 ind = 3*p-5:3*p-3;
 
 J(ind,ind-3)   = id*(0 - lmp(p-1,1)/h + lm(p)/h^2) - om           ;
 J(ind,ind)     = id*(1 + lmp(p-1,1)/h -2*lm(p)/h^2-rho(p-1))    ;
 
 
    
  
  %J(5*N1+2:5*N1+4,3*N1-2:3*N1) = .5*[0 0 0; 0 0 -1;0 1 0];
  %--- Constraint gradient ---
  
  %--- unit vector constraint --
  J(3*N1+p-1,ind) = [bx(p) by(p) bz(p)];
  %--- unit speed constraint --
  J(4*N1+p-1,ind-3) = -1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)];     % b'(N)
  J(4*N1+p-1,ind)   =  1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)];
  
  J(5*N1+1,ind)   = -1/h*[bxp(N+1,1) byp(N+1,1) bzp(N+1,1)];  %b'(N+1)
   
   %---- closure ----
  J(5*N1+2:5*N1+4,ind) =    [ 0                 -(bz(p-1)-bz(p+1))   (by(p-1)-by(p+1))  ;
                               (bz(p-1)-bz(p+1))       0             -(bx(p-1)-bx(p+1))  ;
                              -(by(p-1)-by(p+1))  (bx(p-1)-bx(p+1))           0        ] ;  
 
   %-- differentiation with rho
   
   J(ind,3*N1+p-1) = -[bx(p) ;by(p) ;bz(p)];
  
   %-- differentiation with lambda 
   J(ind,4*N1+p-1)   = -1/h*[bxp(p) ;byp(p) ;bzp(p)]; 
   J(ind,4*N1+p)     =  1/h*[bxp(p) ;byp(p) ;bzp(p)] +  [bx2p(p) ;by2p(p) ;bz2p(p)]; 

   %-- differentiation with u v w 
   J(ind,5*N1+2:5*N1+4)       = - [ 0                 -(bz(p-1)-bz(p+1))   (by(p-1)-by(p+1))  ;
                                    (bz(p-1)-bz(p+1))       0             -(bx(p-1)-bx(p+1))  ;
                                      -(by(p-1)-by(p+1))  (bx(p-1)-bx(p+1))           0        ] ;  

    J(:,5*N1+2:5*N1+4) = -J(:,5*N1+2:5*N1+4);
                                          
    %======================================================================
 
 err = sqrt(sum(F.^2)) ;  
 
   end
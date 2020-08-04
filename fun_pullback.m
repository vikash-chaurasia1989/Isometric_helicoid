      function [F,J] =   fun_pullback(var_in)
  
 %========================================================================
 %====== In this function file, Lagrangian is given by 
 %
 %  L = \sum 1/2*(b(i)-b0(i))^2 - \sum \lambda_i(unit speed) - \sum \rho_i unit vector 
 %      - gamma * closing 
 %   No smoothness of binormal is imposed at the closing point 
 %   This file works for both orientable and non orientable mid-line 
    
 %--- variables -----
 % bi  ... i = 2: N 
 % l 
 % \rho_i ... i = 2:N
 % \lambda_i ...i = 1:N
 % 
 %========================================================================
 %========================================================================
 
 
global  N  N1  u v w   id       count   h   rho lm  tau bx by bz err bx0 by0 bz0
   
  %     var_in = initial_guess2()     ;
  count = count+1                  ;
   
  bx(2:N,1)    = var_in(1:3:3*N1-2,1);
  by(2:N,1)    = var_in(2:3:3*N1-1,1);
  bz(2:N,1)    = var_in(3:3:3*N1,1);

  
  
  %--- lagrange multipliers --
  rho = var_in(3*N1+1:4*N1 ,1)  ;
  lm  = var_in(4*N1+1:5*N1+1,1) ;
  
  u = var_in(5*N1+2,1);
  v = var_in(5*N1+3,1);
  w = var_in(5*N1+4,1);
 
   
 
   %---------  ------------------------------
   om = [ 0  -w  v   ;
          w   0 -u   ;
         -v   u  0 ] ;
  id = [1 0 0;
          0 1 0;
          0 0 1];
 
  f_c = [0;0;0];
  
  f_b =zeros(5*N1+4,1);
    
   for p = 2:N1+1
    
       i = 3*p-5:3*p-3;
       
       bi   = [bx(p)  ;by(p)   ; bz(p)  ];
       bip1 = [bx(p+1);by(p+1) ; bz(p+1)];
       bim1 = [bx(p-1);by(p-1) ; bz(p-1)];
       
    %-------- Function evaluation -----------   
       %--- Euler-- Lagrange multipliers ----
        f_b(i,1) =  bi-[bx0(p,1);by0(p,1);bz0(p,1)]-2*(rho(p-1) + lm(p) + lm(p-1))*bi  + 2*lm(p-1)*bim1 +   2*lm(p)*bip1  + om*(bip1-bim1 ) ;
      
        %---- constraints ---------
       %--  Unity constraint -----
       
       f_b(3*N1+p-1,1) = bx(p)^2+by(p)^2 + bz(p)^2-1;
       
       %--- inextensibility -----
       f_b(4*N1+p-1,1) = sum((bi-bim1).*(bi-bim1)) -tau^2*h^2;
       
       %--- closure -----
       
       f_c  = f_c       + [by(p-1)*bz(p) - bz(p-1)*by(p) ;
                           bz(p-1)*bx(p) - bx(p-1)*bz(p) ;
                           bx(p-1)*by(p) - by(p-1)*bx(p) ] ;  
    %======================================================================
    
    %==================== Hessian calculation =============================
    
    
    W(i,i) = (1 -2*(rho(p-1) + lm(p) + lm(p-1)))*id  ; 
    
    %----- constraint gradient -----
    A(i,p-1) = 2*bi;
    A(i,N1+p-1) = 2*(bi-bim1);
    A(i,2*N1+2:2*N1+4) = [ 0                  -(bim1(3)-bip1(3))   (bim1(2)-bip1(2)) ;
                           (bim1(3)-bip1(3))   0                  -(bim1(1)-bip1(1)) ;
   
                           -(bim1(2)-bip1(2))   (bim1(1)-bip1(1))   0                ]' ;
    if(p>2)
        
  
                           
    W(i,i-3) = (2*lm(p-1))*id - om;
      
    A(i-3,N1+p-1) = -2*(bi-bim1);
    end
    if(p<N1+1)
       
                             
    W(i,i+3) = (2*lm(p))*id + om;
     
    end
    
    
     
   end
   %--- end points ----
   
   p = N1+2;
   
   %=============== Function evaluation =========================
   f_b(4*N1+p-1,1) = (bx(p)-bx(p-1))^2 + (by(p)-by(p-1))^2 + (bz(p)-bz(p-1))^2  -tau^2*h^2;
     
   
   f_b(5*N1+2:5*N1+4,1) = f_c + [by(p-1)*bz(p) - bz(p-1)*by(p) ;
                                 bz(p-1)*bx(p) - bx(p-1)*bz(p) ;
                                 bx(p-1)*by(p) - by(p-1)*bx(p) ] ;  
   
  %==================== Hessian evaluation ================================
  A(3*N1-2:3*N1,2*N1+1) = -2*([bx(p)-bx(p-1);by(p)-by(p-1);bz(p)-bz(p-1)]);
   
 
   
                                
  qd(1:3*N1,1:3*N1)          =    W ;
  qd(1:3*N1,3*N1+1:5*N1+4)   =   -A ;
  qd(3*N1+1:5*N1+4,1:3*N1)   =    A';
 %        
 
 
 
  
  
  F = f_b;
  J = qd;
  
 err = sqrt(sum(f_b.^2))  ;  
 
 
  end
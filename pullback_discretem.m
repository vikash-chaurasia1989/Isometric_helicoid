     function [F,J] =   fun_pullback4(x)
  
 %========================================================================
 %====== In this function file, Lagrangian is given by 
 %
 %  L = l - \sum \lambda_i(unit speed) - \sum \rho_i unit vector 
 %      - gamma * closing 
 
 %--- variables -----
 % bi  ... i = 2: N 
 % l 
 % \rho_i ... i = 2:N
 % \lambda_i ...i = 1:N
 % 
 %========================================================================
 %========================================================================
 
 
global  N  N1  lm  err1 l    count  u v w  bx by bz     A tau
   %     x = initial_guess()     ;
  count = count+1                  ;
  
  %--- variables ----- 
  bx(2:N,1)     = x(1:3:3*N1-2,1)    ;
  by(2:N,1)     = x(2:3:3*N1-1,1)    ;
  bz(2:N,1)     = x(3:3:3*N1,1)      ;
  
  
     
 
  f_c = [0;0;0];
  
    
   for p = 2:N1+1
    
     
       ind = 3*p-5:3*p-3;
        %---- constraints ---------
       %--  Unity constraint -----
       
       f_b(p-1,1) = bx(p,1)^2+by(p,1)^2 + bz(p,1)^2-1;
       
       %--- inextensibility -----
       f_b(N1+p-1,1)=  1/2*(bxp(p,1)^2 + byp(p,1)^2 + bzp(p,1)^2 - tau^2); % sum((bi-bim1).*(bi-bim1)) -h^2;
       
       %--- closure -----
       
       f_c  = f_c       + [by(p-1,1)*bz(p,1) - bz(p-1,1)*by(p,1) ;
                                bz(p-1,1)*bx(p,1) - bx(p-1,1)*bz(p,1) ;
                                bx(p-1,1)*by(p,1) - by(p-1,1)*bx(p,1) ] ;  
    %======================================================================
       %----- constraint gradient -----
 
                          
    A(i,p-1) = 2*bi;                            % unit constraint
    
    A(i,N1+p-1) = 2*(bi-bim1);          % unit speed constraint 
    
        
    A(i,2*N1+2:2*N1+4) =      [   0                           -(bim1(3)-bip1(3))    (bim1(2)-bip1(2)) ;
                                               (bim1(3)-bip1(3))      0                         -(bim1(1)-bip1(1)) ;
                                             -(bim1(2)-bip1(2))      (bim1(1)-bip1(1))    0                        ]' ;
    if(p>2)
        A(i-3,N1+p-1) = -2*(bi-bim1);
    end
  
    
   end
    
   %--- end points ----
   
   p = N1+2;
   
   %=============== Function evaluation =========================
     
   f_b(N1+p-1,1) = 1/2*(bxp(p,1)^2 + byp(p,1)^2 + bzp(p,1)^2 - tau^2);
     
   
   f_b(2*N1+2:2*N1+4,1) = f_c + [by(p-1,1)*bz(p,1) - bz(p-1,1)*by(p,1) ;
                                                  bz(p-1,1)*bx(p,1) - bx(p-1,1)*bz(p,1) ;
                                                  bx(p-1,1)*by(p,1) - by(p-1,1)*bx(p,1) ] ; 
     
  %==================== Hessian evaluation ================================
   A(3*N1-2:3*N1,2*N1+1) = -2*([bx(p,1)-bx(p-1,1);by(p,1)-by(p-1,1);bz(p,1)-bz(p-1,1)]);
    
 
 
     F = f_b;
     J = A'   ;
 
%    
   err1 = sqrt(sum(f_b.^2));
  end
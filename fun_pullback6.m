      function [F,J] =   fun_pullback6(var_in)
  
 %========================================================================
%==== This is the version that gives same results as if discrete
% We use fun_descent7 to determine the next point in gradient descent
% direction. Then we calculate the distance d between this point and the
% previous point. Then this fold finds minimium energy configuration which
% is at distance d from the previous step at distance d. 

% This code uses discretized  spatial formulation. fun_pullback3 is the
% disretized version of this code. 
 % 
 %========================================================================
 %========================================================================
 
 
global  N  N1  u v w   tau     id      count   h   rho lm l bx by bz err bxm1 bym1 bzm1 d 
   
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
  del = var_in(5*N1+5,1);
 
   h = tau/N;
  
 
   %---------  ------------------------------
   om = [ 0  -w  v   ;
          w   0 -u   ;
         -v   u  0 ] ;
  id = [1 0 0;
          0 1 0;
          0 0 1];
 
  f_c = [0;0;0];
  fac = 1     ;
  
  
  F = zeros(5*N1+5,1);
  J = zeros(5*N1+5,5*N1+5);
    
   for p = 2:N1+1
    
       i = 3*p-5:3*p-3;
       
       bi   = [bx(p)  ;by(p)   ; bz(p)  ];
       bip1 = [bx(p+1);by(p+1) ; bz(p+1)];
       bim1 = [bx(p-1);by(p-1) ; bz(p-1)];
       
    %-------- Function evaluation -----------   
       %--- Euler-- Lagrange multipliers ----
        F(i,1) =  fac*(2*bi - (bip1+bim1) )-2*(rho(p-1) + lm(p) + lm(p-1))*bi  + 2*lm(p-1)*bim1 +   2*lm(p)*bip1  + om*(bip1-bim1 )...
                        +del*[bx(p,1)-bxm1(p,1);by(p,1)-bym1(p,1);bz(p,1)-bzm1(p,1)];
      % F(i,1) =  fac*(bi - bip1) -2*(rho(p-1) + lm(p) + lm(p-1))*bi  + 2*lm(p-1)*bim1 +   2*lm(p)*bip1  + om*(bip1-bim1 ) ;

        %---- constraints ---------
       %--  Unity constraint -----
       
       F(3*N1+p-1,1) = bx(p)^2+by(p)^2 + bz(p)^2-1;
       
       %--- inextensibility -----
       F(4*N1+p-1,1) = sum((bi-bim1).*(bi-bim1)) -h^2;
       
       %--- closure -----
       
       f_c  = f_c       + [by(p-1)*bz(p) - bz(p-1)*by(p) ;
                           bz(p-1)*bx(p) - bx(p-1)*bz(p) ;
                           bx(p-1)*by(p) - by(p-1)*bx(p) ] ;  
    %======================================================================
    
    %==================== Hessian calculation =============================
    
    
    W(i,i) = (2*fac -2*(rho(p-1) + lm(p) + lm(p-1))+del)*id  ; 
    
    %----- constraint gradient -----
    A(i,p-1) = 2*bi;
    A(i,N1+p-1) = 2*(bi-bim1);
    A(i,2*N1+2:2*N1+4) = [ 0                  -(bim1(3)-bip1(3))   (bim1(2)-bip1(2)) ;
                           (bim1(3)-bip1(3))   0                  -(bim1(1)-bip1(1)) ;
   
                           -(bim1(2)-bip1(2))   (bim1(1)-bip1(1))   0                ]' ;
    if(p>2)
        
    F(i,1) = F(i,1) + fac*1/2*[bx(p,1)-2*bx(p-1,1)+bx(p-2,1);
                                 by(p,1)-2*by(p-1,1)+by(p-2,1);
                                 bz(p,1)-2*bz(p-1,1)+bz(p-2,1)];  
                           
    W(i,i-3) = (-2*fac + 2*lm(p-1))*id - om;
    W(i,i)   = W(i,i) +fac*1/2*id;
    
    A(i-3,N1+p-1) = -2*(bi-bim1);
    end
    if(p<N1+1)
        F(i,1) = F(i,1) + fac*1/2*[bx(p,1)-2*bx(p+1,1)+bx(p+2,1);
                                     by(p,1)-2*by(p+1,1)+by(p+2,1);
                                     bz(p,1)-2*bz(p+1,1)+bz(p+2,1)];  
                             
    W(i,i+3) = (-2*fac+2*lm(p))*id + om;
    W(i,i)   = W(i,i) +fac*1/2*id;
   
    
    
    end
    
    if(p>3)
    W(i,i-6) = fac*1/2*id;
    end
    
    if(p<N1)
    W(i,i+6) = fac*1/2*id;    
    end
    
    %=== terms for distance from the curve bxm1,bym1,bzm1 ===
    J(i,5*N1+5) = [bx(p,1)-bxm1(p,1);by(p,1)-bym1(p,1);bz(p,1)-bzm1(p,1)];
    J(5*N1+5,i) = J(i,5*N1+5)';
   
   end
   %--- end points ----
   
   p = N1+2;
   
   %=============== Function evaluation =========================
   F(4*N1+p-1,1) = (bx(p)-bx(p-1))^2 + (by(p)-by(p-1))^2 + (bz(p)-bz(p-1))^2  -h^2;
     
   
   F(5*N1+2:5*N1+4,1) = f_c + [by(p-1)*bz(p) - bz(p-1)*by(p) ;
                                 bz(p-1)*bx(p) - bx(p-1)*bz(p) ;
                                 bx(p-1)*by(p) - by(p-1)*bx(p) ] ;  
   
  %==================== Hessian evaluation ================================
  A(3*N1-2:3*N1,2*N1+1) = -2*([bx(p)-bx(p-1);by(p)-by(p-1);bz(p)-bz(p-1)]);
   
   F(5*N1+5,1) = 1/2*(sum((bx-bxm1).^2 + (by-bym1).^2 + (bz-bzm1).^2) - d^2);

  
  %====== including b''_1^2 contribution 
% %   
%   F(1:3,1) =  F(1:3,1) + fac/2*[bx(2,1)-bx(N,1)-2*bx(1,1);
%                                     by(2,1)-by(N,1)-2*by(1,1);
%                                     bz(2,1)-bz(N,1)-2*bz(1,1)] ;
%                                 
%   W(1:3,1:3) =     W(1:3,1:3) + fac/2*id;
%   
%   
%   
%   F(3*N1-2:3*N1,1) =  F(3*N1-2:3*N1,1) - fac/2*[bx(2,1)-bx(N,1)-2*bx(1,1);
%                                                     by(2,1)-by(N,1)-2*by(1,1);
%                                                     bz(2,1)-bz(N,1)-2*bz(1,1)] ;
%                                 
%    W(3*N1-2:3*N1,3*N1-2:3*N1) =     W(3*N1-2:3*N1,3*N1-2:3*N1) + fac/2*id;                              
%                                 
   
   
%   F(1:3,1) =  F(1:3,1) - fac*[bx(1,1)-2*bx(2,1)+bx(3,1);
%                                     by(1,1)-2*by(2,1)+by(3,1);
%                                     bz(1,1)-2*bz(2,1)+bz(3,1);] ;
%                                 
%   W(1:3,1:3) =     W(1:3,1:3) +2*fac*id;
%   
%   
%   
%    F(4:6,1) =  F(4:6,1) + fac/2*[bx(1,1)-2*bx(2,1)+bx(3,1);
%                                     by(1,1)-2*by(2,1)+by(3,1);
%                                     bz(1,1)-2*bz(2,1)+bz(3,1);] ;
%                                 
%   W(4:6,4:6) =     W(4:6,4:6) +fac/2*id;                             
%                                 
   
                                
  J(1:3*N1,1:3*N1)          =    W ;
  J(1:3*N1,3*N1+1:5*N1+4) =   -A ;
  J(3*N1+1:5*N1+4,1:3*N1)   =    A';
 %        
 
 
 
 
  
 err = sqrt(sum(F.^2)) ;  
 
 
  end
  function  F  = fun_curve(x)
 

 %--- in this file, we solve for minimizing bending energy of the ruled surface --
   
   %  b'''' + lm'*b' + lm*b'' - rho b + gm x b' = 0;
   %
   % --- constraints --
   %   |b| = 1 and |b|' = \tau
   
   
global  N  N1  u v w    f_b   uni err  count   h tau rho lm 
global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p bext

 
%--- debugging lines ----

%     x = initial_guess2();
%    global x1
%    x1 = x;
%    rho1 = rho;
%    lm1  = lm;
%  
        
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
 
 
%             
%             rho = (rho(1:end-1,1)+rho(2:end,1))/2;
%             rho(N-1,1) = rho(1,1);
%             
%             lm = (lm(1:end-1,1)+lm(2:end,1))/2;
%             lm(N+1,1) = lm(1,1);
   % interpolation points 
%    
%           bx1 = x(5*N1+6,1) ;
%           by1 = x(5*N1+7,1) ;
%           bz1 = x(5*N1+8,1) ;
%           
%           bxN2= x(5*N1+9,1)  ;
%           byN2= x(5*N1+10,1) ;
%           bzN2= x(5*N1+11,1) ;
          
 %---- derivatives -----
 %  b' --- O(h^2) ---
       
     %-- extrapolated points 
     
   
 %---- derivative calculation for b2 --- bN  ---
 
   ig = 2:N;
   
   
   %--  b'   (O(h^2))----
   %-- boundary points ---
   %- forward Euler --
     bxp(1,1) = (bx(2,1)-bext(1,1))/(2*h);
     byp(1,1) = (by(2,1)-bext(2,1))/(2*h);
     bzp(1,1) = (bz(2,1)-bext(3,1))/(2*h);
   
%      bxp(1,1) = (bx(2,1)-bx(1,1))/(h);
%      byp(1,1) = (by(2,1)-by(1,1))/(h);
%      bzp(1,1) = (bz(2,1)-bz(1,1))/(h);
   
   
%    bxp(1,1) = 1/h*(-25/12*bxg(2,1) + 4*bxg(3,1) -3*bxg(4,1)+4/3*bxg(5,1)-1/4*bxg(6,1));
%    byp(1,1) = 1/h*(-25/12*byg(2,1) + 4*byg(3,1) -3*byg(4,1)+4/3*byg(5,1)-1/4*byg(6,1));
%    bzp(1,1) = 1/h*(-25/12*bzg(2,1) + 4*bzg(3,1) -3*bzg(4,1)+4/3*bzg(5,1)-1/4*bzg(6,1));
%    %-- Central difference
   %---O(h^2)---
     bxp(ig,1) = (bx(ig+1,1)-bx(ig-1,1))/(2*h);
     byp(ig,1) = (by(ig+1,1)-by(ig-1,1))/(2*h);
     bzp(ig,1) = (bz(ig+1,1)-bz(ig-1,1))/(2*h);
     
      %--- Backward Euler --
    bxp(N+1,1) = (bext(4,1)-bx(N,1))/(2*h);
    byp(N+1,1) = (bext(5,1)-by(N,1))/(2*h);
    bzp(N+1,1) = (bext(6,1)-bz(N,1))/(2*h);
      
%     bxp(N+1,1) = (bext(4,1)-bx(N+1,1))/(h);
%     byp(N+1,1) = (bext(5,1)-by(N+1,1))/(h);
%     bzp(N+1,1) = (bext(6,1)-bz(N+1,1))/(h);
   
%    %-- O(h^4) --
%        bxp(1,1) = 1/h*(-25/12*bx(1,1) + 4*bx(2,1)-3*bx(3,1)+4/3*bx(4,1)-1/4*bx(5,1));
%        byp(1,1) = 1/h*(-25/12*by(1,1) + 4*by(2,1)-3*by(3,1)+4/3*by(4,1)-1/4*by(5,1));
%        bzp(1,1) = 1/h*(-25/12*bz(1,1) + 4*bz(2,1)-3*bz(3,1)+4/3*bz(4,1)-1/4*bz(5,1));
%        
%        bxp(2,1) = 1/h*(1/12*bext(1,1)-2/3*bx(1,1) + 2/3*bx(3,1)-1/12*bx(4,1));
%        byp(2,1) = 1/h*(1/12*bext(2,1)-2/3*by(1,1) + 2/3*by(3,1)-1/12*by(4,1));
%        bzp(2,1) = 1/h*(1/12*bext(3,1)-2/3*bz(1,1) + 2/3*bz(3,1)-1/12*bz(4,1));
%      
%        ig = 3:N-1;
%                 
%        bxp(3:N-1,1) = (-bx(ig+2,1)+8*bx(ig+1,1) -8*bx(ig-1,1)+bx(ig-2,1))/(12*h);
%        byp(3:N-1,1) = (-by(ig+2,1)+8*by(ig+1,1) -8*by(ig-1,1)+by(ig-2,1))/(12*h);
%        bzp(3:N-1,1) = (-bz(ig+2,1)+8*bz(ig+1,1) -8*bz(ig-1,1)+bz(ig-2,1))/(12*h);
%       
%        bxp(N,1) = 1/h*(1/12*bx(N-2,1)-2/3*bx(N-1,1) + 2/3*bx(N+1,1)-1/12*bext(4,1));
%        byp(N,1) = 1/h*(1/12*by(N-2,1)-2/3*by(N-1,1) + 2/3*by(N+1,1)-1/12*bext(5,1));
%        bzp(N,1) = 1/h*(1/12*bz(N-2,1)-2/3*bz(N-1,1) + 2/3*bz(N+1,1)-1/12*bext(6,1));
   
           
%    
      ig = 2:N;
   
      
      
      
   
    lmp = (lm(ig+1,1)-lm(ig-1))/(2*h);
   
   %-- O(h^4) --
  %  lmp(1,1) = (lm(3,1)-lm(1,1))/(2*h);
   % lmp(1,1) = 1/h*(-25/12*lm(2,1) + 4*lm(3,1) -3*lm(4,1)+4/3*lm(5,1)-1/4*lm(6,1));
      
   % lmp(2:N-2,1) = (-lm(5:N+1,1)+8*lm(4:N,1) -8*lm(2:N-2,1)+lm(1:N-3,1))/(12*h);
     
   % lmp(N-1,1)  = (lm(N+1,1)-lm(N-1,1))/(2*h);
    
    
   %--  b''  (O(h^2))----
   
   bx2p(1,1) = (bx(2,1) + bext(1,1) -2*bx(1,1))/h^2 ;
   by2p(1,1) = (by(2,1) + bext(2,1) -2*by(1,1))/h^2 ;
   bz2p(1,1) = (bz(2,1) + bext(3,1) -2*bz(1,1))/h^2 ;
   
   %--- forward difference at the first point 
%    % O(h^2) 
  %    bx2p(1,1) = (2*bx(1,1)-5*bx(2,1)+4*bx(3,1)-bx(4,1))/h^2;%;+ bxg(1,1) -2*bxg(2,1))/h^2 ;
   %   by2p(1,1) = (2*by(1,1)-5*by(2,1)+4*by(3,1)-by(4,1))/h^2;
 %   bz2p(1,1) = (2*bz(1,1)-5*bz(2,1)+4*bz(3,1)-bz(4,1))/h^2;
   
     
   bx2p(ig,1) = (bx(ig+1,1) + bx(ig-1,1) - 2*bx(ig,1))/h^2 ;
   by2p(ig,1) = (by(ig+1,1) + by(ig-1,1) - 2*by(ig,1))/h^2 ;
   bz2p(ig,1) = (bz(ig+1,1) + bz(ig-1,1) - 2*bz(ig,1))/h^2 ;
   
   bx2p(N+1,1) = (bext(4,1) + bx(N,1)   - 2*bx(N+1,1))/h^2 ;
   by2p(N+1,1) = (bext(5,1) + by(N,1)   - 2*by(N+1,1))/h^2 ;
   bz2p(N+1,1) = (bext(6,1) + bz(N,1)   - 2*bz(N+1,1))/h^2 ;
%    
   %-- backward difference at the end point 
   % O(h^2) 
%      bx2p(N+1,1) = (2*bx(N+1,1)  -5*bx(N,1) + 4*bx(N-1,1) -bx(N-2,1))/h^2;%bxg(N+1,1) -2*bxg(N+2,1))/h^2 ;
%      by2p(N+1,1) = (2*by(N+1,1)  -5*by(N,1) + 4*by(N-1,1) -by(N-2,1))/h^2;
%      bz2p(N+1,1) = (2*bz(N+1,1)  -5*bz(N,1) + 4*bz(N-1,1) -bz(N-2,1))/h^2;
   %--- b'' O(h^4) 
   
%    bx2p(2:N,1) = (-bxg(ig+2,1) + 16*bxg(ig+1,1) -30*bxg(ig,1) + 16*bxg(ig-1,1) -bxg(ig-2,1))/(12*h^2);
%      by2p(2:N,1) = (-byg(ig+2,1) + 16*byg(ig+1,1) -30*byg(ig,1) + 16*byg(ig-1,1) -byg(ig-2,1))/(12*h^2);
%     bz2p(2:N,1) = (-bzg(ig+2,1) + 16*bzg(ig+1,1) -30*bzg(ig,1) + 16*bzg(ig-1,1) -bzg(ig-2,1))/(12*h^2);
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
   
   
  %--- data smoothening ---
%   
%   bx4p = (bx4p(1:end-1,1)+bx4p(2:end,1))/2;
%   bx4p(N-1,1) = -bx4p(1,1);
%    
%   by4p = (by4p(1:end-1,1)+by4p(2:end,1))/2;
%   by4p(N-1,1) = -by4p(1,1);
%     
%   bz4p = (bz4p(1:end-1,1)+bz4p(2:end,1))/2;
%   bz4p(N-1,1) = -bz4p(1,1);
%             
            
   %---------------------------------------
    
  f_c = [0;0;0];
    
  
    for p = 2:N1+1
          
        
        
        gmx = v*bzp(p,1) - w*byp(p,1);
        gmy = w*bxp(p,1) - u*bzp(p,1);
        gmz = u*byp(p,1) - v*bxp(p,1);
        
        f_b(3*p-5,1) = .001*bx4p(p-1,1) + lmp(p-1)*bxp(p,1) + lm(p)*bx2p(p,1) - rho(p-1)*bx(p,1) + gmx;
        f_b(3*p-4,1) = .001*by4p(p-1,1) + lmp(p-1)*byp(p,1) + lm(p)*by2p(p,1) - rho(p-1)*by(p,1) + gmy;
        f_b(3*p-3,1) = .001*bz4p(p-1,1) + lmp(p-1)*bzp(p,1) + lm(p)*bz2p(p,1) - rho(p-1)*bz(p,1) + gmz;
        
        
        f_b(3*N1+p-1,1) = 1/2*(bx(p,1)^2 + by(p,1)^2 + bz(p,1)^2 -1);
        
        %-- This loop will take values from b(2)' --- b(N)' ---
        f_b(4*N1+p,1) = 1/2*(bxp(p,1)^2 + byp(p,1)^2 + bzp(p,1)^2 - tau^2); 
        
       %--- closure -----
%        
%         f_c  = f_c       + [by(p-1,1)*bz(p,1) - bz(p-1,1)*by(p,1) ;
%                             bz(p-1,1)*bx(p,1) - bx(p-1,1)*bz(p,1) ;
%                             bx(p-1,1)*by(p,1) - by(p-1,1)*bx(p,1) ] ;  
                       
       f_c  = f_c     + [by(p,1)*bzp(p,1) - bz(p,1)*byp(p,1) ;
                              bz(p,1)*bxp(p,1) - bx(p,1)*bzp(p,1) ;
                              bx(p,1)*byp(p,1) - by(p,1)*bxp(p,1) ] ;
%                        
    %======================================================================
    
    
   end
%     len(p) = sqrt(sum((bi-bip1).^2)) ; 
      uni(p)=1/2*(bx(p+1)^2+by(p+1)^2 + bz(p+1)^2-1);
   %--- end points ----
   
   p = N1+2;
   
  %--- |b'(0)| = tau
  f_b(4*N1+1,1) = 1/2*(bxp(1,1)^2 + byp(1,1)^2 + bzp(1,1)^2 - tau^2)   ;
  
  % |b'(1)| = tau
  f_b(5*N1+2,1) = 1/2*(bxp(N+1,1)^2 + byp(N+1,1)^2 + bzp(N+1,1)^2 - tau^2)   ;

     
%  f_b(5*N1+3:5*N1+5,1)  = f_c    +   [by(p-1,1)*bz(p,1) - bz(p-1,1)*by(p,1)   ;
%                                      bz(p-1,1)*bx(p,1) - bx(p-1,1)*bz(p,1)   ;
%                                      bx(p-1,1)*by(p,1) - by(p-1,1)*bx(p,1) ] ;  
 %=============== Function evaluation =========================
        
   
   f_b(5*N1+3:5*N1+5,1) =  h*f_c +   h*[by(p,1)*bzp(p,1)  - bz(p,1)*byp(p,1)   ;
                                                            bz(p,1)*bxp(p,1) - bx(p,1)*bzp(p,1)   ;
                                                            bx(p,1)*byp(p,1) - by(p,1)*bxp(p,1) ] ;  
%    
  %--- boundary terms ---
  
  % |b(0)|=1, |b(N+2)| = 1
  
 % f_b(5*N1+4,1) = bxg(1,1)^2   + byg(1,1)^2   + byg(1,1)^2 - 1   ;
 % f_b(5*N1+5,1) = bxg(N+3,1)^2 + byg(N+3,1)^2 + byg(N+3,1)^2-1   ;
  
  
  
  % ---- b'(1) + b'(N+1) = 0 ----
  
  f_b(5*N1+6,1) = bxp(1,1) + bxp(N+1,1);% bxp1^2  + byp1^2  + bzp1^2  - tau^2 ;
  f_b(5*N1+7,1) = byp(1,1) + byp(N+1,1);%bxpN1^2 + bypN1^2 + bzpN1^2 - tau^2 ;
  f_b(5*N1+8,1) = bzp(1,1) + bzp(N+1,1);
  
  %----- b''(1) + b''(N+1) = 0 ----
 
  
  f_b(5*N1+9,1)   = bx2p(1,1) + bx2p(N+1,1);
  f_b(5*N1+10,1)  = by2p(1,1) + by2p(N+1,1);
  f_b(5*N1+11,1)  = bz2p(1,1) + bz2p(N+1,1);
  
%   
%   f_b(5*N1+12,1) = bxp(2,1) + bxp(N,1);% bxp1^2  + byp1^2  + bzp1^2  - tau^2 ;
%   f_b(5*N1+13,1) = byp(2,1) + byp(N,1);%bxpN1^2 + bypN1^2 + bzpN1^2 - tau^2 ;
%   f_b(5*N1+14,1) = bzp(2,1) + bzp(N,1);
  
  f_b(3,:) = [];
  f_b(3*N1-1,:) = [];
  
 F = f_b;
 err = sqrt(sum(f_b.^2)) ;  
 
   end
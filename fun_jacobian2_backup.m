  function F = fun_jacobian2(x)
 
 
 
global N a r1x r1y r1z r2x r2y r2z r3x r3y r3z lm12 lm13 lm23 rho gmx gmy gmz 
 

 %x = initial_guess2();
    %----- reading input and constructing bx by bz ----
%== position vectors ==


  r1x(1:N,1)     = x(1:9:9*N-8,1)     ;
  r1y(1:N,1)     = x(2:9:9*N-7,1)     ;
  r1z(1:N,1)     = x(3:9:9*N-6,1)     ;
  
  r2x(1:N,1)     = x(4:9:9*N-5,1)     ;
  r2y(1:N,1)     = x(5:9:9*N-4,1)     ;
  r2z(1:N,1)     = x(6:9:9*N-3,1)     ;
  
  r3x(1:N,1)     = x(7:9:9*N-2,1)     ;
  r3y(1:N,1)     = x(8:9:9*N-1,1)     ;
  r3z(1:N,1)     = x(9:9:9*N,1)       ;
  
  %=== Lagrange multipliers ===
  
  rho(1:N,1)     = x(9*N+1:10*N,1);
  
  lm12(1:N,1)    = x(10*N+1:11*N,1);
  lm13(1:N,1)    = x(11*N+1:12*N,1);
  lm23(1:N,1)    = x(12*N+1:13*N,1);
  
  gmx(1:N,1)     = x(13*N+1:14*N,1);
  gmy(1:N,1)     = x(14*N+1:15*N,1);
  gmz(1:N,1)     = x(15*N+1:16*N,1);
  
  
  F = zeros(16*N,1);
  J = zeros(16*N,16*N);
 
  
  for i = 1:N
      
  %=== Euler--Lagrange equations ===
  ind = 9*i-8:9*i;
  
  j = 1:N;
  j = j(j~=i);    % index from 1 to N except i
  
  %--- 1/r^3
  d11 = ((r1x(i,1)-r1x(j,1)).^2 + (r1y(i,1)-r1y(j,1)).^2 + (r1z(i,1)-r1z(j,1)).^2).^1.5;  
  d12 = ((r1x(i,1)-r2x(j,1)).^2 + (r1y(i,1)-r2y(j,1)).^2 + (r1z(i,1)-r2z(j,1)).^2).^1.5; 
  d13 = ((r1x(i,1)-r3x(j,1)).^2 + (r1y(i,1)-r3y(j,1)).^2 + (r1z(i,1)-r3z(j,1)).^2).^1.5; 
    
  d21 = ((r2x(i,1)-r1x(j,1)).^2 + (r2y(i,1)-r1y(j,1)).^2 + (r2z(i,1)-r1z(j,1)).^2).^1.5;  
  d22 = ((r2x(i,1)-r2x(j,1)).^2 + (r2y(i,1)-r2y(j,1)).^2 + (r2z(i,1)-r2z(j,1)).^2).^1.5; 
  d23 = ((r2x(i,1)-r3x(j,1)).^2 + (r2y(i,1)-r3y(j,1)).^2 + (r2z(i,1)-r3z(j,1)).^2).^1.5; 
  
  d31 = ((r3x(i,1)-r1x(j,1)).^2 + (r3y(i,1)-r1y(j,1)).^2 + (r3z(i,1)-r1z(j,1)).^2).^1.5;  
  d32 = ((r3x(i,1)-r2x(j,1)).^2 + (r3y(i,1)-r2y(j,1)).^2 + (r3z(i,1)-r2z(j,1)).^2).^1.5; 
  d33 = ((r3x(i,1)-r3x(j,1)).^2 + (r3y(i,1)-r3y(j,1)).^2 + (r3z(i,1)-r3z(j,1)).^2).^1.5; 
 
    %--- 1/r^5
  d11_5 = ((r1x(i,1)-r1x(j,1)).^2 + (r1y(i,1)-r1y(j,1)).^2 + (r1z(i,1)-r1z(j,1)).^2).^2.5;  
  d12_5 = ((r1x(i,1)-r2x(j,1)).^2 + (r1y(i,1)-r2y(j,1)).^2 + (r1z(i,1)-r2z(j,1)).^2).^2.5; 
  d13_5 = ((r1x(i,1)-r3x(j,1)).^2 + (r1y(i,1)-r3y(j,1)).^2 + (r1z(i,1)-r3z(j,1)).^2).^2.5; 
    
  d21_5 = ((r2x(i,1)-r1x(j,1)).^2 + (r2y(i,1)-r1y(j,1)).^2 + (r2z(i,1)-r1z(j,1)).^2).^2.5;  
  d22_5 = ((r2x(i,1)-r2x(j,1)).^2 + (r2y(i,1)-r2y(j,1)).^2 + (r2z(i,1)-r2z(j,1)).^2).^2.5; 
  d23_5 = ((r2x(i,1)-r3x(j,1)).^2 + (r2y(i,1)-r3y(j,1)).^2 + (r2z(i,1)-r3z(j,1)).^2).^2.5; 
  
  d31_5 = ((r3x(i,1)-r1x(j,1)).^2 + (r3y(i,1)-r1y(j,1)).^2 + (r3z(i,1)-r1z(j,1)).^2).^2.5;  
  d32_5 = ((r3x(i,1)-r2x(j,1)).^2 + (r3y(i,1)-r2y(j,1)).^2 + (r3z(i,1)-r2z(j,1)).^2).^2.5; 
  d33_5 = ((r3x(i,1)-r3x(j,1)).^2 + (r3y(i,1)-r3y(j,1)).^2 + (r3z(i,1)-r3z(j,1)).^2).^2.5; 
  
  
  %== r1 x r2 x r3 =====  
  r123x = (r1y(i,1)*r2z(i,1) - r1z(i,1)*r2y(i,1)) + (r2y(i,1)*r3z(i,1) - r2z(i,1)*r3y(i,1)) + (r3y(i,1)*r1z(i,1) - r3z(i,1)*r1y(i,1));
  r123y = (r1z(i,1)*r2x(i,1) - r1x(i,1)*r2z(i,1)) + (r2z(i,1)*r3x(i,1) - r2x(i,1)*r3z(i,1)) + (r3z(i,1)*r1x(i,1) - r3x(i,1)*r1z(i,1));
  r123z = (r1x(i,1)*r2y(i,1) - r1y(i,1)*r2x(i,1)) + (r2x(i,1)*r3y(i,1) - r2y(i,1)*r3x(i,1)) + (r3x(i,1)*r1y(i,1) - r3y(i,1)*r1x(i,1));

 
    
  r123c = [ 0     -r123z  r123y  ;  
            r123z  0     -r123x  ;
           -r123y  r123x  0     ];
  
  cx = r1x(i,1) + r2x(i,1) + r3x(i,1);
  cy = r1y(i,1) + r2y(i,1) + r3y(i,1);
  cz = r1z(i,1) + r2z(i,1) + r3z(i,1);
  
  
  r23c = (r2x(i,1)-r3x(i,1))*cx +  (r2y(i,1)-r3y(i,1))*cy + (r2z(i,1)-r3z(i,1))*cz;
  r31c = (r3x(i,1)-r1x(i,1))*cx +  (r3y(i,1)-r1y(i,1))*cy + (r3z(i,1)-r1z(i,1))*cz;
  r12c = (r1x(i,1)-r2x(i,1))*cx +  (r1y(i,1)-r2y(i,1))*cy + (r1z(i,1)-r2z(i,1))*cz;

  T1 = (r123c -[(r2x(i,1)-r3x(i,1))*cx-r23c, (r2x(i,1)-r3x(i,1))*cy,       (r2x(i,1)-r3x(i,1))*cz;
                 (r2y(i,1)-r3y(i,1))*cx,      (r2y(i,1)-r3y(i,1))*cy-r23c,  (r2y(i,1)-r3y(i,1))*cz; 
                 (r2z(i,1)-r3z(i,1))*cx,      (r2z(i,1)-r3z(i,1))*cy,       (r2z(i,1)-r3z(i,1))*cz-r23c]')*[gmx(i,1);gmy(i,1);gmz(i,1)] ;
   
  T2 = (r123c - [(r3x(i,1)-r1x(i,1))*cx-r31c   (r3x(i,1)-r1x(i,1))*cy      (r3x(i,1)-r1x(i,1))*cz        ;
                 (r3y(i,1)-r1y(i,1))*cx        (r3y(i,1)-r1y(i,1))*cy-r31c (r3y(i,1)-r1y(i,1))*cz        ; 
                 (r3z(i,1)-r1z(i,1))*cx        (r3z(i,1)-r1z(i,1))*cy      (r3z(i,1)-r1z(i,1))*cz-r31c]')*[gmx(i,1);gmy(i,1);gmz(i,1)] ;   
            
  T3 = (r123c - [(r1x(i,1)-r2x(i,1))*cx-r12c   (r1x(i,1)-r2x(i,1))*cy       (r1x(i,1)-r2x(i,1))*cz       ;
                 (r1y(i,1)-r2y(i,1))*cx        (r1y(i,1)-r2y(i,1))*cy-r12c  (r1y(i,1)-r2y(i,1))*cz        ; 
                 (r1z(i,1)-r2z(i,1))*cx        (r1z(i,1)-r2z(i,1))*cy       (r1z(i,1)-r2z(i,1))*cz-r12c]')*[gmx(i,1);gmy(i,1);gmz(i,1)] ;
  matc = [0  -cz  cy;
          cz  0  -cx;
         -cy  cx  0 ];
  mat1 = matc*[ 0                   -(r2z(i,1)-r3z(i,1))  (r2y(i,1)-r3y(i,1)) ;
                (r2z(i,1)-r3z(i,1))           0          -(r2x(i,1)-r3x(i,1)) ;
               -(r2y(i,1)-r3y(i,1))  (r2x(i,1)-r3x(i,1))            0         ];  
           
  mat2 = matc*[ 0                   -(r3z(i,1)-r1z(i,1))  (r3y(i,1)-r1y(i,1)) ;
                (r3z(i,1)-r1z(i,1))           0          -(r3x(i,1)-r1x(i,1)) ;
               -(r3y(i,1)-r1y(i,1))  (r3x(i,1)-r1x(i,1))            0         ];  
  
  mat3 = matc*[ 0                   -(r1z(i,1)-r2z(i,1))  (r1y(i,1)-r2y(i,1)) ;
                (r1z(i,1)-r2z(i,1))           0          -(r1x(i,1)-r2x(i,1)) ;
               -(r1y(i,1)-r2y(i,1))  (r1x(i,1)-r2x(i,1))            0         ];     
           
           
  %=== with r1 ===          
  F(ind(1),1) =  -sum((r1x(i,1)-r1x(j,1))./d11 + (r1x(i,1)-r2x(j,1))./d12 + (r1x(i,1)-r3x(j,1))./d13)...
                  -rho(i,1)*(r1x(i,1) + r2x(i,1) + r3x(i,1)) - (lm12(i,1)*(r1x(i,1)-r2x(i,1))+ lm13(i,1)*(r1x(i,1)-r3x(i,1))) + T1(1);
                    
  F(ind(2),1) =  -sum((r1y(i,1)-r1y(j,1))./d11 + (r1y(i,1)-r2y(j,1))./d12 + (r1y(i,1)-r3y(j,1))./d13)...
                  -rho(i,1)*(r1y(i,1) + r2y(i,1) + r3y(i,1)) - (lm12(i,1)*(r1y(i,1)-r2y(i,1))+ lm13(i,1)*(r1y(i,1)-r3y(i,1))) + T1(2);
              
  F(ind(3),1) =  -sum((r1z(i,1)-r1z(j,1))./d11 + (r1z(i,1)-r2z(j,1))./d12 + (r1z(i,1)-r3z(j,1))./d13)...
                  -rho(i,1)*(r1z(i,1) + r2z(i,1) + r3z(i,1)) - (lm12(i,1)*(r1z(i,1)-r2z(i,1))+ lm13(i,1)*(r1z(i,1)-r3z(i,1))) + T1(3);   
  %== with r2 ===
  F(ind(4),1) =  -sum((r2x(i,1)-r1x(j,1))./d21 + (r2x(i,1)-r2x(j,1))./d22 + (r2x(i,1)-r3x(j,1))./d23)...
                  -rho(i,1)*(r1x(i,1) + r2x(i,1) + r3x(i,1)) - (lm12(i,1)*(r2x(i,1)-r1x(i,1))+ lm23(i,1)*(r2x(i,1)-r3x(i,1))) + T2(1);  
              
  F(ind(5),1) =  -sum((r2y(i,1)-r1y(j,1))./d21 + (r2y(i,1)-r2y(j,1))./d22 + (r2y(i,1)-r3y(j,1))./d23)...
                  -rho(i,1)*(r1y(i,1) + r2y(i,1) + r3y(i,1)) - (lm12(i,1)*(r2y(i,1)-r1y(i,1))+ lm23(i,1)*(r2y(i,1)-r3y(i,1))) + T2(2); 
              
  F(ind(6),1) =  -sum((r2z(i,1)-r1z(j,1))./d21 + (r2z(i,1)-r2z(j,1))./d22 + (r2z(i,1)-r3z(j,1))./d23)...
                  -rho(i,1)*(r1z(i,1) + r2z(i,1) + r3z(i,1)) - (lm12(i,1)*(r2z(i,1)-r1z(i,1))+ lm23(i,1)*(r2z(i,1)-r3z(i,1))) + T2(3);  
              
  %=== with r3 ==
  F(ind(7),1) =  -sum((r3x(i,1)-r1x(j,1))./d31 + (r3x(i,1)-r2x(j,1))./d32 + (r3x(i,1)-r3x(j,1))./d33)...
                  -rho(i,1)*(r1x(i,1) + r2x(i,1) + r3x(i,1)) - (lm13(i,1)*(r3x(i,1)-r1x(i,1))+ lm23(i,1)*(r3x(i,1)-r2x(i,1))) + T3(1);  
  
  F(ind(8),1) =  -sum((r3y(i,1)-r1y(j,1))./d31 + (r3y(i,1)-r2y(j,1))./d32 + (r3y(i,1)-r3y(j,1))./d33)...
                  -rho(i,1)*(r1y(i,1) + r2y(i,1) + r3y(i,1)) - (lm13(i,1)*(r3y(i,1)-r1y(i,1))+ lm23(i,1)*(r3y(i,1)-r2y(i,1))) + T3(2);
              
  F(ind(9),1) =  -sum((r3z(i,1)-r1z(j,1))./d31 + (r3z(i,1)-r2z(j,1))./d32 + (r3z(i,1)-r3z(j,1))./d33)...
                  -rho(i,1)*(r1z(i,1) + r2z(i,1) + r3z(i,1)) - (lm13(i,1)*(r3z(i,1)-r1z(i,1))+ lm23(i,1)*(r3z(i,1)-r2z(i,1))) + T3(3);            
              
              
  
  
  %==== jacobian evaluation === 
  %== diagonal term ===
  %-- f_r1_r1
  
  
  %==============================================================
  %=== Differentiation of r_i1 terms 
  
  fri1xrj1x = (r1x(i,1)-r1x(j,1)).*(r1x(i,1)-r1x(j,1))./d11_5;
  fri1xrj2x = (r1x(i,1)-r2x(j,1)).*(r1x(i,1)-r2x(j,1))./d12_5;
  fri1xrj3x = (r1x(i,1)-r3x(j,1)).*(r1x(i,1)-r3x(j,1))./d13_5;
  
  fri1xrj1y = (r1x(i,1)-r1x(j,1)).*(r1y(i,1)-r1y(j,1))./d11_5;
  fri1xrj2y = (r1x(i,1)-r2x(j,1)).*(r1y(i,1)-r2y(j,1))./d12_5;
  fri1xrj3y = (r1x(i,1)-r3x(j,1)).*(r1y(i,1)-r3y(j,1))./d13_5;
  
  fri1xrj1z = (r1x(i,1)-r1x(j,1)).*(r1z(i,1)-r1z(j,1))./d11_5;
  fri1xrj2z = (r1x(i,1)-r2x(j,1)).*(r1z(i,1)-r2z(j,1))./d12_5;
  fri1xrj3z = (r1x(i,1)-r3x(j,1)).*(r1z(i,1)-r3z(j,1))./d13_5;
  %==============================================================
  
  fri1yrj1x = (r1y(i,1)-r1y(j,1)).*(r1x(i,1)-r1x(j,1))./d11_5;
  fri1yrj2x = (r1y(i,1)-r2y(j,1)).*(r1x(i,1)-r2x(j,1))./d12_5;
  fri1yrj3x = (r1y(i,1)-r3y(j,1)).*(r1x(i,1)-r3x(j,1))./d13_5;
  
  fri1yrj1y = (r1y(i,1)-r1y(j,1)).*(r1y(i,1)-r1y(j,1))./d11_5;
  fri1yrj2y = (r1y(i,1)-r2y(j,1)).*(r1y(i,1)-r2y(j,1))./d12_5;
  fri1yrj3y = (r1y(i,1)-r3y(j,1)).*(r1y(i,1)-r3y(j,1))./d13_5;
  
  fri1yrj1z = (r1y(i,1)-r1y(j,1)).*(r1z(i,1)-r1z(j,1))./d11_5;
  fri1yrj2z = (r1y(i,1)-r2y(j,1)).*(r1z(i,1)-r2z(j,1))./d12_5;
  fri1yrj3z = (r1y(i,1)-r3y(j,1)).*(r1z(i,1)-r3z(j,1))./d13_5;
  
 %==============================================================
 
  fri1zrj1x = (r1z(i,1)-r1z(j,1)).*(r1x(i,1)-r1x(j,1))./d11_5;
  fri1zrj2x = (r1z(i,1)-r2z(j,1)).*(r1x(i,1)-r2x(j,1))./d12_5;
  fri1zrj3x = (r1z(i,1)-r3z(j,1)).*(r1x(i,1)-r3x(j,1))./d13_5;
  
  fri1zrj1y = (r1z(i,1)-r1z(j,1)).*(r1y(i,1)-r1y(j,1))./d11_5;
  fri1zrj2y = (r1z(i,1)-r2z(j,1)).*(r1y(i,1)-r2y(j,1))./d12_5;
  fri1zrj3y = (r1z(i,1)-r3z(j,1)).*(r1y(i,1)-r3y(j,1))./d13_5;
  
  fri1zrj1z = (r1z(i,1)-r1z(j,1)).*(r1z(i,1)-r1z(j,1))./d11_5;
  fri1zrj2z = (r1z(i,1)-r2z(j,1)).*(r1z(i,1)-r2z(j,1))./d12_5;
  fri1zrj3z = (r1z(i,1)-r3z(j,1)).*(r1z(i,1)-r3z(j,1))./d13_5;
  
  %
  J(ind(1),ind(1)) = -sum(1./d11  + 1./d12  + 1./d13) +3*sum(fri1xrj1x + fri1xrj2x +fri1xrj3x) -rho(i,1) -(lm12(i,1)+lm13(i,1));
  J(ind(1),ind(2)) =                                  +3*sum(fri1xrj1y + fri1xrj2y +fri1xrj3y);
  J(ind(1),ind(3)) =                                  +3*sum(fri1xrj1z + fri1xrj2z +fri1xrj3z); 
  
  J(ind(2),ind(1)) =                                  +3*sum(fri1yrj1x + fri1yrj2x +fri1yrj3x);                                                     
  J(ind(2),ind(2)) = -sum(1./d11  + 1./d12  + 1./d13) +3*sum(fri1yrj1y + fri1yrj2y +fri1yrj3y) -rho(i,1) -(lm12(i,1)+lm13(i,1)); 
  J(ind(2),ind(3)) =                                  +3*sum(fri1yrj1z + fri1yrj2z +fri1yrj3z);  
 
   %
  J(ind(3),ind(1)) =                                  +3*sum(fri1zrj1x + fri1zrj2x +fri1zrj3x);
  J(ind(3),ind(2)) =                                  +3*sum(fri1zrj1y + fri1zrj2y +fri1zrj3y); 
  J(ind(3),ind(3)) = -sum(1./d11  + 1./d12  + 1./d13) +3*sum(fri1zrj1z + fri1zrj2z +fri1zrj3z) -rho(i,1) -(lm12(i,1)+lm13(i,1)); 
  
  %===========  Terms due to edge length constraint $(ri-rj)^2-a^2 = 0 $
  %=========
  J(ind(1:3),ind(4:6)) = J(ind(1:3),ind(4:6))  + [lm12(i,1) 0 0;0 lm12(i,1) 0;0 0 lm12(i,1)]; % (r1-r2)^2-a^2 
  J(ind(1:3),ind(7:9)) = J(ind(1:3),ind(7:9))  + [lm13(i,1) 0 0;0 lm13(i,1) 0;0 0 lm13(i,1)]; % (r1-r3)^2-a^2 

 %=================================================================================================================================================
 %=================================================================================================================================================
 % .  f_r2_r2
 %==============================================================
  %=== Differentiation of r_i1 terms 
  
  fri2xrj1x = (r2x(i,1)-r1x(j,1)).*(r2x(i,1)-r1x(j,1))./d21_5;
  fri2xrj2x = (r2x(i,1)-r2x(j,1)).*(r2x(i,1)-r2x(j,1))./d22_5;
  fri2xrj3x = (r2x(i,1)-r3x(j,1)).*(r2x(i,1)-r3x(j,1))./d23_5;
  
  fri2xrj1y = (r2x(i,1)-r1x(j,1)).*(r2y(i,1)-r1y(j,1))./d21_5;
  fri2xrj2y = (r2x(i,1)-r2x(j,1)).*(r2y(i,1)-r2y(j,1))./d22_5;
  fri2xrj3y = (r2x(i,1)-r3x(j,1)).*(r2y(i,1)-r3y(j,1))./d23_5;
  
  fri2xrj1z = (r2x(i,1)-r1x(j,1)).*(r2z(i,1)-r1z(j,1))./d21_5;
  fri2xrj2z = (r2x(i,1)-r2x(j,1)).*(r2z(i,1)-r2z(j,1))./d22_5;
  fri2xrj3z = (r2x(i,1)-r3x(j,1)).*(r2z(i,1)-r3z(j,1))./d23_5;
  %==============================================================
  
  fri2yrj1x = (r2y(i,1)-r1y(j,1)).*(r2x(i,1)-r1x(j,1))./d21_5;
  fri2yrj2x = (r2y(i,1)-r2y(j,1)).*(r2x(i,1)-r2x(j,1))./d22_5;
  fri2yrj3x = (r2y(i,1)-r3y(j,1)).*(r2x(i,1)-r3x(j,1))./d23_5;
  
  fri2yrj1y = (r2y(i,1)-r1y(j,1)).*(r2y(i,1)-r1y(j,1))./d21_5;
  fri2yrj2y = (r2y(i,1)-r2y(j,1)).*(r2y(i,1)-r2y(j,1))./d22_5;
  fri2yrj3y = (r2y(i,1)-r3y(j,1)).*(r2y(i,1)-r3y(j,1))./d23_5;
  
  fri2yrj1z = (r2y(i,1)-r1y(j,1)).*(r2z(i,1)-r1z(j,1))./d21_5;
  fri2yrj2z = (r2y(i,1)-r2y(j,1)).*(r2z(i,1)-r2z(j,1))./d22_5;
  fri2yrj3z = (r2y(i,1)-r3y(j,1)).*(r2z(i,1)-r3z(j,1))./d23_5;
  
 %==============================================================
 
  fri2zrj1x = (r2z(i,1)-r1z(j,1)).*(r2x(i,1)-r1x(j,1))./d21_5;
  fri2zrj2x = (r2z(i,1)-r2z(j,1)).*(r2x(i,1)-r2x(j,1))./d22_5;
  fri2zrj3x = (r2z(i,1)-r3z(j,1)).*(r2x(i,1)-r3x(j,1))./d23_5;
  
  fri2zrj1y = (r2z(i,1)-r1z(j,1)).*(r2y(i,1)-r1y(j,1))./d21_5;
  fri2zrj2y = (r2z(i,1)-r2z(j,1)).*(r2y(i,1)-r2y(j,1))./d22_5;
  fri2zrj3y = (r2z(i,1)-r3z(j,1)).*(r2y(i,1)-r3y(j,1))./d23_5;
  
  fri2zrj1z = (r2z(i,1)-r1z(j,1)).*(r2z(i,1)-r1z(j,1))./d21_5;
  fri2zrj2z = (r2z(i,1)-r2z(j,1)).*(r2z(i,1)-r2z(j,1))./d22_5;
  fri2zrj3z = (r2z(i,1)-r3z(j,1)).*(r2z(i,1)-r3z(j,1))./d23_5;
  
    %
  J(ind(4),ind(4)) = -sum(1./d21  + 1./d22  + 1./d23) +3*sum(fri2xrj1x + fri2xrj2x +fri2xrj3x)-rho(i,1)-(lm12(i,1)+lm23(i,1));
  J(ind(4),ind(5)) =                                  +3*sum(fri2xrj1y + fri2xrj2y +fri2xrj3y);
  J(ind(4),ind(6)) =                                  +3*sum(fri2xrj1z + fri2xrj2z +fri2xrj3z); 
  
  J(ind(5),ind(4)) =                                  +3*sum(fri2yrj1x + fri2yrj2x +fri2yrj3x);                                                     
  J(ind(5),ind(5)) = -sum(1./d21  + 1./d22  + 1./d23) +3*sum(fri2yrj1y + fri2yrj2y +fri2yrj3y)-rho(i,1)-(lm12(i,1)+lm23(i,1)); 
  J(ind(5),ind(6)) =                                  +3*sum(fri2yrj1z + fri2yrj2z +fri2yrj3z);  
 
   %
  J(ind(6),ind(4)) =                                  +3*sum(fri2zrj1x + fri2zrj2x +fri2zrj3x);
  J(ind(6),ind(5)) =                                  +3*sum(fri2zrj1y + fri2zrj2y +fri2zrj3y); 
  J(ind(6),ind(6)) = -sum(1./d21  + 1./d22  + 1./d23) +3*sum(fri2zrj1z + fri2zrj2z +fri2zrj3z)-rho(i,1)-(lm12(i,1)+lm23(i,1)); 
                                                         
 
  %===========  Terms due to edge length constraint $(ri-rj)^2-a^2 = 0 $

  J(ind(4:6),ind(1:3)) = J(ind(4:6),ind(1:3))  + [lm12(i,1) 0 0;0 lm12(i,1) 0;0 0 lm12(i,1)];    % (r2-r1)^2-a^2=0
  J(ind(4:6),ind(7:9)) = J(ind(4:6),ind(7:9))  + [lm23(i,1) 0 0;0 lm23(i,1) 0;0 0 lm23(i,1)];    % (r2-r3)^2-a^2=0
 %=================================================================================================================================================
 %=================================================================================================================================================
 % .  f_r3_r3
 %==============================================================
  %=== Differentiation of r_i1 terms 
  
  fri3xrj1x = (r3x(i,1)-r1x(j,1)).*(r3x(i,1)-r1x(j,1))./d31_5;
  fri3xrj2x = (r3x(i,1)-r2x(j,1)).*(r3x(i,1)-r2x(j,1))./d32_5;
  fri3xrj3x = (r3x(i,1)-r3x(j,1)).*(r3x(i,1)-r3x(j,1))./d33_5;
  
  fri3xrj1y = (r3x(i,1)-r1x(j,1)).*(r3y(i,1)-r1y(j,1))./d31_5;
  fri3xrj2y = (r3x(i,1)-r2x(j,1)).*(r3y(i,1)-r2y(j,1))./d32_5;
  fri3xrj3y = (r3x(i,1)-r3x(j,1)).*(r3y(i,1)-r3y(j,1))./d33_5;
  
  fri3xrj1z = (r3x(i,1)-r1x(j,1)).*(r3z(i,1)-r1z(j,1))./d31_5;
  fri3xrj2z = (r3x(i,1)-r2x(j,1)).*(r3z(i,1)-r2z(j,1))./d32_5;
  fri3xrj3z = (r3x(i,1)-r3x(j,1)).*(r3z(i,1)-r3z(j,1))./d33_5;
  %==============================================================
  
  fri3yrj1x = (r3y(i,1)-r1y(j,1)).*(r3x(i,1)-r1x(j,1))./d31_5;
  fri3yrj2x = (r3y(i,1)-r2y(j,1)).*(r3x(i,1)-r2x(j,1))./d32_5;
  fri3yrj3x = (r3y(i,1)-r3y(j,1)).*(r3x(i,1)-r3x(j,1))./d33_5;
  
  fri3yrj1y = (r3y(i,1)-r1y(j,1)).*(r3y(i,1)-r1y(j,1))./d31_5;
  fri3yrj2y = (r3y(i,1)-r2y(j,1)).*(r3y(i,1)-r2y(j,1))./d32_5;
  fri3yrj3y = (r3y(i,1)-r3y(j,1)).*(r3y(i,1)-r3y(j,1))./d33_5;
  
  fri3yrj1z = (r3y(i,1)-r1y(j,1)).*(r3z(i,1)-r1z(j,1))./d31_5;
  fri3yrj2z = (r3y(i,1)-r2y(j,1)).*(r3z(i,1)-r2z(j,1))./d32_5;
  fri3yrj3z = (r3y(i,1)-r3y(j,1)).*(r3z(i,1)-r3z(j,1))./d33_5;
  
 %==============================================================
 
  fri3zrj1x = (r3z(i,1)-r1z(j,1)).*(r3x(i,1)-r1x(j,1))./d31_5;
  fri3zrj2x = (r3z(i,1)-r2z(j,1)).*(r3x(i,1)-r2x(j,1))./d32_5;
  fri3zrj3x = (r3z(i,1)-r3z(j,1)).*(r3x(i,1)-r3x(j,1))./d33_5;
  
  fri3zrj1y = (r3z(i,1)-r1z(j,1)).*(r3y(i,1)-r1y(j,1))./d31_5;
  fri3zrj2y = (r3z(i,1)-r2z(j,1)).*(r3y(i,1)-r2y(j,1))./d32_5;
  fri3zrj3y = (r3z(i,1)-r3z(j,1)).*(r3y(i,1)-r3y(j,1))./d33_5;
  
  fri3zrj1z = (r3z(i,1)-r1z(j,1)).*(r3z(i,1)-r1z(j,1))./d31_5;
  fri3zrj2z = (r3z(i,1)-r2z(j,1)).*(r3z(i,1)-r2z(j,1))./d32_5;
  fri3zrj3z = (r3z(i,1)-r3z(j,1)).*(r3z(i,1)-r3z(j,1))./d33_5;
  
    %
  J(ind(7),ind(7)) = -sum(1./d31  + 1./d32  + 1./d33) +3*sum(fri3xrj1x + fri3xrj2x +fri3xrj3x) -rho(i,1)-(lm13(i,1)+lm23(i,1));
  J(ind(7),ind(8)) =                                  +3*sum(fri3xrj1y + fri3xrj2y +fri3xrj3y);
  J(ind(7),ind(9)) =                                  +3*sum(fri3xrj1z + fri3xrj2z +fri3xrj3z); 
  
  J(ind(8),ind(7)) =                                  +3*sum(fri3yrj1x + fri3yrj2x +fri3yrj3x);                                                     
  J(ind(8),ind(8)) = -sum(1./d31  + 1./d32  + 1./d33) +3*sum(fri3yrj1y + fri3yrj2y +fri3yrj3y) -rho(i,1)-(lm13(i,1)+lm23(i,1));
  J(ind(8),ind(9)) =                                  +3*sum(fri3yrj1z + fri3yrj2z +fri3yrj3z);  
 
   %
  J(ind(9),ind(7)) =                                  +3*sum(fri3zrj1x + fri3zrj2x +fri3zrj3x);
  J(ind(9),ind(8)) =                                  +3*sum(fri3zrj1y + fri3zrj2y +fri3zrj3y); 
  J(ind(9),ind(9)) = -sum(1./d31  + 1./d32  + 1./d33) +3*sum(fri3zrj1z + fri3zrj2z +fri3zrj3z) -rho(i,1)-(lm13(i,1)+lm23(i,1)); 
                                                         
 
  
 
  %===========  Terms due to edge length constraint $(ri-rj)^2-a^2 = 0 $

  J(ind(7:9),ind(1:3)) = J(ind(7:9),ind(1:3))  + [lm13(i,1) 0 0;0 lm13(i,1) 0;0 0 lm13(i,1)];    % (r3-r1)^2-a^2=0
  J(ind(7:9),ind(4:6)) = J(ind(7:9),ind(4:6))  + [lm23(i,1) 0 0;0 lm23(i,1) 0;0 0 lm23(i,1)];    % (r3-r2)^2-a^2=0 
 %=================================================================================================================================================
 %=================================================================================================================================================
                                                        
 %=================================================================================================================================================
 %=================================================================================================================================================
 
 %=============================  Nonlocal terms =============
 
 
 for p1 = 1:N-1
     
         
     indj = 9*j(p1)-8:9*j(p1);   % nonlocal indices 
     
     %j1 = j0(1:3);   % index of j_1
     %j2 = j0(4:6);   % index of j_2
     %j3 = j0(7:9);   % index of j_3
      %=================================================================================================================================================
     %=================================================================================================================================================

     %rri1_rj1
     J(ind(1:3),indj(1:3)) =  [1./d11(p1)-3*fri1xrj1x(p1)        -3*fri1xrj1y(p1)           -3*fri1xrj1z(p1);
                                                                 -3*fri1yrj1x(p1)  1./d11(p1)-3*fri1yrj1y(p1)           -3*fri1yrj1z(p1);
                                                                 -3*fri1zrj1x(p1)            -3*fri1zrj1y(p1) 1./d11(p1)-3*fri1zrj1z(p1)];
     %rri1_rj2                               
     J(ind(1:3),indj(4:6)) =  [1./d12(p1)-3*fri1xrj2x(p1)        -3*fri1xrj2y(p1)           -3*fri1xrj2z(p1);
                                                                 -3*fri1yrj2x(p1)  1./d12(p1)-3*fri1yrj2y(p1)           -3*fri1yrj2z(p1);
                                                                 -3*fri1zrj2x(p1)            -3*fri1zrj2y(p1) 1./d12(p1)-3*fri1zrj2z(p1)];                          
                                  
     %rri1_rj3                               
     J(ind(1:3),indj(7:9)) =  [1./d13(p1)-3*fri1xrj3x(p1)        -3*fri1xrj3y(p1)           -3*fri1xrj3z(p1);
                                                                 -3*fri1yrj3x(p1)  1./d13(p1)-3*fri1yrj3y(p1)           -3*fri1yrj3z(p1);
                                                                 -3*fri1zrj3x(p1)            -3*fri1zrj3y(p1) 1./d13(p1)-3*fri1zrj3z(p1)];   
                                    
    %=================================================================================================================================================
    %=================================================================================================================================================
    %=================================================================================================================================================
    %=================================================================================================================================================

     %rri2_rj1
     J(ind(4:6),indj(1:3)) = [1./d21(p1)-3*fri2xrj1x(p1)         -3*fri2xrj1y(p1)           -3*fri2xrj1z(p1);
                                                                 -3*fri2yrj1x(p1)  1./d21(p1)-3*fri2yrj1y(p1)           -3*fri2yrj1z(p1);
                                                                 -3*fri2zrj1x(p1)            -3*fri2zrj1y(p1) 1./d21(p1)-3*fri2zrj1z(p1)];
     %rri2_rj2                               
     J(ind(4:6),indj(4:6)) = [1./d22(p1)-3*fri2xrj2x(p1)         -3*fri2xrj2y(p1)           -3*fri2xrj2z(p1);
                                                                 -3*fri2yrj2x(p1)  1./d22(p1)-3*fri2yrj2y(p1)           -3*fri2yrj2z(p1);
                                                                 -3*fri2zrj2x(p1)            -3*fri2zrj2y(p1) 1./d22(p1)-3*fri2zrj2z(p1)];                          
                                  
     %rri2_rj3                               
     J(ind(4:6),indj(7:9)) = [1./d23(p1)-3*fri2xrj3x(p1)         -3*fri2xrj3y(p1)           -3*fri2xrj3z(p1);
                                                                 -3*fri2yrj3x(p1)  1./d23(p1)-3*fri2yrj3y(p1)           -3*fri2yrj3z(p1);
                                                                 -3*fri2zrj3x(p1)            -3*fri2zrj3y(p1) 1./d23(p1)-3*fri2zrj3z(p1)];   
                                    
    %=================================================================================================================================================
    %=================================================================================================================================================                                 
           
     %rri3_rj1
     J(ind(7:9),indj(1:3)) = [1./d31(p1)-3*fri3xrj1x(p1)        -3*fri3xrj1y(p1)           -3*fri3xrj1z(p1);
                                                                -3*fri3yrj1x(p1)  1./d31(p1)-3*fri3yrj1y(p1)           -3*fri3yrj1z(p1);
                                                                -3*fri3zrj1x(p1)            -3*fri3zrj1y(p1) 1./d31(p1)-3*fri3zrj1z(p1)];
     %rri3_rj2                               
     J(ind(7:9),indj(4:6)) = [1./d32(p1)-3*fri3xrj2x(p1)        -3*fri3xrj2y(p1)           -3*fri3xrj2z(p1);
                                                                -3*fri3yrj2x(p1)  1./d32(p1)-3*fri3yrj2y(p1)           -3*fri3yrj2z(p1);
                                                                -3*fri3zrj2x(p1)            -3*fri3zrj2y(p1) 1./d32(p1)-3*fri3zrj2z(p1)];                          
                                  
     %rri3_rj3                               
     J(ind(7:9),indj(7:9)) = [1./d33(p1)-3*fri3xrj3x(p1)        -3*fri3xrj3y(p1)           -3*fri3xrj3z(p1);
                                                                -3*fri3yrj3x(p1)  1./d33(p1)-3*fri3yrj3y(p1)           -3*fri3yrj3z(p1);
                                                                -3*fri3zrj3x(p1)            -3*fri3zrj3y(p1) 1./d33(p1)-3*fri3zrj3z(p1)];   
                                    
    %=================================================================================================================================================
    %=================================================================================================================================================                                   
           
                         
                                                        
 end
  
 %== Jacobian calculation completed. Now computing constraints and constraint gradient 
  
 %==== Constraint 
  %==== constraints ===
  %--- centroid of the triangle on the sphere 
  F(9*N+i,1) = 1/2*(cx^2 + cy^2 + cz^2 - 9);
  
  %--- sides of the triangle
  F(10*N+i,1) = 1/2*((r1x(i,1)-r2x(i,1))^2 + (r1y(i,1)-r2y(i,1))^2 +(r1z(i,1)-r2z(i,1))^2   -a^2);
  
  F(11*N+i,1) = 1/2*((r1x(i,1)-r3x(i,1))^2 + (r1y(i,1)-r3y(i,1))^2 +(r1z(i,1)-r3z(i,1))^2   -a^2);
  
  F(12*N+i,1) = 1/2*((r2x(i,1)-r3x(i,1))^2 + (r2y(i,1)-r3y(i,1))^2 +(r2z(i,1)-r3z(i,1))^2   -a^2);

 
  %--- traingle in tangential plane of the sphere
  F(13*N+i,1) = r123y*cz - r123z*cy;  
  F(14*N+i,1) = r123z*cx - r123x*cz; 
  F(15*N+i,1) = r123x*cy - r123y*cx;
  
 %=================================================================================================================================================
 %=================================================================================================================================================                                   
 %--- Centroid on sphere 
 
 J(9*N+i,ind(1:3)) = [cx,cy,cz];
 J(9*N+i,ind(4:6)) = [cx,cy,cz]; 
 J(9*N+i,ind(7:9)) = [cx,cy,cz];
 
 %--- Side of the triangle = a 
 %-- r12 
 J(10*N+i,ind(1:3)) =  [(r1x(i,1)-r2x(i,1)) (r1y(i,1)-r2y(i,1)) (r1z(i,1)-r2z(i,1))];
 J(10*N+i,ind(4:6)) = -[(r1x(i,1)-r2x(i,1)) (r1y(i,1)-r2y(i,1)) (r1z(i,1)-r2z(i,1))]; 
 
 %-- r13 
 J(11*N+i,ind(1:3)) =  [(r1x(i,1)-r3x(i,1)) (r1y(i,1)-r3y(i,1)) (r1z(i,1)-r3z(i,1))];
 J(11*N+i,ind(7:9)) = -[(r1x(i,1)-r3x(i,1)) (r1y(i,1)-r3y(i,1)) (r1z(i,1)-r3z(i,1))];
 
  %-- r23 
 J(12*N+i,ind(4:6)) =  [(r2x(i,1)-r3x(i,1)) (r2y(i,1)-r3y(i,1)) (r2z(i,1)-r3z(i,1))];
 J(12*N+i,ind(7:9)) = -[(r2x(i,1)-r3x(i,1)) (r2y(i,1)-r3y(i,1)) (r2z(i,1)-r3z(i,1))];

       
       
         
       
 J(13*N+i,ind(1:3)) =   [ 0     -r123z  r123y] + mat1(1,1:3);  
 J(14*N+i,ind(1:3)) =   [ r123z  0     -r123x] + mat1(2,1:3);  
 J(15*N+i,ind(1:3)) =   [-r123y  r123x  0    ] + mat1(3,1:3);  
 
 J(13*N+i,ind(4:6)) =   [ 0     -r123z  r123y] + mat2(1,1:3);   
 J(14*N+i,ind(4:6)) =   [ r123z  0     -r123x] + mat2(2,1:3); 
 J(15*N+i,ind(4:6)) =   [-r123y  r123x  0    ] + mat2(3,1:3); 

 J(13*N+i,ind(7:9)) =   [ 0     -r123z  r123y] + mat3(1,1:3);   
 J(14*N+i,ind(7:9)) =   [ r123z  0     -r123x] + mat3(2,1:3); 
 J(15*N+i,ind(7:9)) =   [-r123y  r123x  0    ] + mat3(3,1:3);

 
 
%   J = jac;
%  err = sqrt(sum(F.^2)) ;  
 
  end
  
  J(1:16*N,9*N+1:16*N) = -J(9*N+1:16*N,1:16*N)';
  
  
  end
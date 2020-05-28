 function F = fun_lagrangian6(x)



global  N  N1  u v w    len id uni err qd bx2 E count bxN h tau rho lm bxg byg bzg sig
global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p bext J2 fac


 % x = initial_guess_evolution2();
 
     count = count+1 ;
  %----- reading input and constructing bx by bz ----

  bx(2:N,1)     = x(1:3:3*N1-2,1)    ;
  by(2:N,1)     = x(2:3:3*N1-1,1)    ;
  bz(2:N,1)     = x(3:3:3*N1,1)      ;




  bext(1:3,1) = -[bx(N,1);by(N,1);bz(N,1)];
  bext(4:6,1) = -[bx(2,1);by(2,1);bz(2,1)];

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
  fac = 0.001;

  %== For studying the effect of aspect ratio on the lagrange multipliers
  %===

  %fac = asinh(tau*pi*sig)/(4*pi^3*tau^3);


  p = 2;




        f_b(3*p-5,1) = fac*bx4p(p-1,1) ;
        f_b(3*p-4,1) = fac*by4p(p-1,1) ;
        f_b(3*p-3,1) = fac*bz4p(p-1,1)  ;



%
  %==================== Jacobian entry ================================

  ind = 3*p-5:3*p-3;


  %--- Constraint gradient ---

  %--- unit vector constraint --
  jac(p-1,ind) = [bx(p) by(p) bz(p)];
  %--- unit speed constraint --
  jac(N1+p-1,ind)  =  1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)];

   %---- closure ----
  jac(2*N1+2:2*N1+4,ind) =   [ 0                 -(bz(p-1)-bz(p+1))   (by(p-1)-by(p+1))  ;
                               (bz(p-1)-bz(p+1))       0             -(bx(p-1)-bx(p+1))  ;
                              -(by(p-1)-by(p+1))  (bx(p-1)-bx(p+1))           0        ] ;



  p = 3;



        f_b(3*p-5,1) = fac*bx4p(p-1,1);
        f_b(3*p-4,1) = fac*by4p(p-1,1);
        f_b(3*p-3,1) = fac*bz4p(p-1,1);



  %--- Constraint gradient ---

  ind = 3*p-5:3*p-3;
  %--- unit vector constraint --
  jac(p-1,ind) = [bx(p) by(p) bz(p)];
  %--- unit speed constraint --
  jac(N1+p-1,ind-3) = -1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)];
  jac(N1+p-1,ind)   =  1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)];

   %---- closure ----
  jac(2*N1+2:2*N1+4,ind) =   [ 0                 -(bz(p-1)-bz(p+1))   (by(p-1)-by(p+1))  ;
                               (bz(p-1)-bz(p+1))       0             -(bx(p-1)-bx(p+1))  ;
                              -(by(p-1)-by(p+1))  (bx(p-1)-bx(p+1))           0        ] ;
      %======================================================================



    for p = 4:N-2



        f_b(3*p-5,1) = fac*bx4p(p-1,1)  ;
        f_b(3*p-4,1) = fac*by4p(p-1,1)  ;
        f_b(3*p-3,1) = fac*bz4p(p-1,1)  ;




  %==================== Jacobian entry ================================

  ind = 3*p-5:3*p-3;


  %--- Constraint gradient ---

  %--- unit vector constraint --
  jac(p-1,ind) = [bx(p) by(p) bz(p)];
  %--- unit speed constraint --
  jac(N1+p-1,ind-3) = -1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)];
  jac(N1+p-1,ind)   =  1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)];

   %---- closure ----
  jac(2*N1+2:2*N1+4,ind) =   [ 0                 -(bz(p-1)-bz(p+1))   (by(p-1)-by(p+1))  ;
                               (bz(p-1)-bz(p+1))       0             -(bx(p-1)-bx(p+1))  ;
                              -(by(p-1)-by(p+1))  (bx(p-1)-bx(p+1))           0        ] ;



    %======================================================================
    end

 %----- Euler--Lagrange and Jacobian at the boundary points

  p = N-1;




        f_b(3*p-5,1) = fac*bx4p(p-1,1)  ;
        f_b(3*p-4,1) = fac*by4p(p-1,1)  ;
        f_b(3*p-3,1) = fac*bz4p(p-1,1)  ;



  %==================== Jacobian entry ================================

  ind = 3*p-5:3*p-3;


  %--- Constraint gradient ---

  %--- unit vector constraint --
  jac(p-1,ind) = [bx(p) by(p) bz(p)];
  %--- unit speed constraint --
  jac(N1+p-1,ind-3) = -1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)];
  jac(N1+p-1,ind)   =  1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)];

   %---- closure ----
  jac(2*N1+2:2*N1+4,ind) =   [ 0                 -(bz(p-1)-bz(p+1))   (by(p-1)-by(p+1))  ;
                               (bz(p-1)-bz(p+1))       0             -(bx(p-1)-bx(p+1))  ;
                              -(by(p-1)-by(p+1))  (bx(p-1)-bx(p+1))           0        ] ;



    %======================================================================

 %----- Euler--Lagrange and Jacobian at the boundary points

 p = N;




 f_b(3*p-5,1) = fac*bx4p(p-1,1)  ;
 f_b(3*p-4,1) = fac*by4p(p-1,1)  ;
 f_b(3*p-3,1) = fac*bz4p(p-1,1)  ;



 %==================== Jacobian entry ================================

 ind = 3*p-5:3*p-3;


  %--- Constraint gradient ---

  %--- unit vector constraint --
  jac(p-1,ind) = [bx(p) by(p) bz(p)];
  %--- unit speed constraint --
  jac(N1+p-1,ind-3) = -1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)];     % b'(N)
  jac(N1+p-1,ind)   =  1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)];

  jac(2*N1+1,ind)   = -1/h*[bxp(N+1,1) byp(N+1,1) bzp(N+1,1)];  %b'(N+1)

   %---- closure ----
  jac(2*N1+2:2*N1+4,ind) =    [ 0                 -(bz(p-1)-bz(p+1))   (by(p-1)-by(p+1))  ;
                               (bz(p-1)-bz(p+1))       0             -(bx(p-1)-bx(p+1))  ;
                              -(by(p-1)-by(p+1))  (bx(p-1)-bx(p+1))           0        ] ;


%== pseudo inverse and lagrange multipliers at equilibrium
 
F = -(pinv(jac'))*f_b;



end

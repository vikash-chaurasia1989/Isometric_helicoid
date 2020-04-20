   function  stbl = fun_stability8(x)

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

      %=== Difference between fun_jacobian6 and this file is that here we are
      %allowing for the end point b(1) to be free and solving for
      %the Euler Lagrange equation at it
      % Also, we use b(N+1) = -b(1)


     global  N  N1  u v w    f_b    len id uni err qd bx2 E count bxN h tau rho lm bxg byg bzg sig
     global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p bext J2 fac


      %  x = initial_guess3();

          count = count+1 ;
       %----- reading input and constructing bx by bz ----

       bx(1:N,1)     = x(1:3:3*N-2,1)   ;
       by(1:N,1)     = x(2:3:3*N-1,1)   ;
       bz(1:N,1)     = x(3:3:3*N,1)     ;


       rho           = x(3*N+1:4*N,1)   ;
       lm            = x(4*N+1:5*N,1)   ;

       u             = x(5*N+1,1)       ;
       v             = x(5*N+2,1)       ;
       w             = x(5*N+3,1)       ;



       bext(1:3,1) = -[bx(N,1);by(N,1);bz(N,1)];
       bext(4:6,1) = -[bx(2,1);by(2,1);bz(2,1)];

       bx(N+1,1) = -bx(1,1);
       by(N+1,1) = -by(1,1);
       bz(N+1,1) = -bz(1,1);


       om =  [ 0  -w   v   ;
               w   0  -u   ;
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
          lmp(1,1)   = (lm(1,1)-lm(N,1))/h;
          lmp(2:N,1) = (lm(ig,1)-lm(ig-1,1))/(h);


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





        bx4p(1,1) = (bx(3,1) - 4*bx(2,1) + 6*bx(1,1)+4*bx(N,1)-bx(N-1,1))/h^4;
        by4p(1,1) = (by(3,1) - 4*by(2,1) + 6*by(1,1)+4*by(N,1)-by(N-1,1))/h^4;
        bz4p(1,1) = (bz(3,1) - 4*bz(2,1) + 6*bz(1,1)+4*bz(N,1)-bz(N-1,1))/h^4;

        bx4p(2,1) = (bx(4,1) - 4*bx(3,1) + 6*bx(2,1)-4*bx(1,1)+bext(1,1))/h^4;
        by4p(2,1) = (by(4,1) - 4*by(3,1) + 6*by(2,1)-4*by(1,1)+bext(2,1))/h^4;
        bz4p(2,1) = (bz(4,1) - 4*bz(3,1) + 6*bz(2,1)-4*bz(1,1)+bext(3,1))/h^4;

        ig = 3:N-1;

        bx4p(ig,1) = (bx(ig+2,1) - 4*bx(ig+1,1) + 6*bx(ig,1) - 4*bx(ig-1,1) + bx(ig-2,1))/(h^4);
        by4p(ig,1) = (by(ig+2,1) - 4*by(ig+1,1) + 6*by(ig,1) - 4*by(ig-1,1) + by(ig-2,1))/(h^4);
        bz4p(ig,1) = (bz(ig+2,1) - 4*bz(ig+1,1) + 6*bz(ig,1) - 4*bz(ig-1,1) + bz(ig-2,1))/(h^4);

        bx4p(N,1)   = (bext(4,1)   -4*bx(N+1,1)  + 6*bx(N,1)  - 4*bx(N-1,1) + bx(N-2,1) )/h^4;
        by4p(N,1)   = (bext(5,1)  - 4*by(N+1,1)  + 6*by(N,1)  - 4*by(N-1,1) + by(N-2,1) )/h^4;
        bz4p(N,1)   = (bext(6,1)  - 4*bz(N+1,1)  + 6*bz(N,1)  - 4*bz(N-1,1) + bz(N-2,1) )/h^4;


        %---------------------------------------

        f_c = [0;0;0];


       %----- Euler--Lagrange and Jacobian at the boundary points
        fac = 0.001;

       %== For studying the effect of aspect ratio on the lagrange multipliers
       %===

       %fac = asinh(tau*pi*sig)/(4*pi^3*tau^3);
       p = 1;

             gm = om*[bx(2,1)+bx(N,1);by(2,1)+by(N,1);bz(2,1)+bz(N,1)];

             f_b(3*p-2,1) = fac*bx4p(p,1) + lmp(p,1)*bxp(p,1) + lm(p)*bx2p(p,1) - rho(p,1)*bx(p,1) + gm(1);
             f_b(3*p-1,1) = fac*by4p(p,1) + lmp(p,1)*byp(p,1) + lm(p)*by2p(p,1) - rho(p,1)*by(p,1) + gm(2);
             f_b(3*p,1)   = fac*bz4p(p,1) + lmp(p,1)*bzp(p,1) + lm(p)*bz2p(p,1) - rho(p,1)*bz(p,1) + gm(3);


             f_b(3*N+p,1) = 1/2*(bx(p,1)^2 + by(p,1)^2 + bz(p,1)^2 -1)        ;  % 3*N+1


     %
       %==================== Jacobian entry ================================

       ind = 3*p-2:3*p;

       jac(ind,ind)         = id*(6*fac/h^4 + lmp(p,1)/h -2*lm(p)/h^2-rho(p,1))  ;
       jac(ind,ind+3)       = id*(-4*fac/h^4   +lm(p)/h^2)  +om                  ;
       jac(ind,ind+6)       = id*1*fac/h^4                                       ;

       jac(ind,3*N-2:3*N)   = id*(4*fac/h^4  -lm(p,1)/h^2 +lmp(p,1)/h) + om      ;
       jac(ind,3*N-5:3*N-3) = -id*1*fac/h^4                                      ;


       %--- Constraint gradient ---

       %--- unit vector constraint --
       jac(3*N+p,ind) = [bx(p) by(p) bz(p)];
       %--- unit speed constraint --

        %---- closure ----
       jac(5*N+1:5*N+3,ind) =    -[ 0               -(bz(2,1)+ bz(N,1))   (by(2,1)+ by(N,1)) ;
                                  (bz(2,1)+bz(N,1))        0             -(bx(2,1)+bx(N,1))  ;
                                 -(by(2,1)+by(N,1))  (bx(2,1)+bx(N,1))            0        ] ;

       %====== Differentiation with lagrange multipliers

        %-- differentiation with rho

        jac(ind,3*N+p) = -[bx(p) ;by(p) ;bz(p)];

        %-- differentiation with lambda

        jac(ind,4*N+p)     =  1/h*[bxp(p) ;byp(p) ;bzp(p)] +  [bx2p(p) ;by2p(p) ;bz2p(p)]; % with lambda_1
        jac(ind,5*N)       = -1/h*[bxp(p) ;byp(p) ;bzp(p)];                                % with lambda_N
       %-- differentiation with u v w
        jac(ind,5*N+1:5*N+3)       =  -[ 0                 -(bz(2,1)+ bz(N,1))   (by(2,1)+ by(N,1)) ;
                                        (bz(2,1)+bz(N,1))       0              -(bx(2,1)+bx(N,1))  ;
                                       -(by(2,1)+by(N,1))  (bx(2,1)+bx(N,1))           0        ]  ;



         %======================================================================


       p = 2;

             gm = om*[bx(p+1,1)-bx(p-1,1);by(p+1,1)-by(p-1,1);bz(p+1,1)-bz(p-1,1)];


             f_b(3*p-2,1) = fac*bx4p(p,1) + lmp(p,1)*bxp(p,1) + lm(p)*bx2p(p,1) - rho(p,1)*bx(p,1) + gm(1);
             f_b(3*p-1,1) = fac*by4p(p,1) + lmp(p,1)*byp(p,1) + lm(p)*by2p(p,1) - rho(p,1)*by(p,1) + gm(2);
             f_b(3*p,1)   = fac*bz4p(p,1) + lmp(p,1)*bzp(p,1) + lm(p)*bz2p(p,1) - rho(p,1)*bz(p,1) + gm(3);


             f_b(3*N+p,1) = 1/2*(bx(p,1)^2 + by(p,1)^2 + bz(p,1)^2 -1);

             f_b(4*N+p-1,1) = 1/2*(bxp(p,1)^2 + byp(p,1)^2 + bzp(p,1)^2 - tau^2);

            %--- closure -----
     %
             f_c  = f_c       + [by(p-1,1)*bz(p,1) - bz(p-1,1)*by(p,1)    ;
                                 bz(p-1,1)*bx(p,1) - bx(p-1,1)*bz(p,1)    ;
                                 bx(p-1,1)*by(p,1) - by(p-1,1)*bx(p,1) ] ;
     %
       %==================== Jacobian entry ================================

       ind = 3*p-2:3*p;

       jac(ind,ind-3)   = id*(-4*fac/h^4 - lmp(p,1)/h + lm(p)/h^2) - om        ;
       jac(ind,ind)     = id*(6*fac/h^4 + lmp(p,1)/h -2*lm(p)/h^2-rho(p,1))    ;
       jac(ind,ind+3)   = id*(-4*fac/h^4   +lm(p)/h^2)  +om                    ;
       jac(ind,ind+6)   = id*1*fac/h^4                                         ;

       jac(ind,3*N-2:3*N) = -id*1*fac/h^4                                      ;
         %--- Constraint gradient ---

       %--- unit vector constraint --
       jac(3*N+p,ind)      = [bx(p) by(p) bz(p)]                  ;
       %--- unit speed constraint --
       jac(4*N+p-1,ind-3)  =  -1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)] ;
       jac(4*N+p-1,ind)    =   1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)] ;

        %---- closure ----
       jac(5*N+1:5*N+3,ind) =   [ 0                     -(bz(p-1,1)-bz(p+1,1))  (by(p-1,1)-by(p+1,1))  ;
                                  (bz(p-1,1)-bz(p+1,1))       0                -(bx(p-1,1)-bx(p+1,1))  ;
                                 -(by(p-1,1)-by(p+1,1))  (bx(p-1,1)-bx(p+1,1))           0           ] ;


       % ======== Differentiation with Lagrange multipliers ===
        %-- differentiation with rho

        jac(ind,3*N+p) = -[bx(p) ;by(p) ;bz(p)];

        %-- differentiation with lambda
        jac(ind,4*N+p-1)   = -1/h*[bxp(p) ;byp(p) ;bzp(p)];
        jac(ind,4*N+p)     =  1/h*[bxp(p) ;byp(p) ;bzp(p)] +  [bx2p(p) ;by2p(p) ;bz2p(p)];

       %-- differentiation with u v w
        jac(ind,5*N+1:5*N+3)       =   [ 0                 -(bz(p-1,1)-bz(p+1,1))   (by(p-1,1)-by(p+1,1))  ;
                                         (bz(p-1,1)-bz(p+1,1))       0             -(bx(p-1,1)-bx(p+1,1))  ;
                                        -(by(p-1,1)-by(p+1,1))  (bx(p-1,1)-bx(p+1,1))           0        ] ;


         %======================================================================

       %----- Euler--Lagrange and Jacobian at the boundary points

       p = 3;

             gm = om*[bx(p+1,1)-bx(p-1,1);by(p+1,1)-by(p-1,1);bz(p+1,1)-bz(p-1,1)];

             f_b(3*p-2,1) = fac*bx4p(p,1) + lmp(p,1)*bxp(p,1) + lm(p)*bx2p(p,1) - rho(p,1)*bx(p,1) + gm(1);
             f_b(3*p-1,1) = fac*by4p(p,1) + lmp(p,1)*byp(p,1) + lm(p)*by2p(p,1) - rho(p,1)*by(p,1) + gm(2);
             f_b(3*p,1) = fac*bz4p(p,1) + lmp(p,1)*bzp(p,1) + lm(p)*bz2p(p,1) - rho(p,1)*bz(p,1) + gm(3);


             f_b(3*N+p,1) = 1/2*(bx(p,1)^2 + by(p,1)^2 + bz(p,1)^2 -1);

             f_b(4*N+p-1,1) = 1/2*(bxp(p,1)^2 + byp(p,1)^2 + bzp(p,1)^2 - tau^2);

            %--- closure -----

             f_c  = f_c       + [by(p-1,1)*bz(p,1) - bz(p-1,1)*by(p,1) ;
                                 bz(p-1,1)*bx(p,1) - bx(p-1,1)*bz(p,1) ;
                                 bx(p-1,1)*by(p,1) - by(p-1,1)*bx(p,1) ] ;

       %==================== Jacobian entry ================================

       ind = 3*p-2:3*p;

       jac(ind,ind-6)   = id*1*fac/h^4                                             ;
       jac(ind,ind-3)   = id*(-4*fac/h^4 - lmp(p,1)/h + lm(p)/h^2) - om        ;
       jac(ind,ind)     = id*(6*fac/h^4 + lmp(p,1)/h -2*lm(p)/h^2-rho(p,1))        ;
       jac(ind,ind+3)   = id*(-4*fac/h^4   +lm(p)/h^2)    +om  ;
       jac(ind,ind+6)   = id*1*fac/h^4                                             ;



       %--- Constraint gradient ---

       %--- unit vector constraint --
       jac(3*N+p,ind) = [bx(p) by(p) bz(p)];
       %--- unit speed constraint --
       jac(4*N+p-1,ind-3) = -1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)];
       jac(4*N+p-1,ind)   =  1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)];

        %---- closure ----
       jac(5*N+1:5*N+3,ind) =   [ 0                    -(bz(p-1,1)-bz(p+1,1))   (by(p-1,1)-by(p+1,1))  ;
                                 (bz(p-1,1)-bz(p+1,1))            0            -(bx(p-1,1)-bx(p+1,1))  ;
                                -(by(p-1,1)-by(p+1,1))  (bx(p-1,1)-bx(p+1,1))              0        ]  ;

        %-- differentiation with rho

        jac(ind,3*N+p) = -[bx(p) ;by(p) ;bz(p)];

        %-- differentiation with lambda
        jac(ind,4*N+p-1)   = -1/h*[bxp(p) ;byp(p) ;bzp(p)];
        jac(ind,4*N+p)     =  1/h*[bxp(p) ;byp(p) ;bzp(p)] +  [bx2p(p) ;by2p(p) ;bz2p(p)];

       %-- differentiation with u v w
        jac(ind,5*N+1:5*N+3)       =  [ 0                 -(bz(p-1,1)-bz(p+1,1))   (by(p-1,1)-by(p+1,1))  ;
                                         (bz(p-1,1)-bz(p+1,1))       0             -(bx(p-1,1)-bx(p+1,1))  ;
                                           -(by(p-1,1)-by(p+1,1))  (bx(p-1,1)-bx(p+1,1))           0        ] ;


         %======================================================================



         for p = 4:N-2

             gm = om*[bx(p+1,1)-bx(p-1,1);by(p+1,1)-by(p-1,1);bz(p+1,1)-bz(p-1,1)];


             f_b(3*p-2,1) = fac*bx4p(p,1) + lmp(p,1)*bxp(p,1) + lm(p)*bx2p(p,1) - rho(p,1)*bx(p,1) + gm(1);
             f_b(3*p-1,1) = fac*by4p(p,1) + lmp(p,1)*byp(p,1) + lm(p)*by2p(p,1) - rho(p,1)*by(p,1) + gm(2);
             f_b(3*p,1) = fac*bz4p(p,1) + lmp(p,1)*bzp(p,1) + lm(p)*bz2p(p,1) - rho(p,1)*bz(p,1) + gm(3);


             f_b(3*N+p,1) = 1/2*(bx(p,1)^2 + by(p,1)^2 + bz(p,1)^2 -1);

             f_b(4*N+p-1,1) = 1/2*(bxp(p,1)^2 + byp(p,1)^2 + bzp(p,1)^2 - tau^2);

            %--- closure -----

               f_c  = f_c       + [by(p-1,1)*bz(p,1) - bz(p-1,1)*by(p,1) ;
                                   bz(p-1,1)*bx(p,1) - bx(p-1,1)*bz(p,1) ;
                                   bx(p-1,1)*by(p,1) - by(p-1,1)*bx(p,1) ] ;


       %==================== Jacobian entry ================================

       ind = 3*p-2:3*p;

       jac(ind,ind-6)   = id*fac/h^4                                         ;
       jac(ind,ind-3)   = id*(-4*fac/h^4 - lmp(p,1)/h + lm(p)/h^2) - om      ;
       jac(ind,ind)     = id*(6*fac/h^4 + lmp(p,1)/h -2*lm(p)/h^2-rho(p,1))  ;
       jac(ind,ind+3)   = id*(-4*fac/h^4   +lm(p)/h^2)    + om               ;
       jac(ind,ind+6)   = id*1*fac/h^4                                       ;



       %--- Constraint gradient ---

       %--- unit vector constraint --
       jac(3*N+p,ind) = [bx(p) by(p) bz(p)];
       %--- unit speed constraint --
       jac(4*N+p-1,ind-3) = -1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)];
       jac(4*N+p-1,ind)   =  1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)];

        %---- closure ----
       jac(5*N+1:5*N+3,ind) =   [ 0                 -(bz(p-1,1)-bz(p+1,1))   (by(p-1,1)-by(p+1,1))  ;
                                    (bz(p-1,1)-bz(p+1,1))       0             -(bx(p-1,1)-bx(p+1,1))  ;
                                   -(by(p-1,1)-by(p+1,1))  (bx(p-1,1)-bx(p+1,1))           0        ] ;

        %-- differentiation with rho

        jac(ind,3*N+p) = -[bx(p) ;by(p) ;bz(p)];

        %-- differentiation with lambda
        jac(ind,4*N+p-1)   = -1/h*[bxp(p) ;byp(p) ;bzp(p)];
        jac(ind,4*N+p)     =  1/h*[bxp(p) ;byp(p) ;bzp(p)] +  [bx2p(p) ;by2p(p) ;bz2p(p)];

       %-- differentiation with u v w
          jac(ind,5*N+1:5*N+3)       =   [ 0                 -(bz(p-1,1)-bz(p+1,1))   (by(p-1,1)-by(p+1,1))  ;
                                         (bz(p-1,1)-bz(p+1,1))       0             -(bx(p-1,1)-bx(p+1,1))  ;
                                           -(by(p-1,1)-by(p+1,1))  (bx(p-1,1)-bx(p+1,1))           0        ] ;


         %======================================================================
         end

      %----- Euler--Lagrange and Jacobian at the boundary points

       p = N-1;

            gm = om*[bx(p+1,1)-bx(p-1,1);by(p+1,1)-by(p-1,1);bz(p+1,1)-bz(p-1,1)];



             f_b(3*p-2,1) = fac*bx4p(p,1) + lmp(p,1)*bxp(p,1) + lm(p)*bx2p(p,1) - rho(p,1)*bx(p,1) + gm(1);
             f_b(3*p-1,1) = fac*by4p(p,1) + lmp(p,1)*byp(p,1) + lm(p)*by2p(p,1) - rho(p,1)*by(p,1) + gm(2);
             f_b(3*p,1)   = fac*bz4p(p,1) + lmp(p,1)*bzp(p,1) + lm(p)*bz2p(p,1) - rho(p,1)*bz(p,1) + gm(3);


             f_b(3*N+p,1) = 1/2*(bx(p,1)^2 + by(p,1)^2 + bz(p,1)^2 -1);

             f_b(4*N+p-1,1) = 1/2*(bxp(p,1)^2 + byp(p,1)^2 + bzp(p,1)^2 - tau^2);

            %--- closure -----

             f_c  = f_c       + [by(p-1,1)*bz(p,1) - bz(p-1,1)*by(p,1) ;
                                 bz(p-1,1)*bx(p,1) - bx(p-1,1)*bz(p,1) ;
                                 bx(p-1,1)*by(p,1) - by(p-1,1)*bx(p,1) ] ;

       %==================== Jacobian entry ================================

       ind = 3*p-2:3*p;

       jac(ind,ind-6)   = id*fac/h^4                                                  ;
       jac(ind,ind-3)   = id*(-4*fac/h^4 - lmp(p,1)/h + lm(p)/h^2) - om            ;
       jac(ind,ind)     = id*(6*fac/h^4 + lmp(p,1)/h -2*lm(p)/h^2-rho(p,1))        ;
       jac(ind,ind+3)   = id*(-4*fac/h^4   +lm(p)/h^2)             +om                   ;

       jac(ind,1:3)     = id*(-1*fac/h^4 )  ;


       %--- Constraint gradient ---

       %--- unit vector constraint --
       jac(3*N+p,ind) = [bx(p) by(p) bz(p)];
       %--- unit speed constraint --
       jac(4*N+p-1,ind-3) = -1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)];
       jac(4*N+p-1,ind)   =  1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)];

        %---- closure ----
       jac(5*N+1:5*N+3,ind) =   [ 0                 -(bz(p-1,1)-bz(p+1,1))   (by(p-1,1)-by(p+1,1))  ;
                                    (bz(p-1,1)-bz(p+1,1))       0             -(bx(p-1,1)-bx(p+1,1))  ;
                                   -(by(p-1,1)-by(p+1,1))  (bx(p-1,1)-bx(p+1,1))           0        ] ;

        %-- differentiation with rho

        jac(ind,3*N+p) = -[bx(p) ;by(p) ;bz(p)];

        %-- differentiation with lambda
        jac(ind,4*N+p-1)   = -1/h*[bxp(p) ;byp(p) ;bzp(p)];
        jac(ind,4*N+p)     =  1/h*[bxp(p) ;byp(p) ;bzp(p)] +  [bx2p(p) ;by2p(p) ;bz2p(p)];

        %-- differentiation with u v w
        jac(ind,5*N+1:5*N+3)       =  [ 0                 -(bz(p-1,1)-bz(p+1,1))   (by(p-1,1)-by(p+1,1))  ;
                                         (bz(p-1,1)-bz(p+1,1))       0             -(bx(p-1,1)-bx(p+1,1))  ;
                                           -(by(p-1,1)-by(p+1,1))  (bx(p-1,1)-bx(p+1,1))           0        ] ;


         %======================================================================

      %----- Euler--Lagrange and Jacobian at the boundary points

      p = N;

      gm = om*[bx(p+1,1)-bx(p-1,1);by(p+1,1)-by(p-1,1);bz(p+1,1)-bz(p-1,1)];

      f_b(3*p-2,1) = fac*bx4p(p,1) + lmp(p,1)*bxp(p,1) + lm(p)*bx2p(p,1) - rho(p,1)*bx(p,1) + gm(1);
      f_b(3*p-1,1) = fac*by4p(p,1) + lmp(p,1)*byp(p,1) + lm(p)*by2p(p,1) - rho(p,1)*by(p,1) + gm(2);
      f_b(3*p,1) = fac*bz4p(p,1) + lmp(p,1)*bzp(p,1) + lm(p)*bz2p(p,1) - rho(p,1)*bz(p,1) + gm(3);


      f_b(3*N+p,1) = 1/2*(bx(p,1)^2 + by(p,1)^2 + bz(p,1)^2 -1)         ;

      f_b(4*N+p-1,1) = 1/2*(bxp(p,1)^2 + byp(p,1)^2 + bzp(p,1)^2 - tau^2) ;

      %--- closure -----

      f_c  = f_c       + [  by(p-1,1)*bz(p,1) - bz(p-1,1)*by(p,1)   ;
                            bz(p-1,1)*bx(p,1) - bx(p-1,1)*bz(p,1)   ;
                            bx(p-1,1)*by(p,1) - by(p-1,1)*bx(p,1) ] ;


      f_b(5*N+1:5*N+3,1) =  f_c +     [by(N,1)*bz(N+1,1) - bz(N,1)*by(N+1,1)   ;
                                       bz(N,1)*bx(N+1,1) - bx(N,1)*bz(N+1,1)   ;
                                       bx(N,1)*by(N+1,1) - by(N,1)*bx(N+1,1) ] ;


      f_b(5*N,1) = 1/2*(bxp(N+1,1)^2 + byp(N+1,1)^2 + bzp(N+1,1)^2 - tau^2)   ;

      %==================== Jacobian entry ================================

      ind = 3*p-2:3*p;

      jac(ind,ind-6)   = id*fac/h^4                                          ;
      jac(ind,ind-3)   = id*(-4*fac/h^4 - lmp(p,1)/h + lm(p)/h^2) - om       ;
      jac(ind,ind)     = id*(6*fac/h^4 + lmp(p,1)/h -2*lm(p)/h^2-rho(p,1))   ;

      jac(ind,1:3)     = id*(4*fac/h^4-lm(p)/h^2) - om                       ;
      jac(ind,4:6)     =-id*fac/h^4                                          ;

       %jac(5*N+1:5*N+3,3*N-2:3*N) = .5*[0 0 0; 0 0 -1;0 1 0];
       %--- Constraint gradient ---

       %--- unit vector constraint --
       jac(3*N+p,ind) = [bx(p) by(p) bz(p)];
       %--- unit speed constraint --
       jac(4*N+p-1,ind-3) = -1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)];     % b'(N)
       jac(4*N+p-1,ind)   =  1/(h)*[bxp(p,1) byp(p,1) bzp(p,1)];


       jac(5*N,ind)   = -1/h*[bxp(N+1,1) byp(N+1,1) bzp(N+1,1)];  %b'(N+1) differentiated with bN
       jac(5*N,1:3)   = -1/h*[bxp(N+1,1) byp(N+1,1) bzp(N+1,1)];  %b'(N+1) differentiated with b1

        %---- closure ----
       jac(5*N+1:5*N+3,ind) =    [ 0                 -(bz(p-1,1)-bz(p+1,1))   (by(p-1,1)-by(p+1,1))  ;
                                    (bz(p-1,1)-bz(p+1,1))       0             -(bx(p-1,1)-bx(p+1,1))  ;
                                   -(by(p-1,1)-by(p+1,1))  (bx(p-1,1)-bx(p+1,1))           0        ] ;

        %-- differentiation with rho

        jac(ind,3*N+p) = -[bx(p) ;by(p) ;bz(p)];

        %-- differentiation with lambda
        jac(ind,4*N+p-1)   = -1/h*[bxp(p) ;byp(p) ;bzp(p)];
        jac(ind,4*N+p)     =  1/h*[bxp(p) ;byp(p) ;bzp(p)] +  [bx2p(p) ;by2p(p) ;bz2p(p)];

        %-- differentiation with u v w
        jac(ind,5*N+1:5*N+3)       =  [ 0                 -(bz(p-1,1)-bz(p+1,1))   (by(p-1,1)-by(p+1,1))  ;
                                         (bz(p-1,1)-bz(p+1,1))       0             -(bx(p-1,1)-bx(p+1,1))  ;
                                           -(by(p-1,1)-by(p+1,1))  (bx(p-1,1)-bx(p+1,1))           0        ] ;

        % jac(:,5*N+1:5*N+3) = -jac(:,5*N+1:5*N+3);



 %============== So far we have calculated jacobian of the function
 %== Now we compute the minimum eigen value and check if its positive or
 %negative ---

 % % %
% % %============================= Stability check ============================
%
num_eq = 3*N;    % number of equilibrium equations
num_C  = 2*N+3;   % number of constraints
%
A = (jac(num_eq+1:end,1:num_eq))';   % num_C x num_eq constraint gradient matrix

%--- Qr decomposition of the constraint matrix to obtain null space

[Q R] = qr(A)                ;
Z     = Q(:,num_C+1:num_eq)  ;


%   %---- Projected hessian ----
%
      Ac = Z'*jac(1:num_eq,1:num_eq)*Z      ;
%
%   %--- eigsen values and eigsen vectors ----
%
     [V,D] =  eigs(Ac,10,'smallestabs') ;
%
     [B1, I ] = sort(diag(D),'ascend');

     B1
     B1 = sort(real(B1))  ;

     B = real(B1(1));

     tolB = 10^-8;

     if(abs(B)<tolB)
         stbl=1;
     else

         if(B>0)
             stbl = 1;
         else
             stbl =-1;

         end

     end



   end

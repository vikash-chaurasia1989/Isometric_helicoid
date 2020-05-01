  function f =   initial_evolution3()


global   N1  N    id   h    lm  tau   p1 fac p q

global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p nfold b0x b0y b0z b2x b2y b2z sig E
format longE



path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_evolution/';


nfold = 3;

%N = 72;
tau = 16;
p1 = 1;
N = 105;
sig = 0.01;

tau2 = tau;

str1 =['5pi_knot_N' num2str(N) '_tau_' num2str(10^10*tau2) '.txt'];
str2 = ['3fold_N' num2str(N) '_tau_' num2str(10^10*tau2) '.txt'];

strb1 =  [path 'b_' str1] ;
strb2 =  [path 'b_' str2] ;
%==== load b ===========
temp = load(strb1);

bx1 = temp(:,1);
by1 = temp(:,2);
bz1 = temp(:,3);
    
%==== load b ===========
temp = load(strb2);

bx2 = temp(:,1);
by2 = temp(:,2);
bz2 = temp(:,3);


%=== Isotopy between 1 and 2 ====
t = 0.5;

bx = (1-t)*bx1 + t*bx2;
by = (1-t)*by1 + t*by2;
bz = (1-t)*bz1 + t*bz2;

norm = sqrt(bx.^2+by.^2+bz.^2);

% bx = bx./norm;
% by = by./norm;
% bz = bz./norm;
% 


%== energy of the configuration ===

E = .5*(energy_b(bx1,by1,bz1) + energy_b(bx2,by2,bz2));%1.9*energy_b(bx,by,bz);

n1 = tau/2/pi;
del = 0*0.5*(asinh(n1*pi*sig)/(4*pi^3*n1^3))   ;

 

b(1:3:3*N1-2,1) = bx(2:N,1)  ;
b(2:3:3*N1-1,1) = by(2:N,1)  ;
b(3:3:3*N1,1)   = bz(2:N,1)  ;
 
 f = b(1:3*N1,1);

%--- boundary points --

bx(1,1) = 1;
by(1,1) = 0;
bz(1,1) = 0;

bx(N+1,1) = -1;
by(N+1,1) =  0;
bz(N+1,1) =  0;


%------------ Array initialization ----------------

bxp = zeros(N+1,1);
byp = zeros(N+1,1);
bzp = zeros(N+1,1);

bx2p = zeros(N+1,1);
by2p = zeros(N+1,1);
bz2p = zeros(N+1,1);

bx4p = zeros(N-1,1);
by4p = zeros(N-1,1);
bz4p = zeros(N-1,1);



%-- hessian
id = [1 0 0;0 1 0;0 0 1];


 

  end
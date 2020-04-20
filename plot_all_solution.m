%--- 
clear all
clc

format longE
%-- This is the final file for calculating Linking number
%-- Here, I have used codes developed by Zin Arai and modified it for my
%purpose. The accuracy is impressive. Way better than the double integrals
%I was using to calculate Writhe and then the linking number using the
%Caligerano theorem  Lk = 2*(\tau/(2\pi) + Wr);

%---- read r and b file and create ruling data ----
  col1 = [0,     0.976, 0.968]  ;      % -- color for mode n = 2
  col2 = [0.831, 0.545, 0.247]  ;      % -- color for mode n = 3
  col3 = [0.341, 0.505, 0.819]  ;      % -- color for mode n = 4
  col4 = [0.705, 0.701, 0.070]  ;      % -- color for mode n = 5
  col5 = [0.301, 0.811, 0.498]  ;      % -- color for mode n = 6

  col = {col1,col2,col3,col4,col5} ;

%temp = load('05_fold.txt');


%===== tau for 3 fold ====
tau1 = [8.1:.1:9.5 10:1:19];

tau1 = [8.1:.1:9.5 10:1:11.9 11.97871 11.98071 11.98571 11.99071 11.99571  12:.1:13.9 14:.5:21];

N = 72;
h = 1/N;

N1 = N-1;
%==================

%===== tau for 5 fold ====
% % % 
%     tau1 = [12:.1:13.9 14:.5:25];
% % % 
%     tau1 = [11.97871 11.98071 11.98571 11.99071 11.99571  tau1];
% % % 
%  %   tau1 =  12.1;%11.97871 ;
% % % %  
%       N = 120      ;
% %==================

h = 1/N      ;

% N = length(temp)-1;
% 
%  N = N/2;
% h = 1/N;
 

for p1 = 1:2:length(tau1)-1
    
    tau = tau1(p1);
    
    
    
    %-- 3 fold --
    
    %-- 3 fold --
    
   str0 = ['3fold_N' num2str(N) '_tau_' num2str(10*tau) '.txt'];

     if(tau>11.9)
          str0 = ['3fold_N' num2str(N) '_tau_' num2str(10*tau) '_1.txt'];
     end
    
 

%     %-- 5 fold ----    
%       str0 = ['5pi_knot_N120_tau_' num2str(10*tau) '_2.txt'];
% %  
%       if(tau>13.9)
%               str0 = ['5pi_knot_N120_tau_' num2str(10*tau) '.txt'];
%      end
% %     
   % str0  = '3fold_N72_tau_170.txt';
    
    %-- 5pi --
    %  str0 = ['5pi_knot_N120_tau_' num2str(10*tau1(p1)) '.txt'];
    
    
    strb = ['b_' str0];
    strr = ['r_' str0];
    
    strlm = ['lm_' str0];
    strrho = ['rho_' str0];
    struvw = ['uvw_' str0];
    strbext = ['bext_' str0];
    
    
    %--- load b
    temp = load(strb);

 bx  = temp(1:N+1,1);
 by = temp(1:N+1,2);
 bz = temp(1:N+1,3);
 
  
 %--- rotating about x axis by pi/2 --
 R = [ cos(pi/2) -sin(pi/2);
         sin(pi/2)   cos(pi/2)];
    
    
 for i = 1:N+1
   
     %-- Rotation about x axis by pi/2
     temp = R*[by(i,1);
                     bz(i,1)];
                 
      by(i,1) = temp(1);
      bz(i,1) = temp(2);
      
      %-- Rotation about z axis by pi/2 
      temp = R'*[bx(i,1);
                     by(i,1)];
                 
      bx(i,1) = temp(1);
      by(i,1) = temp(2);
 end
 
 
 dl = sqrt((bx(1:end-1,1)-bx(2:end,1)).^2 + (by(1:end-1,1)-by(2:end,1)).^2 +(bz(1:end-1,1)-bz(2:end,1)).^2);

tau = sum(dl);
 
 
i = 1:N;

tx = by(i,1).*bz(i+1,1) - bz(i,1).*by(i+1,1);
ty = bz(i,1).*bx(i+1,1) - bx(i,1).*bz(i+1,1);
tz = bx(i,1).*by(i+1,1) - by(i,1).*bx(i+1,1);



tx(N+1,1) = tx(1,1);
ty(N+1,1) = ty(1,1);
tz(N+1,1) = tz(1,1);

tx = tx/tau/h;
ty = ty/tau/h;
tz = tz/tau/h;

 
rx(1,1) = 0;
ry(1,1) = 0;
rz(1,1) = 0;

for i=1:N 

    rx(i+1,1) = rx(i,1) + h*tx(i,1);
    ry(i+1,1) = ry(i,1) + h*ty(i,1);
    rz(i+1,1) = rz(i,1) + h*tz(i,1);
end

 nfold = 3;
 ind = [1 1 + (N/nfold:N/nfold:((nfold-1)*N/nfold))];

 
 
  rx = rx -  mean(rx(ind,1)); 
  ry = ry  -  mean(ry(ind,1));  
  rz = rz -   mean(rz(ind,1));
  
 
 w = 0.02;


nx = by.*tz - bz.*ty;
ny = bx.*tz - bz.*tx;
nz = bx.*ty - by.*tx;


% 
 

%-- edge curve and its tangent 

rx0 = rx + w*bx;
ry0 = ry + w*by;
rz0 = rz + w*bz;



tx0 = tx -w*tau*nx;
ty0 = ty -w*tau*ny;
tz0 = tz -w*tau*nz;
 
gamma0(1,:) = rx0'; 
gamma0(2,:) = ry0'; 
gamma0(3,:) = rz0'; 

gamma1(1,:) = rx'; 
gamma1(2,:) = ry'; 
gamma1(3,:) = rz'; 

%-- creating array for linking number subroutine --
 
n= 0;
 

length_0 = N+1;
length_1 = N+1;

 

for i = 1:length_1
   
        for j = 1:length_1
        a = gamma1(:, j)                    + 0     - gamma0(:, i);
        b = gamma1(:, j)                    - gamma0(:, mod(i, length_0) + 1);
        c = gamma1(:, mod(j, length_1) + 1) - gamma0(:, mod(i, length_0) + 1);
        d = gamma1(:, mod(j, length_1) + 1) + 0     - gamma0(:, i);
        n = n + solid_angle_quadrilateral(a, b, c, d);
        
       
        end
               
        
end

 

n = n/(4*pi)*2;

Lk(p1) = n ;

 
%===== plotting the solutions 

figure(1)
plotbrowser on 
plot(rx,ry)
set(gca,'FontSize',25,'LineWidth',1)
legend(num2str(n))
box on
grid on
hold on
  
end
 
 
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

tau1 = 8.1;
N = 72;
%==================

%===== tau for 5 fold ====
% % 
    tau1 = [12:.1:13.9 14:.5:25];
% % 
    tau1 = [11.97871 11.98071 11.98571 11.99071 11.99571  tau1];
% % 
 %   tau1 =  12.1;%11.97871 ;
% % %  
      N = 120      ;
%==================

h = 1/N      ;

% N = length(temp)-1;
% 
%  N = N/2;
% h = 1/N;
 

for p1 = 17:17%length(tau1)
    
    tau = tau1(p1);
    
    
    
    %-- 3 fold --
    
%      str0 = ['3fold_N' num2str(N) '_tau_' num2str(10*tau) '.txt'];

       
    %-- 5 fold ----    
      str0 = ['5pi_knot_N120_tau_' num2str(10*tau) '_1.txt'];
%  
      if(tau>13.9)
              str0 = ['5pi_knot_N120_tau_' num2str(10*tau) '.txt'];
     end
%     
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

 nfold = 5;
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

% 
%   plot3(rx,ry,rz)
%   hold on 
% %   plot3(rx(1,1),ry(1,1),rz(1,1),'ok')
% %   hold on 
% %   plot3(rx(3,1),ry(3,1),rz(3,1),'or')
% %   hold on 
%   plot3(rx0,ry0,rz0,'r')
%   hold on 
%   
%   
%   set(gca,'FontSize',25)
%   set(gcf,'color','w');
%   box on
%   grid on
%   

%-- creating array for linking number subroutine --

for i = 1:N+1
    
    gamma0(1,i) = rx0(i,1);
    gamma0(2,i) = ry0(i,1);
    gamma0(3,i) = rz0(i,1);
    
    gamma1(1,i) = rx(i,1);
    gamma1(2,i) = ry(i,1);
    gamma1(3,i) = rz(i,1);
end
n= 0;

 gamma1 =  gamma1;
 %gamma1 = gamma0;

length_0 = size(gamma0, 2);
length_1 = size(gamma1, 2);

% gamma0 = intval(gamma0);
% gamma1 = intval(gamma1);

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

 %==================================================
 %-------------- Intersection point using InterX subroutine
 %-----------------
  
 
 P = InterX( [rx';ry'],[rx0';ry0']);
 
 xin = P(1,:);   % these are interpolated points. May not exist in our data 
 yin = P(2,:);
 
 %-- extracting points in our data nearest to interpolated data xin and yin

 
 lenP = length(P(2,:));
 indi = zeros(1,lenP);
 indj = indi;
 
 for m1 = 1:lenP
 
     
 temp = (rx-xin(m1)).^2 + (ry-yin(m1)).^2;
 
 [temp2,ind] =min(temp);
 
 indi(m1) = ind ;
 
  temp = (rx0-xin(m1)).^2 + (ry0-yin(m1)).^2;
 
 [temp2,ind] =min(temp);
 
 indj(m1) = ind ;
     
 end
 
 %===== plotting the sign of difference of the vertical height between z
 %coordinates of midline and edge 

 temp = zeros(lenP,1);
 
 temp = (rz(indi,1)-rz0(indj,1))./abs(rz(indi,1)-rz0(indj,1));
% 
%temp = [temp' temp(1,1)]';
figure(p1)
plotbrowser on

plot(temp,'o')
legend(num2str(round(n)))
%   
  plot(temp,'o')
  set(gca,'FontSize',25)
  set(gcf,'color','w');
  box on
  grid on
%   
count = 0;
 
for i = 1:length(temp)-1
    if(temp(i)*temp(i+1)<0)
        count = count+.5;
    end
end

num(p1) = count;
  
end
 
% 
% hold on 
% plot(xin,yin,'o')
% hold on 
% %-- plotting the nearest points in our data to the interpolated points xin
% %and yin 
% 
 
figure(100)
plot(num,'xk')
hold on 
plot(Lk,'or')
hold on 
  set(gca,'FontSize',25)
  set(gcf,'color','w');
  box on
  grid on

 %--- 
clear all
clc
  addpath('/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_3pi')
  addpath('/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_5pi')
  addpath('/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_7pi')


 global    lm  N1 p1 count tau  rho f_b  h tx ty tz  u v w tau1 p1

  global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p  N   len rx ry rz Lk err qd   E bext

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
  col6 = [1, 0.929, 0.278];
  col = {col1,col2,col3,col4,col5,col6} ;
  colw = [41,44,52]/255;

%temp = load('05_fold.txt');
%N = 72;
tau = 16;
 
N = 105;
h = 1/N
N1 = N-1;
sig = 0.01;

%==== Isotopy =====
t1 = linspace(-0.2,0,50);
t2 = linspace(0,.5,50);
t3 = linspace(.5,1,50);
t4 = linspace(1,1.2,50);

t = [t1 t2 t3 t4];


tau2 = tau;

path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_evolution/';

%=== loading stability results ===
 fid = figure;
  
 str = ['animation_' num2str(213)];
 
  % Set up the movie.
writerObj = VideoWriter([path str]); % Name it.
writerObj.FrameRate = 5; % How many frames per second.

open(writerObj); 
 
 
 
 
 rx = zeros(N+1,1);
 ry = zeros(N+1,1);
 rz = zeros(N+1,1);
 
 
 bx = zeros(N+1,1);
 by = zeros(N+1,1);
 bz = zeros(N+1,1);

for p1 = 1:length(t)
   %==== saving data =====
    
   str0 = ['branch_213_N_' num2str(N) '_step_' num2str(p1) '_tau_' num2str(10^10*tau) '.txt'];
    
   %== for evolution ==
    strb = [path 'b_' str0];
    
    
    
    %--- load b
    temp = load(strb);

 bx  = temp(1:N+1,1);
 by = temp(1:N+1,2);
 bz = temp(1:N+1,3);
 
   %--- boundary points --
  
  bx(1,1) = 1;
  by(1,1) = 0;
  bz(1,1) = 0;
  
  bx(N+1,1) = -1;
  by(N+1,1) =  0;
  bz(N+1,1) =  0;
  
  
 
 %dl = sqrt((bx(1:end-1,1)-bx(2:end,1)).^2 + (by(1:end-1,1)-by(2:end,1)).^2 +(bz(1:end-1,1)-bz(2:end,1)).^2);

%tau = sum(dl);
    ig = 2:N;
   
 
   
    bext(1:3,1) = -[bx(N,1);by(N,1);bz(N,1)];
    bext(4:6,1) = -[bx(2,1);by(2,1);bz(2,1)];
    
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
 
i = 1:N+1   ;

% % 
%    tx(i,1) = by(i,1).*bz(i+1,1) - bz(i,1).*by(i+1,1);
%    ty(i,1) = bz(i,1).*bx(i+1,1) - bx(i,1).*bz(i+1,1);
%    tz(i,1) = bx(i,1).*by(i+1,1) - by(i,1).*bx(i+1,1);
% 
tx(i,1) = by(i,1).*bzp(i,1) - bz(i).*byp(i,1);
ty(i,1) = bz(i,1).*bxp(i,1) - bx(i).*bzp(i,1);
tz(i,1) = bx(i,1).*byp(i,1) - by(i).*bxp(i,1);

%
%    tx(N+1,1) = tx(1,1);
%    ty(N+1,1) = ty(1,1);
%    tz(N+1,1) = tz(1,1);
% % %
% %
 tx = tx/tau;
 ty = ty/tau;
 tz = tz/tau;


 
rx(1,1) = 0;
ry(1,1) = 0;
rz(1,1) = 0;

for i=1:N 

    rx(i+1,1) = rx(i,1) + h*tx(i,1);
    ry(i+1,1) = ry(i,1) + h*ty(i,1);
    rz(i+1,1) = rz(i,1) + h*tz(i,1);
end

rx = rx-mean(rx);
ry = ry-mean(ry);
rz = rz-mean(rz);

 

%== rescaling the midline 

l = sum(sqrt((rx(1:end-1,1)-rx(2:end,1)).^2 + (ry(1:end-1,1)-ry(2:end,1)).^2+(rz(1:end-1,1)-rz(2:end,1)).^2));

rx = rx/l;
ry = ry/l;
rz = rz/l;
% %==== reorienting the midline for animation 
% 
%   
% ind = [1 1 + (N/nfold:N/nfold:((nfold-1)*N/nfold))];
%   
%  x1 = rx(ind );
%  y1 = ry(ind );
%  z1 = rz(ind );
%  
%  xc = mean(x1);
%  yc = mean(y1);
%  zc = mean(z1);
% 
% 
%  rx = rx-xc;
%  ry = ry-yc;
%  rz = rz-zc;
%  
%  x1 = x1-xc;
%  y1 = y1-yc;
%  z1 = z1-zc;
% 
%  
% n1 = [x1(1) y1(1) z1(1)];
% n2 = [x1(2) y1(2) z1(2)];
% n3 = [x1(3) y1(3) z1(3)];
% 
% normal = cross(n1 - n2, n1 - n3);
%  
% norm = sum(sqrt(normal.*normal));
% 
% a = normal(1)/norm;
% b = normal(2)/norm;
% c = normal(3)/norm;
% 
% %-- rotating the mid-plane to align with x-y plane 
% %--- axis of rotation u = (u1, u2, u3)
% %  u = n x k  = (b,-a,0)
% u1 =  b;
% u2 = -a;
% 
% th = acos(c);
% 
% Ru = [(cos(th) + u1^2*(1-cos(th))) u1*u2*(1-cos(th))                u2*sin(th);
%        u1*u2*(1-cos(th))           (cos(th) + u2^2*(1-cos(th)))    -u1*sin(th);
%        -u2*sin(th)                 u1*sin(th)                       cos(th)   ] ;
%    
%   
% %-- rotating the mid-line such that the midplane is x-y plane 
% for i = 1:N+1
%     
% temp =  Ru*[rx(i);ry(i);rz(i)] ;
%  
% rx(i) = temp(1);
% ry(i) = temp(2);
% rz(i) = temp(3);
% 
% end
% 

%================================================


figure(fid)

set(fid, 'Position', [100, 100, 1049, 895]);


%plotbrowser on
title('midline')
%hold on
plot3(rx,ry,rz,'color',col{1},'LineWidth',2)

set(gca, 'DataAspectRatio',[1,1,1]);
daspect([1,1,1]);
axis off
set(gcf,'color',colw);

xlim auto
ylim auto
zlim auto

%axis vis3d
%view(normal)
view([1,1,1])


hold off




pause(.01);

frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
writeVideo(writerObj, frame);

             
           
 
end
close(writerObj); 
 

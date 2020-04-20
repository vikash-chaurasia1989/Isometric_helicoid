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


for branch =3:3
 %===== branch 1   
  
 if(branch==1)
%           tau_3 = [8.09575 8.0958  8.09585 8.096 8.097 8.09975 8.1:.1:9.5 9.6:.1:11.9 12:.1:12.5 12.55 12.552 12.554 12.555 12.556 12.56 12.57 12.6:.1:13.2 13.21 13.215 ...
%               13.2165 13.21675 13.216865 13.21686525 13.2168653 13.2168654 13.21686545 13.21686546 13.21686547 13.2168654725    ] ;
%           tau_5 = [13.21686547250001  13.21686547255 13.2168654726 13.21686547275 13.216865473 13.216865474355 13.216865474 13.2168655 13.21686575 13.216866...
%               13.21687 13.2169 13.217 13.2175 13.22 13.2225 13.2235 13.2236 13.223675 13.22375 13.2238];
%      
%           tau_7 = [13.224  13.225 13.23 13.3:.1:13.9 14.5:.5:21 21.5:.5:45];
%      
%      
%      
%           tau1 =      [tau_3 tau_5 tau_7];
%      
   
%         
%           tau_3 = [ 8.0957344694550208E+00  8.09574 8.095745 8.09575 8.0958  8.09585 8.096 8.097 8.09975 8.1:.05:8.5 8.6:.1:9.5 9.6:.1:11.9 12:.1:12.5 12.55 12.552 12.554 12.555 12.556 12.56 12.57 12.6:.1:13.2 13.21 13.215 ...
%          13.2165 13.21675 13.216865 13.21686525 13.2168653 13.2168654 13.21686545 13.21686546 13.21686547 13.2168654725    ] ;
%        tau_5 = [13.21686547250001  13.21686547255 13.2168654726 13.21686547275 13.216865473 13.216865474355 13.216865474 13.2168655 13.21686575 13.216866...
%          13.21687 13.2169 13.217 13.2175 13.22 13.2225 13.2235 13.2236 13.223675 13.22375 13.2238];
%      
%        tau_7 = [13.224  13.225 13.23 13.3 13.31 13.32  13.33 13.3 13.34 13.35  13.36 13.37 13.38 13.39 13.395  13.4:.1:13.9 14.5:.5:21 21.5:.5:33.5 33.6:.1:33.9 33.95 33.96 33.965 33.966 33.9665 33.96675 33.96676 33.96677 ];
% %      
%     
%      tau1 =      [tau_3 tau_5 tau_7];
%      
      


     N = 72;
     N1 = N-1;
     h = 1/N;
     nfold=3;
     
    str1 =  '3fold_N' ;
    path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_3pi_3/';
    
    strtau = 'tau_branch_1.txt';

    tau1 = load([path strtau]);

 elseif(branch==2)
     %================ Branch 2 =====================
     
     %tau1 = [ 11.984235 11.98425  11.9845 11.985 12:.1:25 25.5:.5:34 34.1 34.15 34.155 34.157  34.158 34.159];
     
     
     N = 120;
     N1 = N-1;
     h = 1/N;
     nfold=5;
     
    str1 = '5pi_knot_N' ;
    path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_5pi_3/';
    
     strtau = 'tau_branch_2.txt';

    tau1 = load([path strtau]);
    
     % =============== Branch 3 =====================
     
 else
     tau1 = [15.45:.1:25 25.1:.1:40];
     
     N = 84;
     N1 = N-1;
     h = 1/N;
     
     nfold=7;
    str1 = '7pi_knot_N' ;
    path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_7pi_3/';
      strtau = 'tau_branch_3.txt';

    tau1 = load([path strtau]);
 
 end
%==========================================================================  
 
  
 
 
%==================
 
%==================

h = 1/N      ;

 
  
 Lk = zeros(1,length(tau1));
 st = zeros(1,length(tau1));
 
 
 
 rx = zeros(N+1,1);
 ry = zeros(N+1,1);
 rz = zeros(N+1,1);
 
 
 bx = zeros(N+1,1);
 by = zeros(N+1,1);
 bz = zeros(N+1,1);
 
 
 %=== loading stability results ===
 fid = figure;
  
 str = ['animation_' num2str(branch)];
 
  % Set up the movie.
writerObj = VideoWriter([path str]); % Name it.
writerObj.FrameRate = 5; % How many frames per second.

open(writerObj); 

 
 
for p1 =  1:length(tau1)
    
    sv = 0;
    tau =  tau1(p1);
    
   
    
    %-- 3 fold --
    
      str0 = [str1  num2str(N) '_tau_' num2str((10^10*tau)) '.txt'];
    
   
 
    
    strb = [path 'b_' str0];
    strr = [path  'r_' str0];
    
    strlm = [path  'lm_' str0];
    strrho = [path 'rho_' str0];
    struvw = [path 'uvw_' str0];
    strbext = [path 'bext_' str0];
    
    
    
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


%==== reorienting the midline for animation 

  
ind = [1 1 + (N/nfold:N/nfold:((nfold-1)*N/nfold))];
  
 x1 = rx(ind );
 y1 = ry(ind );
 z1 = rz(ind );
 
 xc = mean(x1);
 yc = mean(y1);
 zc = mean(z1);


 rx = rx-xc;
 ry = ry-yc;
 rz = rz-zc;
 
 x1 = x1-xc;
 y1 = y1-yc;
 z1 = z1-zc;

 
n1 = [x1(1) y1(1) z1(1)];
n2 = [x1(2) y1(2) z1(2)];
n3 = [x1(3) y1(3) z1(3)];

normal = cross(n1 - n2, n1 - n3);
 
norm = sum(sqrt(normal.*normal));

a = normal(1)/norm;
b = normal(2)/norm;
c = normal(3)/norm;

%-- rotating the mid-plane to align with x-y plane 
%--- axis of rotation u = (u1, u2, u3)
%  u = n x k  = (b,-a,0)
u1 =  b;
u2 = -a;

th = acos(c);

Ru = [(cos(th) + u1^2*(1-cos(th))) u1*u2*(1-cos(th))                u2*sin(th);
       u1*u2*(1-cos(th))           (cos(th) + u2^2*(1-cos(th)))    -u1*sin(th);
       -u2*sin(th)                 u1*sin(th)                       cos(th)   ] ;
   
  
%-- rotating the mid-line such that the midplane is x-y plane 
for i = 1:N+1
    
temp =  Ru*[rx(i);ry(i);rz(i)] ;
 
rx(i) = temp(1);
ry(i) = temp(2);
rz(i) = temp(3);

end


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
end

%=== In this file we compare solutions from 3 different branches for same value of \tau and N
%== data to be used in this file are saved in data_evolution 
clear all
clc



  col1 = [0,     0.976, 0.968]  ;      % -- color for mode n = 2
  col2 = [0.831, 0.545, 0.247]  ;      % -- color for mode n = 3
  col3 = [0.341, 0.505, 0.819]  ;      % -- color for mode n = 4
  col4 = [0.705, 0.701, 0.070]  ;      % -- color for mode n = 5
  col5 = [0.301, 0.811, 0.498]  ;      % -- color for mode n = 6
  col6 = [1, 0.929, 0.278];
  col = {col1,col2,col3,col4,col5,col6} ;
  colw = [41,44,52]/255;
  
  
  
  
  N = 105;
  N1 = N-1;
  h = 1/N;
  
  %------------ Array initialization ----------------
  
  bxp = zeros(N+1,1);
  byp = zeros(N+1,1);
  bzp = zeros(N+1,1);
  
    
  tx = zeros(N+1,1);
  ty = zeros(N+1,1);
  tz = zeros(N+1,1);
  
    
  rx = zeros(N+1,1);
  ry = zeros(N+1,1);
  rz = zeros(N+1,1);
  
  


  m1 = [1 2 3];
  
 

  
  figure(1)
  plot_sphere()
  hold on 
  
  for p1 = 1:length(m1)
      
      tau = 16; 
      
  
  
  if(p1==1)
      str0 = ['3fold_N' num2str(N) '_tau_' num2str(10^10*tau) '.txt'];
  elseif(p1==2)
      str0 =['5pi_knot_N' num2str(N) '_tau_' num2str(10^10*tau) '.txt'];
      %
  else
      str0 = ['7pi_knot_N' num2str(N) '_tau_' num2str(10^10*tau) '.txt'];
  end
  
   path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_evolution/';
  strb =  [path 'b_' str0] ;
   
  temp = load(strb);
  bx  = temp(1:N+1,1);
  by  = temp(1:N+1,2);
  bz  = temp(1:N+1,3);
  
  
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
  
  %-- curvature 
  
  k2(p1,:) = 1/tau^2*(bx2p.^2 + by2p.^2 + bz2p.^2)-tau^2;
  
  
  %============= constructing tangent and midline ===
  i = 1:N+1   ;
  %
  tx(i,1) = by(i,1).*bzp(i,1) - bz(i).*byp(i,1);
  ty(i,1) = bz(i,1).*bxp(i,1) - bx(i).*bzp(i,1);
  tz(i,1) = bx(i,1).*byp(i,1) - by(i).*bxp(i,1);
  
  tx = tx/tau;
  ty = ty/tau;
  tz = tz/tau;

%--------------------

% %--- position vector using integration of tangent
% % initialization
%h=1;
rx(1,1) = 0; ry(1,1) = 0; rz(1,1) = 0;

 

for i = 1:N

     rx(i+1,1) = rx(i,1) + h*tx(i+1,1);
     ry(i+1,1) = ry(i,1) + h*ty(i+1,1);
     rz(i+1,1) = rz(i,1) + h*tz(i+1,1);

end

rx = rx-mean(rx);
ry = ry-mean(ry);
rz = rz-mean(rz);


%== rescaling the midline 

l = sum(sqrt((rx(1:end-1,1)-rx(2:end,1)).^2 + (ry(1:end-1,1)-ry(2:end,1)).^2+(rz(1:end-1,1)-rz(2:end,1)).^2));

rx = rx/l;
ry = ry/l;
rz = rz/l;

l_b(p1) = sum(sqrt((bx(1:end-1,1)-bx(2:end,1)).^2 + (by(1:end-1,1)-by(2:end,1)).^2+(bz(1:end-1,1)-bz(2:end,1)).^2));


figure(1)
plotbrowser on
title('binormal')

hold on
plot3(bx,by,bz,'color',col{p1},'LineWidth',2)
hold on
set(gca, 'DataAspectRatio',[1,1,1]);
daspect([1,1,1]);
axis off
set(gcf,'color',colw);
axis vis3d
hold on 


figure(2)
plotbrowser on
title('midline')
hold on
plot3(rx,ry,rz,'color',col{p1},'LineWidth',2)
hold on
set(gca, 'DataAspectRatio',[1,1,1]);
daspect([1,1,1]);
axis off
set(gcf,'color',colw);
axis vis3d
hold on 


%=== Storing midline in the master array 

Rx(p1,:) = rx;
Ry(p1,:) = ry;
Rz(p1,:) = rz;


%=== saving midline for geodesics =====

path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/Geodesics/data/';


str0 = ['branch' num2str(m1(p1)) '_N' num2str(N) '_tau_' num2str(10^10*tau) '.txt'];

strr = ['r_' str0];

 % 
fileID = fopen([path strr],'w');
fprintf(fileID,'%30.16E %30.16E %30.16E \r\n',[rx'    ;ry'    ; rz'    ] );
fclose(fileID);


  end
  
  
%==== Isotopy =====
t1 = linspace(-0.2,0,20);
t2 = linspace(0,.5,20);
t3 = linspace(.5,1,20);
t4 = linspace(1,1.2,10);

t = [t1 t2 t3 t4];

 t = linspace(-.1,1.1,1000 );
t = 1;
for i = 1:length(t)
    
    %=== for isotopy between 1,2,3 curves 
    Iso_x(i,:) = (1-t(i))*(0.5-t(i))/0.5*Rx(2,:)  + (1-t(i))*t(i)/0.25*Rx(1,:) - t(i)*(0.5-t(i))/0.5*Rx(3,:);
    Iso_y(i,:) = (1-t(i))*(0.5-t(i))/0.5*Ry(2,:)  + (1-t(i))*t(i)/0.25*Ry(1,:) - t(i)*(0.5-t(i))/0.5*Ry(3,:);
    Iso_z(i,:) = (1-t(i))*(0.5-t(i))/0.5*Rz(2,:)  + (1-t(i))*t(i)/0.25*Rz(1,:) - t(i)*(0.5-t(i))/0.5*Rz(3,:);
   
    %=== for isotopy between two curves 
%     
%     Iso_x(i,:) = (1-t(i))*Rx(1,:)  + t(i)*Rx(2,:) ;
%     Iso_y(i,:) = (1-t(i))*Ry(1,:)  + t(i)*Ry(2,:) ;
%     Iso_z(i,:) = (1-t(i))*Rz(1,:)  + t(i)*Rz(2,:) ;
    
    rx = Iso_x(i,:);
    ry = Iso_y(i,:);
    rz = Iso_z(i,:);
    
    %=== Derivatives and curvature ==
    rxg = [rx(N) rx rx(2)];
    ryg = [ry(N) ry ry(2)];
    rzg = [rz(N) rz rz(2)];
    
    ig = 2:N+2;
    
    rxp = (rxg(ig)-rxg(ig-1))/h;
    ryp = (ryg(ig)-ryg(ig-1))/h;
    rzp = (rzg(ig)-rzg(ig-1))/h;
    
    rx2p = (rxg(ig+1) +rxg(ig-1) -2*rxg(ig))/h^2;
    ry2p = (ryg(ig+1) +ryg(ig-1) -2*ryg(ig))/h^2;
    rz2p = (rzg(ig+1) +rzg(ig-1) -2*rzg(ig))/h^2;

    %= curvature square ==
    k2 = ((ryp.*rz2p-rzp.*ry2p).^2 + (rzp.*rx2p-rxp.*rz2p).^2 + (rxp.*ry2p-ryp.*rx2p).^2)./(rxp.^2+ryp.^2 + rzp.^2).^3;
    
    E(i) = sum(k2);
% figure(4)
% plotbrowser on
% title('midline')
% hold on
% plot3(rx,ry,rz,'color',col{p1},'LineWidth',2)
% hold on
% set(gca, 'DataAspectRatio',[1,1,1]);
% daspect([1,1,1]);
% axis off
% set(gcf,'color',colw);
% axis vis3d
% hold on 

end
    
figure(4)
plotbrowser on
title('Energy')
hold on
plot(t,E/max(E),'--o','LineWidth',.5)
set(gca,'LineWidth',.5)
  

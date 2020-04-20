 clear all
 clc


%---- read r and b file and create ruling data ----
  col1 = [0,     0.976, 0.968]  ;      % -- color for mode n = 2
  col2 = [0.831, 0.545, 0.247]  ;      % -- color for mode n = 3
  col3 = [0.341, 0.505, 0.819]  ;      % -- color for mode n = 4
  col4 = [0.705, 0.701, 0.070]  ;      % -- color for mode n = 5
  col5 = [0.301, 0.811, 0.498]  ;      % -- color for mode n = 6
  
     
  % col = {'','','r','b','g','k','m'} ;
  col = {col1,col2,col3,col4,col5} ;
%--- midline ---  

tau = 8.1;

tau = 12.1 ;

  str0 = ['3fold_N72_tau_' num2str(10*tau) '.txt'];
    str0 = ['5pi_knot_N120_tau_' num2str(10*tau) '_1.txt'];

  strr = ['r_' str0];
  strb = ['b_' str0];
  


% 
% temp = load(strr);
% 
%   rx = temp(:,1);
%   ry = temp(:,2);
%   rz = temp(:,3);
  

 
    %--- load b
    temp = load(strb);

    N = length(temp)-1;
    h = 1/N;
 bx  = temp(1:N+1,1);
 by = temp(1:N+1,2);
 bz = temp(1:N+1,3);
 

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



 w = 0.02;


nx = by.*tz - bz.*ty;
ny = bx.*tz - bz.*tx;
nz = bx.*ty - by.*tx;


% 

%====== Transforming the midline such that its midlplane is in x-y plane
%and its geometric center is origin 

%--- 
 nfold = 5;
 ind = [1 1 + (N/nfold:N/nfold:((nfold-1)*N/nfold))];

  rx = rx -  mean(rx(ind,1)); 
  ry = ry  -  mean(ry(ind,1));  
  rz = rz -   mean(rz(ind,1));
  
%p1 = [rx(ind(1),1) ry(ind(1),1) rz(ind(1),1)];
%p2 = [rx(ind(2),1) ry(ind(2),1) rz(ind(2),1)];
%p3 = [rx(ind(3),1) ry(ind(3),1) rz(ind(3),1)];


normal = cross( [rx(ind(1),1) ry(ind(1),1) rz(ind(1),1)] - [rx(ind(2),1) ry(ind(2),1) rz(ind(2),1)], [rx(ind(1),1) ry(ind(1),1) rz(ind(1),1)] - [rx(ind(3),1) ry(ind(3),1) rz(ind(3),1)]);
 
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
for i = 1:length(rx) 
    
temp =  Ru*[rx(i,1);
                   ry(i,1);
                   rz(i,1)] ;
 
rx(i,1) = temp(1);
ry(i,1) = temp(2);
rz(i,1) = temp(3);

temp =  Ru*[bx(i,1);
                   by(i,1);
                   bz(i,1)] ;
 
bx(i,1) = temp(1);
by(i,1) = temp(2);
bz(i,1) = temp(3);

temp =  Ru*[tx(i,1);
                   ty(i,1);
                   tz(i,1)] ;
 
tx(i,1) = temp(1);
ty(i,1) = temp(2);
tz(i,1) = temp(3);


temp =  Ru*[nx(i,1);
                   ny(i,1);
                   nz(i,1)] ;
 
nx(i,1) = temp(1);
ny(i,1) = temp(2);
nz(i,1) = temp(3);
 
end

 

%-- edge curve and its tangent 

rx0 = rx + w*bx;
ry0 = ry + w*by;
rz0 = rz + w*bz;



tx0 = tx -w*tau*nx;
ty0 = ty -w*tau*ny;
tz0 = tz -w*tau*nz;

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%                  CLOSED CONFIGURATION H 
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%--- parameter for ruling --
%--- 5pi twist 
wd = .002;

%--- 3pi twist 
%wd = .1;
 v = linspace(-wd,wd,10);

%--- midline ---


 figure(1); 
 plotbrowser on 
 plot3(rx,ry,rz,'color',col{5},'LineWidth',2);
 hold on 
 %-- edges---
 
 rx0 = rx + v(1)*bx;
 ry0 = ry + v(1)*by;
 rz0 = rz + v(1)*bz;
 
 
  plot3(rx0,ry0,rz0,'color',col{2},'LineWidth',2);
 hold on
  
 %-- edges---
 
 rx1 = rx + v(end)*bx;
 ry1 = ry + v(end)*by;
 rz1 = rz + v(end)*bz;
 
 plot3(rx1,ry1,rz1,'color',col{3},'LineWidth',2);

%--- ruling --- 
 
 for i = 1:N 
     
     x = rx(i) + v*bx(i);
     y = ry(i) + v*by(i);
     z = rz(i) + v*bz(i);
     
     plot3(x,y,z,'color',col{1})
     hold on 
 end
 
 set(gca,'FontSize',25)
 set(gcf,'color','w');
box on
grid on
hold on 
%  %---nonorientable edges --
 bx0 = rx(1) + 1.1*v*bx(1);
 by0 = ry(1) + 1.1*v*by(1);
 bz0 = rz(1) + 1.1*v*bz(1);
 
 plot3(bx0,by0,bz0,'y','LineWidth',3.5);
 
 bxN = rx(N) + 1.1*v*bx(N);
 byN = ry(N) + 1.1*v*by(N);
 bzN = rz(N) + 1.1*v*bz(N);
 
 plot3(bxN,byN,bzN,'k','LineWidth',3.5);
 
 
 
%  
%  %-------- ploting helicoid -----
%  
%  figure(1);
%  plotbrowser on 
%  om = 8.0935;
%  n = 8.0935/(2*pi);
%  
%  s = linspace(0,1,N+1);
%  Nz = ones(length(v),1);
%  
%  
%  %--- plotting the edges --
%  
%  %-- midline --
%  
%  plot3(0*s,0*s,s,'color',col{5},'LineWidth',4);
%  hold on 
%  %--- edges
%  
%  hx0 = v(1)*cos(om*s);
%  hy0 = v(1)*sin(om*s);
%  hz0 = s             ;
%  
%  plot3(hx0,hy0,hz0,'color',col{2},'LineWidth',4);
%  hold on
%  
%   %--- edges
%  
%  hx1 = v(end)*cos(om*s);
%  hy1 = v(end)*sin(om*s);
%  hz1 = s             ;
%  
%  plot3(hx1,hy1,hz1,'color',col{3},'LineWidth',4);
%  hold on
% 
%  
%  
%  for i = 1:2:N+1
%      
%      hx = v*cos(om*s(i));
%      hy = v*sin(om*s(i));
%      hz = Nz*s(i); 
%      
%      plot3(hx,hy,hz,'color',col{1})
%      hold on 
%  end 
%  
 
 figure(2); 
 plotbrowser on 
 plot3(rx,ry,rz,'color',col{5},'LineWidth',2);
   set(gca,'FontSize',25)
 set(gcf,'color','w');
box on
grid on
hold on 
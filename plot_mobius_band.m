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

%tau = 12.1 ;

  str0 = ['3fold_N72_tau_' num2str(10*tau) '.txt'];
 %   str0 = ['5pi_knot_N120_tau_' num2str(10*tau) '_1.txt'];

  strr = ['r_' str0];
  strb = ['b_' str0];
  


% 
% temp = load(strr);
% 
%   rx = temp(:,1);
%   ry = temp(:,2);
%   rz = temp(:,3);
  

%--- binormal ---

temp = load(strb);

  bx = temp(:,1);
  by = temp(:,2);
  bz = temp(:,3);

  N = length(bx)-1;

  dl = sqrt((bx(1:end-1,1)-bx(2:end,1)).^2 + (by(1:end-1,1)-by(2:end,1)).^2 +(bz(1:end-1,1)-bz(2:end,1)).^2);

tau = sum(dl);
 h = 1/N;
 
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
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%                  CLOSED CONFIGURATION H 
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%--- parameter for ruling --
%--- 5pi twist 
wd = .005;

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
 



%=== saving edge data == 
 
  str1 =  ['3fold_N72_tau_' num2str(10*tau) '_edge1.txt'] ;
  str2 =  ['3fold_N72_tau_' num2str(10*tau) '_edge2.txt'] ;

  str12 = ['3fold_N72_tau_' num2str(10*tau) '_edge.txt'] ;
  
  
  for i = 1:N 
      
      p1 = 1+ 4*(i-1):4*i;
      
      %=====================
      data(p1(1),:) = [rx0(i,1) ry0(i,1) rz0(i,1)];
      data(p1(2),:) = [rx1(i,1) ry1(i,1) rz1(i,1)];
      
      data(p1(3),:) = [rx1(i+1,1) ry1(i+1,1) rz1(i+1,1)];
      data(p1(4),:) = [rx0(i+1,1) ry0(i+1,1) rz0(i+1,1)];
    
      
      
  end

  
  
  plot3(data(:,1),data(:,2),data(:,3))
   
  fileID = fopen(str12,'w');
  fprintf(fileID,'%30.16E %30.16E %30.16E \r\n',[data(:,1)'    ;data(:,2)'    ; data(:,3)'    ] );
  fclose(fileID);
%   
%   
%   fileID = fopen(str2,'w');
%   fprintf(fileID,'%30.16E %30.16E %30.16E \r\n',[rx1'    ;ry1'    ; rz1'   ] );
%   fclose(fileID);
%   
%   
%   
%  figure(2); 
%  plotbrowser on 
%  plot3(rx,ry,rz,'color',col{5},'LineWidth',2);
%    set(gca,'FontSize',25)
%  set(gcf,'color','w');
% box on
% grid on
% hold on 
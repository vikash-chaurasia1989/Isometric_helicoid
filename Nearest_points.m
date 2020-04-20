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
% 
  tau1 = [12:.1:13.9 14:.5:25];
% 
  tau1 = [11.97871 11.98071 11.98571 11.99071 11.99571  tau1];
% 
  % tau1 = 12.1 ;
%  
  N = 120      ;
%==================

h = 1/N      ;

% N = length(temp)-1;
% 
%  N = N/2;
% h = 1/N;
 

for p1 = 1:1%length(tau1)
    
    tau = tau1(p1);
    
    
    
    %-- 3 fold --
    
    % str0 = ['3fold_N' num2str(N) '_tau_' num2str(10*tau) '.txt'];

       
    %-- 5 fold ----    
    str0 = ['5pi_knot_N120_tau_' num2str(10*tau) '_1.txt'];
 
    if(tau>13.9)
            str0 = ['5pi_knot_N120_tau_' num2str(10*tau) '.txt'];
    end
    
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

% %====== Transforming the midline such that its midlplane is in x-y plane
% %and its geometric center is origin 
% 
% %--- 
%  nfold = 3;
%  ind = [1 1 + (N/nfold:N/nfold:((nfold-1)*N/nfold))];
% 
%   rx = rx -  mean(rx(ind,1)); 
%   ry = ry  -  mean(ry(ind,1));  
%   rz = rz -   mean(rz(ind,1));
%   
% %p1 = [rx(ind(1),1) ry(ind(1),1) rz(ind(1),1)];
% %p2 = [rx(ind(2),1) ry(ind(2),1) rz(ind(2),1)];
% %p3 = [rx(ind(3),1) ry(ind(3),1) rz(ind(3),1)];
% 
% 
% normal = cross( [rx(ind(1),1) ry(ind(1),1) rz(ind(1),1)] - [rx(ind(2),1) ry(ind(2),1) rz(ind(2),1)], [rx(ind(1),1) ry(ind(1),1) rz(ind(1),1)] - [rx(ind(3),1) ry(ind(3),1) rz(ind(3),1)]);
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
% for i = 1:length(rx) 
%     
% temp =  Ru*[rx(i,1);
%                    ry(i,1);
%                    rz(i,1)] ;
%  
% rx(i,1) = temp(1);
% ry(i,1) = temp(2);
% rz(i,1) = temp(3);
% 
% temp =  Ru*[bx(i,1);
%                    by(i,1);
%                    bz(i,1)] ;
%  
% bx(i,1) = temp(1);
% by(i,1) = temp(2);
% bz(i,1) = temp(3);
% 
% temp =  Ru*[tx(i,1);
%                    ty(i,1);
%                    tz(i,1)] ;
%  
% tx(i,1) = temp(1);
% ty(i,1) = temp(2);
% tz(i,1) = temp(3);
% 
% 
% temp =  Ru*[nx(i,1);
%                    ny(i,1);
%                    nz(i,1)] ;
%  
% nx(i,1) = temp(1);
% ny(i,1) = temp(2);
% nz(i,1) = temp(3);
%  
% end
% 
%  

%-- edge curve and its tangent 

rx0 = rx + w*bx;
ry0 = ry + w*by;
rz0 = rz + w*bz;



tx0 = tx -w*tau*nx;
ty0 = ty -w*tau*ny;
tz0 = tz -w*tau*nz;


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

 
 
% %===============================================
% figure(p1)
% plot3(rx,ry,rz,'color',col{2},'LineWidth',2)
% hold on 
% plot3(rx0,ry0,rz0,'color',col{4},'LineWidth',2)
% set(gca,'FontSize',25)
%  set(gcf,'color','w');
% box on
% grid on
% hold on 
% 
% %===============================================
 
% %===================================================
%n = n / (4 * midrad(pi, 1e-14))

n = n/(4*pi)*2;

Lk(p1) = n ;

figure(p1)
plotbrowser on
plot(rz)
legend(num2str(round(n)))
set(gca,'FontSize',25)
set(gcf,'color','w');
box on
grid on
hold on
end
% 
% %-- polynomial fit trial --
% 
% % Q = (1:N+1)';
% % sx = polyfit(Q,rx,5);
% % rxf = polyval(sx,Q);
% % 
% % sy = polyfit(Q,ry,5);
% % ryf = polyval(sy,Q);
% % 
% % sz = polyfit(Q,rz,5);
% % rzf = polyval(sz,Q);
% 
% 
% hold on 
% plot(rx,ry)
% hold on 
% plot(rx0,ry0)
% 
% 
% [indj2,indx]  = sort(indj);
% 
% indi2 = indi(indx);
% 
% tempi = indi2;
% tempj = indj2;
% 
% ne_tol = 4;
% 
% for i = 1:length(indj2)-1
%     for j = (i+1):length(indj2)
%         if(indj2(j)-indj2(i)<ne_tol )
%             
%             %-- comparing distance between the pair(indi2(i),indj2(i)) and (indi2(i),indj2(i))
%             
%             dsi = sqrt((rx(indi2(i),1)-rx0(indj2(i),1))^2 + (ry(indi2(i),1)-ry0(indj2(i),1))^2);
%             
%             dsj = sqrt((rx(indi2(j),1)-rx0(indj2(j),1))^2 + (ry(indi2(j),1)-ry0(indj2(j),1))^2);
%             
%             if(dsi<=dsj)
% %                 tempi(j) = tempi(i);
% %                 tempj(j) = tempj(i);
%                 
%                  tempi(i:j) = ones(j-i+1,1)*tempi(i);
%                 tempj(i:j) = ones(j-i+1,1)*tempj(i);
%             else
%                 tempi(i:j) = ones(j-i+1,1)*tempi(j);
%                 tempj(i:j) = ones(j-i+1,1)*tempj(j);
%             end
%         
%         else
%             break
%         end
%     end
% end
%  
% tempi = unique(tempi);
% tempj = unique(tempj);
% 
% for i = 1:length(tempi) 
%     
%     plot(rx(tempi(i),1),ry(tempi(i),1),'o','MarkerSize',10)
%     hold on 
%     plot(rx0(tempj(i),1),ry0(tempj(i),1),'o','MarkerSize',10)
%     hold on 
% end
% 
% hold on 
% plot(rx,ry)
% hold on 
% plot(rx0,ry0)


% 
% [indi,indj] = intersection_2curve(rx,ry,rx0,ry0);
% 
% 
% figure(p1)
% plot(rx,ry,'color',col{2},'LineWidth',2)
% hold on 
% plot(rx0,ry0,'color',col{4},'LineWidth',2)
% set(gca,'FontSize',25)
%  set(gcf,'color','w');
% box on
% grid on
% hold on


  
[indi,indj] = intersection_2curve(rx,ry,rz,rx0,ry0,rz0);


figure(p1)
plot3(rx,ry,rz,'color',col{2},'LineWidth',2)
hold on 
plot3(rx0,ry0,rz0,'color',col{4},'LineWidth',2)
hold on 
plot3(rx(indi,1),ry(indi,1),rz(indi,1),'o')
hold on 
plot3(rx0(indj,1),ry0(indj,1),rz0(indj,1),'o')

set(gca,'FontSize',25)
set(gcf,'color','w');
box on
grid on
hold on


% %--- ploting the intersection points 
% 
% plot(rx(indi),ry(indi),'o')
% hold on 
% 
% plot(rx0(indj),ry0(indj),'o')
% 
% [indi2,indj2] = intersection_1curve(rx,ry);
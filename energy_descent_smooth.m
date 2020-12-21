clear all
clc

global N N1 h tau sig fac branch   path  nstep strpath


smop = {'moving','lowess','lowess','sgolay','rlowess','rloess'};
opid  = 3;


col1 = [0,     0.976, 0.968]  ;      % -- color for mode n = 2
col2 = [0.831, 0.545, 0.247]  ;      % -- color for mode n = 3
col3 = [0.341, 0.505, 0.819]  ;      % -- color for mode n = 4
col4 = [0.705, 0.701, 0.070]  ;      % -- color for mode n = 5
col5 = [0.301, 0.811, 0.498]  ;      % -- color for mode n = 6
col6 = [1, 0.929, 0.278];
col = {col1,col2,col3,col4,col5,col6} ;
colw = [41,44,52]/255;



parameters5(); 

s = linspace(0,1,N+1);

N2 = 100;
s1 = linspace(0,1,N2+1);

%== prefactor for aspect ratio
 

 n = tau/2/pi;
         
for p1=  1:500%1:500%1:4000 %1:253%length(t1)

    str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*16)) '_step_' num2str(p1-1)   '.txt'];


    x = load([path str0]);

    bx(2:N,1)     = x(1:3:3*N1-2,1)    ;
    by(2:N,1)     = x(2:3:3*N1-1,1)    ;
    bz(2:N,1)     = x(3:3:3*N1,1)      ;

    bx(1,1) = 1;
    by(1,1) = 0;
    bz(1,1) = 0;

    bx(N+1,1) = -1;
    by(N+1,1) =  0;
    bz(N+1,1) =  0;

%     bx = smooth(bx,smop{opid});
%     by = smooth(by,smop{opid});
%     bz = smooth(bz,smop{opid});
% %     
    bx(1,1) = 1;
    by(1,1) = 0;
    bz(1,1) = 0;
    bx(N+1,1) = -1;
    by(N+1,1) =  0;
    bz(N+1,1) =  0;
    
    E(p1) = energy_b(bx,by,bz);
    
    lenE(p1,:) = sqrt(bx.^2+by.^2 + bz.^2)-1;
    
    
    bext(1:3,1) = -1*[bx(N,1);by(N,1);bz(N,1)];
  bext(4:6,1) = -1*[bx(2,1);by(2,1);bz(2,1)];
  
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
    
%     


% 
%     rx1 = interp1(s, rx, s1, 'spline', 'extrap');
%     rx1 = rx1'; 
%     
%     ry1 = interp1(s, ry, s1, 'spline', 'extrap');
%     ry1 = ry1'; 
%     rz1 = interp1(s, rz, s1, 'spline', 'extrap');
%     rz1 = rz1'; 
%     
%     

      rx1 = movmean(rx,4);
      ry1 = movmean(ry,4);
      rz1 = movmean(rz,4);
      
      rx1 = rx;
      ry1 = ry;
      rz1 = rz;
      
      rx1(end,1) = rx1(1,1);
      ry1(end,1) = ry1(1,1);
      rz1(end,1) = rz1(1,1);
      
      dl = sum(sqrt((rx1(1:end-1,1)-rx1(2:end,1)).^2 + (ry1(1:end-1,1)-ry1(2:end,1)).^2+(rz1(1:end-1,1)-rz1(2:end,1)).^2));
      
      rx1 = rx1/dl;
      ry1 = ry1/dl;
      rz1 = rz1/dl;
      
      
    rx1g  = [rx1(N2-1,1) rx1(N2,1) rx1' rx1(2,1) rx1(3,1)];
    ry1g  = [ry1(N2-1,1) ry1(N2,1) ry1' ry1(2,1) ry1(3,1)];
    rz1g  = [rz1(N2-1,1) rz1(N2,1) rz1' rz1(2,1) rz1(3,1)];
    
    
     
    
    
    h1 = 1/N2;
    
    ig = 3:N2+3;
    
    
    
    
    rx1p= (rx1g(ig+1)-rx1g(ig-1))/(2*h1);
    ry1p= (ry1g(ig+1)-ry1g(ig-1))/(2*h1);
    rz1p= (rz1g(ig+1)-rz1g(ig-1))/(2*h1);

    
    
    rx12p = (rx1g(ig+1) + rx1g(ig-1) -2*rx1g(ig))/h1^2;
    ry12p = (ry1g(ig+1) + ry1g(ig-1) -2*ry1g(ig))/h1^2;
    rz12p = (rz1g(ig+1) + rz1g(ig-1) -2*rz1g(ig))/h1^2;
    
    rx13p = (-1/2*rx1g(ig-2) + rx1g(ig-1) -rx1g(ig+1) +1/2*rx1g(ig+2))/h1^3;
    ry13p = (-1/2*ry1g(ig-2) + ry1g(ig-1) -ry1g(ig+1) +1/2*ry1g(ig+2))/h1^3;
    rz13p = (-1/2*rz1g(ig-2) + rz1g(ig-1) -rz1g(ig+1) +1/2*rz1g(ig+2))/h1^3;

    
    
  %  mod2 = (ty'.*rz12p - ry12p.*tz').^2 + (rx12p.*tz' - tx'.*rz12p).^2 + (tx'.*ry12p - rx12p.*ty').^2;
  %  tau2(p1,:) = (rx13p.*(ty'.*rz12p - ry12p.*tz') + ry13p.*(rx12p.*tz' - tx'.*rz12p) + rz13p.*(tx'.*ry12p - rx12p.*ty'))./mod2;
    
    
    mod2 = (ry1p.*rz12p - ry12p.*rz1p).^2 + (rx12p.*rz1p - rx1p.*rz12p).^2 + (rx1p.*ry12p - rx12p.*ry1p).^2;
    tau2(p1,:) = (rx13p.*(ry1p.*rz12p - ry12p.*rz1p) + ry13p.*(rx12p.*rz1p - rx1p.*rz12p) + rz13p.*(rx1p.*ry12p - rx12p.*ry1p))./mod2;
  
    
    
    
    
    
    E_r(p1) =  (asinh(n*pi*sig))/(8*pi^3*n^3)*sum(rx12p.^2+ry12p.^2+rz12p.^2);


end

% branch = 1;
%  path1 = [strpath num2str(branch) '/'];
% 
% strtau = ['tau_branch_' num2str(branch) '.txt'];
% 
% tau1 = load([path1 strtau]);
% 
% [val,ind] = min(abs(tau1-tau));
% 
% tau2 = tau1(ind);
% 
% 
% str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau2)) '.txt'];
% 
% 
% x = load([path1 str0]);
% 
% bx(2:N,1)     = x(1:3:3*N1-2,1)    ;
% by(2:N,1)     = x(2:3:3*N1-1,1)    ;
% bz(2:N,1)     = x(3:3:3*N1,1)      ;
% 
% bx(1,1) = 1;
% by(1,1) = 0;
% bz(1,1) = 0;
% 
% bx(N+1,1) = -1;
% by(N+1,1) =  0;
% bz(N+1,1) =  0;
% 
% %E(p1+1) = energy_b(bx,by,bz);

%
%==== E1 = 1.091317488162173e+00;  ====
%

figure(10)
plot(E,'-o')
hold on
plot(E(1),'-or')
hold on
plot(E_r  ,'-x');
%plot(p1+1,E(end),'-or')


plot([E(1:222) E_r(223:293) E(294:end)],'or')

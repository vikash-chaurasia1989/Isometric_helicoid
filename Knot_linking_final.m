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
  col1 = [0,     0.976, 0.968]    ;      % -- color for mode n = 2
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
%=============== 
% h = 1/N;
branch = 4;%
     path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];
    
   % path = [ '/Users/rtodres/Documents/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];

         
    
    strtau = ['tau_branch_' num2str(branch) '.txt'];
    
    tau1 = load([path strtau]);
    

for p1 = 618:618%1:length(tau1)
    
       N = 120      ; 
       N1 = N-1;
       h = 1/N      ;   
%==========================================================================
%==========================================================================

    bx = zeros(N+1,1);
    by = zeros(N+1,1);
    bz = zeros(N+1,1);
   
    rx = zeros(N+1,1);
    ry = zeros(N+1,1);
    rz = zeros(N+1,1);
    
    tx = zeros(N+1,1);
    ty = zeros(N+1,1);
    tz = zeros(N+1,1);
    
    rx0 = zeros(N+1,1);
    ry0 = zeros(N+1,1);
    rz0 = zeros(N+1,1);
    
    rx1 = zeros(N+1,1);
    ry1 = zeros(N+1,1);
    rz1 = zeros(N+1,1);
    
    
    %============================
    
    
    tau = tau1(p1);
    
%     
%     
%     %-- 3 fold --
%     
% %      str0 = ['3fold_N' num2str(N) '_tau_' num2str(10*tau) '.txt'];
% 
%        
%     %-- 5 fold ----    
%       str0 = ['5pi_knot_N120_tau_' num2str(10*tau) '_1.txt'];
% %  
%       if(tau>13.9)
%               str0 = ['5pi_knot_N120_tau_' num2str(10*tau) '.txt'];
%      end
% %     
%    % str0  = '3fold_N72_tau_170.txt';
%     
%     %-- 5pi --
%     %  str0 = ['5pi_knot_N120_tau_' num2str(10*tau1(p1)) '.txt'];
%     
%     
%     strb = ['b_' str0];
%     strr = ['r_' str0];
%     
%     strlm = ['lm_' str0];
%     strrho = ['rho_' str0];
%     struvw = ['uvw_' str0];
%     strbext = ['bext_' str0];
%     
%     
%   
%     %==========================================================================
%     N = 105;%140;
%     N1 = N-1;
%     h = 1/N;
%     
    
    %--- load b
 %   temp = load(strb);
 
 
   str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '.txt'];
        
        temp = load([path str0]);
        
        x = temp(1:3*N1,1);
        
        bx(2:N,1)     = x(1:3:3*N1-2,1)    ;
        by(2:N,1)     = x(2:3:3*N1-1,1)    ;
        bz(2:N,1)     = x(3:3:3*N1,1)      ;
        %--- boundary points --
        
        bx(1,1) = 1;
        by(1,1) = 0;
        bz(1,1) = 0;
        
        bx(N+1,1) = -1;
        by(N+1,1) =  0;
        bz(N+1,1) =  0;
 
 
% 
%  bx  = temp(1:N+1,1);
%  by = temp(1:N+1,2);
%  bz = temp(1:N+1,3);
 
    %-- load r --
    
    
%  s1 = (linspace(0,1,N+1))';
%  %
%  N = 1000;
%  N1 = N-1;
%  h = 1/N;
%  
%  s2 = (linspace(0,1,N+1))';
%  
%  bx = interp1q(s1,bx,s2)  ;
%  by = interp1q(s1,by,s2)  ;
%  bz = interp1q(s1,bz,s2)  ;
 
%  %--- rotating about x axis by pi/2 --
 R = [ cos(pi/2) -sin(pi/2);
         sin(pi/2)   cos(pi/2)];
    
    
%  for i = 1:N+1
%    
%      %-- Rotation about x axis by pi/2
%      temp = R*[by(i,1);
%                      bz(i,1)];
%                  
%       by(i,1) = temp(1);
%       bz(i,1) = temp(2);
%       
% %       %-- Rotation about z axis by pi/2 
%       temp = R'*[bx(i,1);
%                      by(i,1)];
%                  
%       bx(i,1) = temp(1);
%       by(i,1) = temp(2);
%  end
%  
 
%  dl = sqrt((bx(1:end-1,1)-bx(2:end,1)).^2 + (by(1:end-1,1)-by(2:end,1)).^2 +(bz(1:end-1,1)-bz(2:end,1)).^2);
% 
% tau = sum(dl);
%  
%  
% i = 1:N;
% 
% tx = by(i,1).*bz(i+1,1) - bz(i,1).*by(i+1,1);
% ty = bz(i,1).*bx(i+1,1) - bx(i,1).*bz(i+1,1);
% tz = bx(i,1).*by(i+1,1) - by(i,1).*bx(i+1,1);
% 
% 
% 
% tx(N+1,1) = tx(1,1);
% ty(N+1,1) = ty(1,1);
% tz(N+1,1) = tz(1,1);
% 
% tx = tx/tau/h;
% ty = ty/tau/h;
% tz = tz/tau/h;
% 
%  
% rx(1,1) = 0;
% ry(1,1) = 0;
% rz(1,1) = 0;
% 
% for i=1:N 
% 
%     rx(i+1,1) = rx(i,1) + h*tx(i,1);
%     ry(i+1,1) = ry(i,1) + h*ty(i,1);
%     rz(i+1,1) = rz(i,1) + h*tz(i,1);
% end


    
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
        
        
        
 nfold = 5;
 ind = [1 1 + (N/nfold:N/nfold:((nfold-1)*N/nfold))];

 
 
  rx = rx -  mean(rx(ind,1)); 
  ry = ry  -  mean(ry(ind,1));  
  rz = rz -   mean(rz(ind,1));
  
  

  
 

  

 w = 0.002;


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

gamma1(1,:) = rx' ; 
gamma1(2,:) = ry' ; 
gamma1(3,:) = rz' ; 

%-- creating array for linking number subroutine --
 
n= 0;
 

length_0 = N+1;
length_1 = N+1;

 

for i = 1:length_1
   
        for j = 1:length_1
        a = gamma1(:, j)                    + 0     - gamma0(:, i)            ;
        b = gamma1(:, j)                    - gamma0(:, mod(i, length_0) + 1) ;
        c = gamma1(:, mod(j, length_1) + 1) - gamma0(:, mod(i, length_0) + 1) ;
        d = gamma1(:, mod(j, length_1) + 1) + 0     - gamma0(:, i)            ;
        n = n + solid_angle_quadrilateral(a, b, c, d)                         ;
        
       
        end
               
        
end

n = n/(4*pi)*2;

Lk(p1) = n ;

 %==================================================
 %-------------- Intersection point using InterX subroutine
 %-- self intersection -- 
   P = InterX( [rx';ry'] );
 xin = P(1,:);   % these are interpolated points. May not exist in our data 
 yin = P(2,:);
 
 %-- extracting points in our data nearest to interpolated data xin and yin

 
 lenP = length(P(2,:));
 indi = zeros(1,lenP);
 indj = indi;
 
 num_tol = 3;   % -- avoiding neighborhood points from intersection 
 for m1 = 1:lenP
 
     
 temp = (rx-xin(m1)).^2 + (ry-yin(m1)).^2;
 
 [temp2,ind] =sort(temp);
 
 indi(m1) = ind(1) ;
 
 %-- point with which indi intersects 
 
 indj(m1) = ind(2);
 
 count = 0;
 
 while(abs(indj(m1)-indi(m1))<num_tol)
     count = count+1;
     
     indj(m1) = ind(2+count);
 end
 
% temp = find(abs(indj(1:m1-1)-indj(m1))<num_tol);
% 
% count = 0;
%      while(size(temp,2)~=0)
%            temp = find(abs(indj(1:m1-1)-indj(m1))<num_tol);
%     
%          count = count+1;
%      
%      indj(m1) = ind(2+count);
%      end
 
        

 
%      
 end
 
  %==== Eliminating the nearest neighbors ==
  
 
[indi,indjn] = sort(indi);
indj = indj(indjn);
 
%--- Eliminating nearest i
%  for i = 1:lenP-1
%      if(indi(i+1)-indi(i)<num_tol)
%          indi(i+1) = indi(i);
%          indj(i+1) = indj(i);
%      end
%      
%  end
%  
%  [indi,tempj] = unique(indi);
%  indj = indj(tempj);
         
%  %--- Eliminating nearest j
%  for i = 1:length(indi)-1
%      if(abs(indj(i+1)-indj(i))<num_tol)
%          indi(i+1) = indi(i);
%          indj(i+1) = indj(i);
%      end
%      
%  end
%  
% [indj,tempi] = unique(indj);
%  indi = indi(tempi);
 
%  tempi = indi;
%  tempj = indj;
%   for i = 1:length(indi) 
%      if(abs(indj(i)-indi(i))<num_tol)
%        tempi(i) = [];
%        tempj(i) = [];
%      end
%      
%   end
%    
%  indi = tempi;
%  indj = tempj;
 %===== plotting the sign of difference of the vertical height between z
 %coordinates of midline and edge 

 temp = zeros(length(indi),1);
 
  % temp = (rz(indi,1)-rz0(indj,1))./abs(rz(indi,1)-rz0(indj,1));
   
   %-- for self intersection -00-
  temp = (rz(indi,1)-rz(indj,1))./abs(rz(indi,1)-rz(indj,1));

 
count = 0;
    
%==== Calculating the knot linking number ===
link  = 0;
for i =1:length(indi)
    
    
    
    if(indi(i)==1)
    ar1 = [rx(N-1,1)     , ry(N-1,1)     , rz(N-1,1)];
    ar2 = [rx(indi(i)+2,1) , ry(indi(i)+2,1) , rz(indi(i)+2,1)]; 
   
    %-- if i index is one and j index is N+1 then replace N+1 with N-1 or 3
     if(indj(i)==N+1)
         indj(i) = 3;
     end
     
    elseif (indi(i)==2)
    ar1 = [rx(N,1)     , ry(N,1)     , rz(N,1)];
    ar2 = [rx(indi(i)+2,1) , ry(indi(i)+2,1) , rz(indi(i)+2,1)];
    elseif(indi(i)==N)
    ar1 = [rx(indi(i)-2,1) , ry(indi(i)-2,1) , rz(indi(i)-2,1)];  
    ar2 = [rx(2,1)     , ry(2,1)     , rz(2,1)];
    elseif(indi(i)==N+1)
    ar1 = [rx(indi(i)-2,1) , ry(indi(i)-2,1) , rz(indi(i)-2,1)];  
    ar2 = [rx(3,1)     , ry(3,1)     , rz(3,1)];
    else
    ar1 = [rx(indi(i)-2,1) , ry(indi(i)-2,1) , rz(indi(i)-2,1)]; 
    ar2 = [rx(indi(i)+2,1) , ry(indi(i)+2,1) , rz(indi(i)+2,1)];
    end

    %---
    %---- midline --
    if(indj(i)==1)
    ar3 = [rx(N-1,1)     , ry(N-1,1)     , rz(N-1,1)];
    ar4 = [rx(indj(i)+2,1) , ry(indj(i)+2,1) , rz(indj(i)+2,1)]; 
    elseif (indj(i)==2)
    ar3 = [rx(N,1)     , ry(N,1)     , rz(N,1)];
    ar4 = [rx(indj(i)+2,1) , ry(indj(i)+2,1) , rz(indj(i)+2,1)];
    elseif(indj(i)==N)
    ar3 = [rx(indj(i)-2,1) , ry(indj(i)-2,1) , rz(indj(i)-2,1)];  
    ar4 = [rx(2,1)     , ry(2,1)     , rz(2,1)];
    elseif(indj(i)==N+1)
    ar3 = [rx(indj(i)-2,1) , ry(indj(i)-2,1) , rz(indj(i)-2,1)];  
    ar4 = [rx(3,1)     , ry(3,1)     , rz(3,1)];
     %-- if i index is one and j index is N+1 then replace N+1 with N-1 or 3
     if(indj(i)==1)
         indj(i) = 3;
     end
     
    else
    ar3 = [rx(indj(i)-2,1) , ry(indj(i)-2,1) , rz(indj(i)-2,1)]; 
    ar4 = [rx(indj(i)+2,1) , ry(indj(i)+2,1) , rz(indj(i)+2,1)];
    end
    
    %====== Assigning coordinates to top and bottom arrow 
    
    if(rz(indi(i),1)>rz(indj(i),1))
        
        %--- Assigning ar1--ar2 as the top arrow 
        x1 = ar1(1);
        y1 = ar1(2);
        
        x2 = ar2(1);
        y2 = ar2(2);
        
        x3 = ar3(1);
        y3 = ar3(2);
        
        x4 = ar4(1);
        y4 = ar4(2);
    else
         %--- Assigning ar3--ar4 as the top arrow 
        x1 = ar3(1);
        y1 = ar3(2);
        
        x2 = ar4(1);
        y2 = ar4(2);
        
        x3 = ar1(1);
        y3 = ar1(2);
        
        x4 = ar2(1);
        y4 = ar2(2);
    end
    
    %========= p1 = (x1,y1) .  
        
    x4 = x4-x1;
    y4 = y4-y1;
    
    x3 = x3-x1;
    y3 = y3-y1;
    
    x2 = x2-x1;
    y2 = y2-y1;
    
    x1 = x1-x1;
    y1 = y1-y1; 
    
    %=== rotating vectors such that p1-p2 is positive x axis 
    
   
    
    if(x2>0)
        if(y2>0)
             th = acos(x2/(sqrt(x2^2+y2^2)));
        else
             th = 2*pi-acos(x2/(sqrt(x2^2+y2^2)));
        end
    else
        if(y2>0)
             th = pi/2+acos(-x2/(sqrt(x2^2+y2^2)));
        else
             th = pi + acos(-x2/(sqrt(x2^2+y2^2)));
        end
    end
        
    
    
    rot = [cos(th) sin(th);-sin(th) cos(th)];   
    
    
    tempcor = rot*[x2;y2];
    
    x2 = tempcor(1); y2 = tempcor(2);
    
    tempcor = rot*[x3;y3];
    
    x3 = tempcor(1);y3= tempcor(2);
    
    tempcor = rot*[x4;y4];
    
    x4 = tempcor(1);y4 = tempcor(2);
    
   
    if(y4>0)
        link(i) = 1;
    else 
        link(i) = -1;
    end
    %x1 = 0;y1 = 0;x2=0; y2 = 0; x3 = 0;y3 = 0; x4=0;y4=0; th = 0;
end

%=== This is to eliminate contribution from pseudo crossings 
 temp = abs(indi-indj);

for i = 1:length(temp)
    if(temp(i)<=num_tol)
        link(i) = 0;
    end
end
%================================
        
 Lk_knot(p1) = abs(sum(link)*2);
num_t(p1) = round(Lk(p1)-Lk_knot(p1));
 
% 
% if(tau>=12 && tau<12.8)
%     num_t(p1) =5;
%  
% elseif(tau > 14.5 && tau<16)
%     num_t(p1) = 5;
% elseif(tau > 16 && tau<17.5)
%     num_t(p1) = 7;
% elseif(tau>17.4)
%     num_t(p1) = 11;
% end
%  

%  if(tau>16.1)
%      Lk_knot(p1) = 0;
%  end
 
end
 
 
% 
% %==== Plotting the curve and arrows at the intersection points ===
% 
% figure(1)
% plotbrowser on 
% plot3(rx,ry,rz)
% hold on 
% %plot3(rx0,ry0,rz0)
% hold on 
% 
% for i =1:length(indi)
%     
%     
%     
%     if(indi(i)==1)
%     ar1 = [rx(N-1,1)     , ry(N-1,1)     , rz(N-1,1)];
%     ar2 = [rx(indi(i)+2,1) , ry(indi(i)+2,1) , rz(indi(i)+2,1)]; 
%     elseif (indi(i)==2)
%     ar1 = [rx(N,1)     , ry(N,1)     , rz(N,1)];
%     ar2 = [rx(indi(i)+2,1) , ry(indi(i)+2,1) , rz(indi(i)+2,1)];
%     elseif(indi(i)==N)
%     ar1 = [rx(indi(i)-2,1) , ry(indi(i)-2,1) , rz(indi(i)-2,1)];  
%     ar2 = [rx(2,1)     , ry(2,1)     , rz(2,1)];
%     elseif(indi(i)==N+1)
%     ar1 = [rx(indi(i)-2,1) , ry(indi(i)-2,1) , rz(indi(i)-2,1)];  
%     ar2 = [rx(3,1)     , ry(3,1)     , rz(3,1)];
%     else
%     ar1 = [rx(indi(i)-2,1) , ry(indi(i)-2,1) , rz(indi(i)-2,1)]; 
%     ar2 = [rx(indi(i)+2,1) , ry(indi(i)+2,1) , rz(indi(i)+2,1)];
%     end
% 
%         
%     
%    
%     
%     if(temp(i)>0)
%     mArrow3(ar1,ar2,'stemWidth',0.001,'tipWidth',0.002);
%     else
%     mArrow3(ar1,ar2,'stemWidth',0.001,'tipWidth',0.002,'facealpha',0.2);
%     end
%         
%     hold on 
%     plot3(rx(indi(i),1),ry(indi(i),1),rz(indi(i),1),'o')
%     hold on 
%     hold on 
%     
%     
%     %---
%     %---- midline --
%     if(indj(i)==1)
%     ar1 = [rx(N-1,1)     , ry(N-1,1)     , rz(N-1,1)];
%     ar2 = [rx(indj(i)+2,1) , ry(indj(i)+2,1) , rz(indj(i)+2,1)]; 
%     elseif (indj(i)==2)
%     ar1 = [rx(N,1)     , ry(N,1)     , rz(N,1)];
%     ar2 = [rx(indj(i)+2,1) , ry(indj(i)+2,1) , rz(indj(i)+2,1)];
%     elseif(indj(i)==N)
%     ar1 = [rx(indj(i)-2,1) , ry(indj(i)-2,1) , rz(indj(i)-2,1)];  
%     ar2 = [rx(2,1)     , ry(2,1)     , rz(2,1)];
%     elseif(indj(i)==N+1)
%     ar1 = [rx(indj(i)-2,1) , ry(indj(i)-2,1) , rz(indj(i)-2,1)];  
%     ar2 = [rx(3,1)     , ry(3,1)     , rz(3,1)];
%     else
%     ar1 = [rx(indj(i)-2,1) , ry(indj(i)-2,1) , rz(indj(i)-2,1)]; 
%     ar2 = [rx(indj(i)+2,1) , ry(indj(i)+2,1) , rz(indj(i)+2,1)];
%     end
%     
%     
%     if(temp(i)>0)
%       mArrow3(ar1,ar2,'stemWidth',0.001,'tipWidth',0.002,'facealpha',0.2);
%     else
%         
%     mArrow3(ar1,ar2,'stemWidth',0.001,'tipWidth',0.002);
% 
%     end  
%     hold on 
%     plot3(rx(indj(i),1),ry(indj(i),1),rz(indj(i),1),'x')
%     hold on 
%     
% %     %-- edge --
% %     ar1 = [rx0(indj(i)-4,1)     , ry0(indj(i)-4,1)     , rz0(indj(i)-4,1)];
% %     ar2 = [rx0(indj(i)+1,1) , ry0(indj(i)+1,1) , rz0(indj(i)+1,1)];
% %     
% %     mArrow3(ar1,ar2,'stemWidth',0.002,'facealpha',0.5);
% %     hold on 
% %     plot3(rx0(indj(i),1),ry0(indj(i),1),rz(indj(i),1))
% %     hold on 
% %     hold on 
% end
% 
%     
%   
% 
% 
%     
%  
%     
% 
% 

%--- Figure 
% num_t = round(num_t);
% ind = find(num_t==-1);
% num_t(ind) = 3;
% figure(20)
% plotbrowser on 
% plot(tau1/(2*pi),round(num_t),'--o')
% set(gca,'FontSize',25)
%  set(gcf,'color','w');
% box on
% grid on
% hold on 

plot(num_t,'--o')


%=== saving midline in xyz format ====

str0 = 'midline_knot.xyz';
%
fileID = fopen([path str0],'w');



for i = 1:N+1
    fprintf(fileID,'%d %30.16E %30.16E %30.16E \r\n'  ,[i rx(i,1) ry(i,1) rz(i,1)] );
end
fclose(fileID);



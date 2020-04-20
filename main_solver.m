clear all
clc

  global    lm  N1 p1 count tau  rho f_b  h tx ty tz  u v w p2

  global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p  N   len rx ry rz Lk err qd   E bext

  
  
  %===== add path to data file =====
  addpath('/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_3pi')
  addpath('/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_5pi')
  addpath('/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_7pi')


  format longE
  tic

%******************************** Main file for two chain interaction *********
 
tau_3 = [8.1:.1:9.5 10:1:11.9 11.97871 11.98071 11.98571 11.99071 11.99571  12:.1:12.5 12.55 12.552 12.554 12.555 12.556 12.56 12.57 12.6:.1:13.2 13.21 13.215 ...
             13.2165 13.21675 13.216865 13.21686525 13.2168653 13.2168654 13.21686545 13.21686546 13.21686547 13.2168654725    ] ;
tau_5 = [13.216865473 13.216865474355 13.216865474 13.2168655 13.21686575 13.216866 13.21687 13.2169 13.217 13.2175 13.22 13.23];
tau_7 = [13.3:.1:13.9 14.5:.5:21];%tau1 = [8.1:.1:9.5 10:1:11.9 11.97871 11.98071 11.98571 11.99071 11.99571  12:.1:12.5 12.55 12.552 12.554 12.555 12.556 12.56 12.57 12.6:.1:13.2 13.21 13.1:.1:13.9 14:.5:21];
tau1 =  [tau_3 tau_5 tau_7];


N = 49;
h = 1/N;

N1 = N-1;
%=======
  
  %
  
  %====== Tau for 5pi curves =====
%  tau_5 = [11.97871 11.98071 11.98571 11.99071 11.99571 12:.1:12.4 12.45 12.465 12.5:.1:15.8 15.85   16.06];
%  tau_7 = [16.1:.1:17.4 17.475 17.478 17.4783 1.7479604e+01];
%  tau_9 = [17.482 17.4825 17.48275 17.4829256];
%  tau_11 = [ 17.482926 17.4844 17.485 17.5:.1:19  19.5:.5:25];
%  
%  % % 
%   tau1 = [tau_5 tau_7 tau_9 tau_11];
%===========================
  
  %===   %indx =  11    12    30    38    51    54    57    66    69    70    75    76    79    80   100


 
  p2 = 56;
  tau = 1.2400000000000e+01 ;%tau_tau1(p2);
  
  
  %=== tau for 7fold ===
  
  tau1 = [ 1.545269410042355e+01 16];
  
  tau = 15.8;% 1.539778078168040e+01;%1.538970616264540e+01
  
 %tau1(1);
 
  %=== initial guess for 3pi twist ==
  %var_initial = initial_guess_3pi();
  %=====================
  
  nfold = 7;
  
  if(nfold==3)
      var_initial = initial_guess_3pi();
  elseif(nfold==5)
    %=== initial guess for 5pi twist ==
  var_initial = initial_guess_5pi();
  %=====================
  else
  
   %==== initial guess for 7pi twist ===
    var_initial = initial_guess_7pi();
  end
  
  %=== initial guess with refined data ===
  % var_initial = initial_guess_refine();
 
 %=== initial guess with coarse data ====
   %  var_initial = initial_guess_coarse();
      
      


 count = 1;
%--------------------------- Solver ---------------------------------------

 %options.Algorithm   = 'levenberg-marquardt'  ;    %'trust-region-reflective' ;
   options             = optimset('Display','iter', 'Algorithm','levenberg-marquardt','Jacobian','on', 'TOlFun',10^(-32),'TOlX',10^(-32),'MaxFunEvals',695000) ;
  options.MaxIter     = 100000 ;
    %     [x,fval,exitflag,output,qd1] =  fsolve(@fun_curve2 ,var_initial,options)          ;
      [x,fval,exitflag,output,qd1] =  fsolve(@fun_jacobian ,var_initial,options)          ;

              %  options =  optimset('TolX', 10^(-16));
                    % f = newtonraphson(@fun_curve4,var_initial)
                  %     f = newtonraphson4(var_initial)  ;

%--------------------------------------------------------------------------

%------------------------- Plotting  --------------------------------------
 rx(1,1) = 0; ry(1,1) = 0; rz(1,1) = 0;
 Lk = 3;
%---- total length and total arc length ---

dal = 2*asin(len/2);

alpha = sum(dal);

L = sum(len);




toc

%----------------- Post processing -----
h =1/N;% mean(len);
%---- Tangent ti = bi \times bi+1

i = 1:N+1;


% tx(i,1) = by(i,1).*bz(i+1,1) - bz(i,1).*by(i+1,1);
% ty(i,1) = bz(i,1).*bx(i+1,1) - bx(i,1).*bz(i+1,1);
% tz(i,1) = bx(i,1).*by(i+1,1) - by(i,1).*bx(i+1,1);
%
tx(i,1) = by(i,1).*bzp(i,1) - bz(i).*byp(i,1);
ty(i,1) = bz(i,1).*bxp(i,1) - bx(i).*bzp(i,1);
tz(i,1) = bx(i,1).*byp(i,1) - by(i).*bxp(i,1);

%
%   tx(N+1,1) = tx(1,1);
%   ty(N+1,1) = ty(1,1);
%   tz(N+1,1) = tz(1,1);
% %
% %
 tx = tx/tau;
 ty = ty/tau;
 tz = tz/tau;

%--------------------

% %--- position vector using integration of tangent
% % initialization
%h=1;
rx(1) = 0; ry(1) = 0; rz(1) = 0;

% i = 1 ;
%
%     rx(i+1) = rx(i) + h*tx(i);
%     ry(i+1) = ry(i) + h*ty(i);
%     rz(i+1) = rz(i) + h*tz(i);
%
% for i = 2:N-1
%     rx(i+1) = rx(i-1) + 2*h*tx(i);
%     ry(i+1) = ry(i-1) + 2*h*ty(i);
%     rz(i+1) = rz(i-1) + 2*h*tz(i);
% end
%


for i = 1:N

     rx(i+1) = rx(i) + h*tx(i);
     ry(i+1) = ry(i) + h*ty(i);
     rz(i+1) = rz(i) + h*tz(i);

end


figure(2)
plotbrowser on
plot3(bx,by,bz,'LineWidth',4)
title('binormal')
set(gca,'FontSize',25,'LineWidth',1)
box on
grid on

figure(3)
plotbrowser on
plot3(rx,ry,rz,'LineWidth',4)
title('Midline')
set(gca,'FontSize',25,'LineWidth',1)
box on
grid on

figure(4)
plotbrowser on
plot3(tx,ty,tz,'LineWidth',4)
title('Tangent')
set(gca,'FontSize',25,'LineWidth',1)
box on
grid on

%   str = ['b_N20_tau_' num2str(p1) '.txt'];
%
%   fileID = fopen(str,'w');
%   fprintf(fileID,'%30.16E %30.16E %30.16E \r\n',[bx'    ;by'    ; bz'    ] );
%   fclose(fileID);

  %===== curvature energy ==

  Ebend  = 1/tau^4*(h*sum(bx2p.^2 + by2p.^2 + bz2p.^2) - tau^2);
%  end
 
%---  tau                   Ebend

%---- 8.094              1.402667800494005e+02
%---- 8.095              1.403194427103710e+02
%---- 8.096              1.403983444455394e+02
%-----8.100              1.412415425473598e+02

%---

% %--- rotation about z axis
%-- nf = 2*5;
% 
%   nf = 2*7;
% for i = 1:N+1
%     
%     bx1(i,1) = cos(pi/nf)*bx(i,1) -sin(pi/nf)*by(i,1);
%     by1(i,1) = sin(pi/nf)*bx(i,1) +cos(pi/nf)*by(i,1);
% end
% 
% bx = bx1;
% by = by1;
% 
% % %-- rotation about x axis --
% 
% th = pi/2-atan(by(2,1)/bz(2,1));
% %th = pi;
% 
% for i = 1:N+1
%     
%     bz1(i,1) = cos(th)*bz(i,1) -sin(th)*by(i,1);
%     by1(i,1) = sin(th)*bz(i,1) +cos(th)*by(i,1);
% end
% 
% by = by1;
% bz = bz1;
% 
% th = asin(bz(2,1));


%========================================================
%==========  Linking number ===================================
%========================================================
 %---------- writhe calculation 
 wd =  1 ; % width of the edge
 
 rx = rx';
 ry = ry';
 rz = rz';
 
 
for i = 1:N     
     % 
     j = 1:N+1     ;             % j array index
     j = j(j~=i);   % indexes not containing i
      
         %--------------------- 1. Curve r(s)  -------------------------------
         
      rij  =  ((rx(j,1)-rx(i,1)).^2 + (ry(j,1)-ry(i,1)).^2 +(rz(j,1)-rz(i,1)).^2).^.5 ;
       
      %-- for self linking number
      tjix = ty(j,1).*tz(i,1) - tz(j,1).*ty(i,1);
      tjiy = tz(j,1).*tx(i,1) - tx(j,1).*tz(i,1);
      tjiz = tx(j,1).*ty(i,1) - ty(j,1).*tx(i,1);
         
      
      temp =   (tjix.*(rx(j,1)-rx(i,1)) + tjiy.*(ry(j,1)-ry(i,1)) + tjiz.*(rz(j,1)-rz(i,1)))./rij.^3 ;
        
      ker(i,1) = trapz(temp);
      
      
      
end

Wr  =  h^2*trapz(ker)/(4*pi);
Lk  = 2*(tau/(2*pi) + Wr )  ;

% tx = tx';
% ty = ty';
% tz = tz';
% 
% for i = 1:N 
%      % 
%      j = 1:N    ;             % j array index
%      j = j(j~=i);   % indexes not containing i
%       
%          %--------------------- 1. Curve r(s)  -------------------------------
%          
%       rij  =  ((rx(i)-rx(j)).^2 + (ry(i)-ry(j)).^2 +(rz(i)-rz(j)).^2).^.5 ;
%        
%       %-- for self linking number
%       tjix = ty(j).*tz(i) - tz(j).*ty(i);
%       tjiy = tz(j).*tx(i) - tx(j).*tz(i);
%       tjiz = tx(j).*ty(i) - ty(j).*tx(i);
%          
%       
%       temp =   (tjix.*(rx(j)-rx(i)) + tjiy.*(ry(j)-ry(i)) + tjiz.*(rz(j)-rz(i)))./rij.^3 ;
%         
%       ker(i) = trapz(temp);
%       
%       
%       
% end
% 
% Wr  = h^2*trapz(ker)/(4*pi);
% Lk  = 2*(tau/(2*pi) + Wr )  ;

%========================================================
%========================================================

%===== saving data ====
  %str0 = ['3fold_N' num2str(N) '_tau_' num2str(10*tau) '_1.txt'];
%        str0 = ['3fold_N' num2str(N) '_tau_' num2str(10*tau) '.txt'];
% 
%      if(tau>11.9)
%           str0 = ['3fold_N' num2str(N) '_tau_' num2str(10*tau) '_1.txt'];
%      end
     
%       %-- 5 fold ----    
%       str0 = ['5pi_knot_N120_tau_' num2str(10*tau) '_1.txt'];
% %  
%       if(tau>13.9)
%               str0 = ['5pi_knot_N120_tau_' num2str(10*tau) '.txt'];
%      end 
%   
%  str0 = ['5pi_knot_N120_tau_' num2str(10*tau) '.txt'];

%=== 7pi 
 %str0 = ['7pi_knot_' num2str(N) '_tau_' num2str(10*tau) '.txt'];

 
 %==== New strings for saving data files =========
 
 %===== 3pi ======
 
%  str0 = ['3fold_N72_tau_' num2str(10^10*tau) '.txt'];
%   path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_3pi/';
    
%=================================






%====== 5pi ====
% 
% if(nfold==3)
%     str0 = ['3fold_N' num2str(N) '_tau_' num2str(10^10*tau) '.txt'];
%     path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_3pi/';
% elseif(nfold==5)
%     str0 = ['5pi_knot_N' num2str(N) '_tau_' num2str(10^5*tau) '.txt'];
%    path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_5pi/';
% %  
% else
%      str0 = ['7pi_knot_N' num2str(N) '_tau_' num2str(10^10*tau) '.txt'];
%    path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_7pi/';
% end
% 
% 
% 
% strb = ['b_' str0];
% strr = ['r_' str0];
% 
% strlm = ['lm_' str0];
% strrho = ['rho_' str0];
% struvw = ['uvw_' str0];
% strbext = ['bext_' str0];
% 
% fileID = fopen([path strb],'w');
% fprintf(fileID,'%30.16E %30.16E %30.16E \r\n',[bx'    ;by'    ; bz'    ] );
% fclose(fileID);
% 
%  fileID = fopen([path strr],'w');
% fprintf(fileID,'%30.16E %30.16E %30.16E \r\n',[rx'    ;ry'    ; rz'   ] );
% fclose(fileID);
% 
% 
% fileID = fopen([path strlm],'w');
% fprintf(fileID,'%30.16E   \r\n',lm );
% fclose(fileID);
% 
% fileID = fopen([path strrho],'w');
% fprintf(fileID,'%30.16E   \r\n',rho ); 
% fclose(fileID);
% 
% fileID = fopen([path strbext],'w');
% fprintf(fileID,'%30.16E   \r\n',bext );
% fclose(fileID);
% 
% fileID = fopen([path struvw],'w');
% fprintf(fileID,'%30.16E   \r\n',[u v w]' );
% fclose(fileID);
% 

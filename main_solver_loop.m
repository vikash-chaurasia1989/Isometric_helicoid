clear all
clc

  global    lm  N1 p1 count tau  rho f_b  h tx ty tz  u v w tau1 p1 fac

  global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p  N   len rx ry rz Lk err   nfold branch sig w1 w2 w3

  format longE
  tic
 

  col1 = [0,     0.976, 0.968]  ;      % -- color for mode n = 2
  col2 = [0.831, 0.545, 0.247]  ;      % -- color for mode n = 3
  col3 = [0.341, 0.505, 0.819]  ;      % -- color for mode n = 4
  col4 = [0.705, 0.701, 0.070]  ;      % -- color for mode n = 5
  col5 = [0.301, 0.811, 0.498]  ;      % -- color for mode n = 6
  col6 = [1, 0.929, 0.278];
  col = {col1,col2,col3,col4,col5,col6} ;
  colw = [41,44,52]/255;
 %===== branch 1   
  
tau_3 = [8.09575 8.0958  8.09585 8.096 8.097 8.09975 8.1:.1:9.5 9.6:.1:11.9 12:.1:12.5 12.55 12.552 12.554 12.555 12.556 12.56 12.57 12.6:.1:13.2 13.21 13.215 ...
             13.2165 13.21675 13.216865 13.21686525 13.2168653 13.2168654 13.21686545 13.21686546 13.21686547 13.2168654725    ] ;          
tau_5 = [13.21686547250001  13.21686547255 13.2168654726 13.21686547275 13.216865473 13.216865474355 13.216865474 13.2168655 13.21686575 13.216866...
             13.21687 13.2169 13.217 13.2175 13.22 13.2225 13.2235 13.2236 13.223675 13.22375 13.2238];
         
tau_7 = [13.224  13.225 13.23 13.3:.1:13.9 14.5:.5:21 21.5:.5:45];



%tau1 = [8.1:.1:9.5 10:1:11.9 11.97871 11.98071 11.98571 11.99071 11.99571  12:.1:12.5 12.55 12.552 12.554 12.555 12.556 12.56 12.57 12.6:.1:13.2 13.21 13.1:.1:13.9 14:.5:21];
tau1 =      [tau_3 tau_5 tau_7];


N = 72;
  N1 = N-1;
  h = 1/N;
  nfold=3;
  
  
%   
% % %   %=========== branch 2 
% %   
%    tau1 = [ 11.984235 11.98425  11.9845 11.985 12:.1:25 25.5:.5:34 34.1 34.15 34.155 34.157  34.158 34.159];
%    
%    
%   N = 120;
%   N1 = N-1;
%   h = 1/N;  
%     
% %   
%   nfold=5;
% % %===== branch 3 
% % 
%  tau1 = [15.45:.1:25 25.1:.1:40];
%    
%   N = 84;
%   N1 = N-1;
%   h = 1/N;  
%   
%  nfold=7;
    % tau1 = 21;
    
 branch =  2;   
    
 if(branch==1)   % 8.0957344694550208E+00   
%          tau_3 = [8.09575 8.0958  8.09585 8.096 8.097 8.09975 8.1:.1:9.5 9.6:.1:11.9 12:.1:12.5 12.55 12.552 12.554 12.555 12.556 12.56 12.57 12.6:.1:13.2 13.21 13.215 ...
%          13.2165 13.21675 13.216865 13.21686525 13.2168653 13.2168654 13.21686545 13.21686546 13.21686547 13.2168654725    ] ;
%      tau_5 = [13.21686547250001  13.21686547255 13.2168654726 13.21686547275 13.216865473 13.216865474355 13.216865474 13.2168655 13.21686575 13.216866...
%          13.21687 13.2169 13.217 13.2175 13.22 13.2225 13.2235 13.2236 13.223675 13.22375 13.2238];
%      
%      tau_7 = [13.224  13.225 13.23 13.3:.1:13.9 14.5:.5:21 21.5:.5:45];
         
       tau_3 = [ 8.0957344694550208E+00  8.09574 8.095745 8.09575 8.0958  8.09585 8.096 8.097 8.09975 8.1:.05:8.5 8.6:.1:9.5 9.6:.1:11.9 12:.1:12.5 12.55 12.552 12.554 12.555 12.556 12.56 12.57 12.6:.1:13.2...
                  13.3 13.31 13.32  13.33 13.335 13.3375 13.33755];
       tau_5 =   [ 13.337575 13.3376 13.33775];     
     
       tau_7 = [ 13.338 13.3385 13.34 13.35  13.36 13.37 13.38 13.39 13.395  13.4:.1:13.9 14.5:.5:18 18.1:.1:22.9 23:.5:27 27.01:.01:27.05 27.0525 27.0535 27.05355 27.053575];  
       tau_9 = [27.05375 27.054 27.0545  27.055 27.06:.01:27.13 27.1325 27.132575 27.133];
       tau_11 = [27.1335 27.135 27.14 27.15 27.2:.1:27.4 27.5:.5:33.5 33.6:.1:33.9 33.95 33.96 33.965 33.966 33.9665 33.96675 33.96676 33.96677 ];
%      
      
     
     tau1 =      [tau_3 tau_5 tau_7 tau_9 tau_11];
     
     
     N = 72;
     N1 = N-1;
     h = 1/N;
     nfold=3;
     
 elseif(branch==2) % 1.1984230303215979E+01
     %================ Branch 2 =====================
     
    % tau1 = [1.1984233303215979E+01   11.984235 11.98425  11.9845 11.985 12:.1:25 25.5:.5:34 34.1 34.15 34.155 34.157  34.158 34.159];
     tau1 = [1.1984233303215979E+01   11.984235 11.98425  11.9845 11.985 12:.1:16.5 16.55 16.56 16.561 16.5625 16.565 16.57 16.575...
              16.6: 16.7 16.8 16.81 16.815 16.8175 16.8176 16.81775 16.818 16.82 16.825 16.85 16.86 16.865 16.8675 16.8685 16.8695 16.869575 16.8696 16.86975 16.8699 16.87 16.875 16.9:.1:25 25.5:.5:34 34.1 34.15 34.155 34.157  34.158 34.159];

     
     N = 120;
     N1 = N-1;
     h = 1/N;
     nfold=5;
     
    % =============== Branch 3 =====================
     
 else
     tau1 = [15.45:.1:20.95 20.96:.01:21.25  21.35:.1:25 25.1:.1:40];
     
     N = 84;
     N1 = N-1;
     h = 1/N;
     
     nfold=7;
 end
%==========================================================================  
  sig = .06;
 
  
   for p1  = 1:3%48:length(tau1)
       sv = 1;
       tau =   tau1(p1);
       
      % n = tau/2/pi;
       
     %  fac = asinh(n*pi*sig)/(4*pi^3*n^3);
 % for p1 = 1:length(tau_data)

 
 % var_initial = initial_guess_loop2();
  
  var_initial = initial_guess_loop8();


    %var_initial = initial_guess_refine_loop();

  %=====================
  
  %=== initial guess with refined data ===
  %  var_initial = initial_guess_refine();
 
 %=== initial guess with coarse data ====
 %  var_initial = initial_guess_coarse();
      
      


 count = 1;
%--------------------------- Solver ---------------------------------------

 %options.Algorithm   = 'levenberg-marquardt'  ;    %'trust-region-reflective' ;
   options             = optimset('Display','iter', 'Algorithm','levenberg-marquardt','Jacobian','on', 'TOlFun',10^(-20),'TOlX',10^(-20),'MaxFunEvals',69500  ) ;
  options.MaxIter     = 50000  ;
    %     [x,fval,exitflag,output,qd1] =  fsolve(@fun_curve2 ,var_initial,options)          ;
        [x,fval,exitflag,output,qd1] =  fsolve(@fun_jacobian8     ,var_initial,options)          ;
        %      [x,fval,exitflag,output,qd1] =  fsolve(@fun_jacobian7     ,var_initial,options)          ;

      % x =  lsqnonlin(@fun_jacobian6 ,var_initial )          ;
  %fun_jacobian6(var_initial);
              %  options =  optimset('TolX', 10^(-16));
                    % f = newtonraphson(@fun_curve4,var_initial)
                  %     f = newtonraphson4(var_initial)  ;
%-- saving convergence error 

err1(p1) = err;
%--------------------------------------------------------------------------

%------------------------- Plotting  --------------------------------------
  rx(1) = 0; ry(1) = 0; rz(1) = 0;
Lk = 3;
%---- total length and total arc length ---

dal = 2*asin(len/2);

alpha = sum(dal);

L = sum(len);




toc

%----------------- Post processing -----
h =1/N;% mean(len);
%---- Tangent ti = bi \times bi+1

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

     rx(i+1) = rx(i) + h*tx(i+1);
     ry(i+1) = ry(i) + h*ty(i+1);
     rz(i+1) = rz(i) + h*tz(i+1);

end

% for i = 2:N+1
% 
%      rx(i ) = rx(i-1) + h*tx(i);
%      ry(i ) = ry(i-1) + h*ty(i);
%      rz(i ) = rz(i-1) + h*tz(i);
% 
% end

 
% 
% figure(2)
% plotbrowser on
% title('binormal')
% plot_sphere()
% hold on
% plot3(bx,by,bz,'color',col{1},'LineWidth',2)
% hold on
% set(gca, 'DataAspectRatio',[1,1,1]);
% daspect([1,1,1]);
% axis off
% set(gcf,'color',colw);
% axis vis3d
% 
% 
% 
% figure(3)
% plotbrowser on
% title('midline')
% hold on
% plot3(rx,ry,rz,'color',col{1},'LineWidth',2)
% hold on
% set(gca, 'DataAspectRatio',[1,1,1]);
% daspect([1,1,1]);
% axis off
% set(gcf,'color',colw);
% axis vis3d

% 
% figure(4)
% plotbrowser on
% title('tangent')
% hold on
% plot3(tx,ty,tz,'color',col{1},'LineWidth',2)
% hold on
% set(gca, 'DataAspectRatio',[1,1,1]);
% daspect([1,1,1]);
% axis off
% set(gcf,'color',colw);
% axis vis3d


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
% for i = 1:N+1
%     
%     bx1(i,1) = cos(pi/10)*bx(i,1) -sin(pi/10)*by(i,1);
%     by1(i,1) = sin(pi/10)*bx(i,1) +cos(pi/10)*by(i,1);
% end
% 
% bx = bx1;
% by = by1;
% 
% %-- rotation about x axis --
% 
% th = pi/2-atan(by(2,1)/bz(2,1));
% th = pi;
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

Lk(p1) = Lk;

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

 rx = rx';
 ry = ry';
 rz = rz';
 
 

%====== 5pi ====
% 


if(sv==1)
    
if(nfold==3)
    str0 = ['3fold_N' num2str(N) '_tau_' num2str(10^10*tau) '.txt'];
  % str0 = ['3fold_N' num2str(N) '_tau_' num2str(10^10*tau) '_sig_' num2str(sig) '.txt'];

    path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_3pi_4/';
elseif(nfold==5)
    str0 =['5pi_knot_N' num2str(N) '_tau_' num2str(10^10*tau) '.txt'];
    path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_5pi_4/';
%  
else
     str0 = ['7pi_knot_N' num2str(N) '_tau_' num2str(10^10*tau) '.txt'];
   path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_7pi_4/';
end



strb = ['b_' str0];
strr = ['r_' str0];

strlm = ['lm_' str0];
strrho = ['rho_' str0];
struvw = ['uvw_' str0];


 % 
fileID = fopen([path strb],'w');
fprintf(fileID,'%30.16E %30.16E %30.16E \r\n',[bx'    ;by'    ; bz'    ] );
fclose(fileID);

 

fileID = fopen([path strlm],'w');
fprintf(fileID,'%30.16E   \r\n',lm );
fclose(fileID);

fileID = fopen([path strrho],'w');
fprintf(fileID,'%30.16E   \r\n',rho ); 
fclose(fileID);

 

fileID = fopen([path struvw],'w');
fprintf(fileID,'%30.16E   \r\n',[u v w]' );
fclose(fileID);


 








end
   end
%=== saving the torsion values for a given branch 

if(sv==1)
strtau = ['tau_branch_' num2str(branch) '.txt'];

fileID = fopen([path strtau],'w');
fprintf(fileID,'%30.16E   \r\n',tau1' ); 
fclose(fileID);
end

%---
clear all
clc
 

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
col6 = [0.341, 0.505, 0.819]/2  ;         % -- color for mode n = 6
col7 = [0.831, 0.545, 0.247]/2  ;      % -- color for mode n = 3


col = {col1,col2,col3,col4,col5,col6,col7} ;
colw = [41,44,52]/255;

%temp = load('05_fold.txt');
% == 253:309


p3 = 4;
tau0 = 32;%
for branch =p3:p3
    path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];
    
   % path = [ '/Users/rtodres/Documents/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];

         
    
    strtau = ['tau_branch_' num2str(branch) '.txt'];
    
    tau1 = load([path strtau]);
    
    %==========================================================================
    if (branch==1||branch==2||branch==3)
            N=105;
        end
        
        if(branch==4||branch==9)
            N = 180; %
        end
        if(branch==5)
            N = 140;
        end
        
        if(branch== 6)
            N = 90;
        end
        
        if(branch== 7)
            N = 150;
        end
        
        
        if(branch== 8)
            N = 150;
        end
        
        if(branch== 10)
            N = 220;
        end
        
        if(branch==11||branch==14||branch==15)
            N=170;
        end
       
        if(branch== 12)
            N = 190;
        end
        
        
        if(branch==13)
             N=210;
            
            %N=300;
            
            
        end
        N1 = N-1;
         h = 1/N;
    
    
    
    Lk = zeros(1,length(tau1));
    st = zeros(1,length(tau1));
    
    
    
    rx = zeros(N+1,1);
    ry = zeros(N+1,1);
    rz = zeros(N+1,1);
    
    
    bx = zeros(N+1,1);
    by = zeros(N+1,1);
    bz = zeros(N+1,1);
    
    
    %=== loading stability results ===
    
    
    strst = ['stability_branch_' num2str(branch) '_.txt'];
    
    %st = load([path strst]);
    [temp,p2] = min(abs(tau0-tau1));  % index of the saved data nearest to the given tau value

     tau  = tau1(p2);


   % p2 = 64;
    for p1 =  p2:p2% 1:length(tau1)
        
        sv = 0;
        tau =  tau1(p1);
        
     %   path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];
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
        
        
        
         
        
        nx = by.*tz - bz.*ty;
        ny = bx.*tz - bz.*tx;
        nz = bx.*ty - by.*tx;
        
        
        %
        
        
     %     plot3(rx,ry,rz)
%     hold on 
    %---- create edge of the binormal scroll using b and r
    wd = .1;
    
  %  wd = .1;
    v = linspace(-wd,wd,10);
    
 %--- edge   ----
 
    
    rx0 = rx - v(1)*bx;
    ry0 = ry - v(1)*by;
    rz0 = rz - v(1)*bz;
    
    tx0 = tx - v(1)*bxp;
    ty0 = ty - v(1)*byp;
    tz0 = tz - v(1)*bzp;
   
    
%==========================================================================

%-------- self linking number ---
 
 
 
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
Lk(p1)  =  round_odd(2*(tau/(2*pi) + Wr ))  ;



%===========================  Edge linking ===============================
 
for i = 1:N     
     % 
     j = 1:N+1     ;             % j array index
   %  j = j(j~=i);   % indexes not containing i
      
         %--------------------- 1. Curve r(s)  -------------------------------
         
      rij  =  ((rx0(j,1)-rx(i,1)).^2 + (ry0(j,1)-ry(i,1)).^2 +(rz0(j,1)-rz(i,1)).^2).^.5 ;
       
      %-- for self linking number
      tjix = ty0(j,1).*tz(i,1) - tz0(j,1).*ty(i,1);
      tjiy = tz0(j,1).*tx(i,1) - tx0(j,1).*tz(i,1);
      tjiz = tx0(j,1).*ty(i,1) - ty0(j,1).*tx(i,1);
         
      
      temp =   (tjix.*(rx0(j,1)-rx(i,1)) + tjiy.*(ry0(j,1)-ry(i,1)) + tjiz.*(rz0(j,1)-rz(i,1)))./rij.^3 ;
        
      ker(i,1) = trapz(temp);
      
      
      
end

Wr  =   h^2*trapz(ker)/(4*pi);
Lk_edge(p1)  = (tau/(2*pi) + Wr )  ;


        
%         figure(1)
%         plotbrowser on
%         
%         icol =  round((Lk(p1)-1)/2);
%         plot(tau1(p1)/(2*pi),Lk(p1),'-ok','LineWidth',.5,'color',col{icol})
%         
%         
%         
%         hold on
%         set(gca,'FontSize',25,'LineWidth',.5)
%         % hold on
%         
    end
    
    plot(tau1/(2*pi),Lk,'-ok','LineWidth',.5)
    
    title('Linking number')
    
    box on
    grid on
    hold on
    set(gca,'FontSize',25,'LineWidth',.5)
    
%     strLk = ['Linking_branch_' num2str(branch) '.txt'];
%   
    if(sv==1)
        fileID = fopen([path strLk],'w');
        fprintf(fileID,'%30.16E   \r\n',Lk' );
        fclose(fileID);
    end
    
    
%     %=== drawing arrows for direction of the curves and the link number 
    figure(3)
    plotbrowser on
    title('midline')
    hold on
    plot3(rx,ry,rz,'color',col{1},'LineWidth',2)
    hold on 
    
    %== plotting the edges ==
    wd = .01;%
    
    
    plot3(rx+wd*bx,ry+wd*by,rz+wd*bz,'color',col{2},'LineWidth',2)
    hold on 
    plot3(rx-wd*bx,ry-wd*by,rz-wd*bz,'color',col{2},'LineWidth',2)
    hold on 
    
    hold on
    set(gca, 'DataAspectRatio',[1,1,1]);
    daspect([1,1,1]);
    axis off
    set(gcf,'color',colw);
    axis vis3d
    
    nfold = 9 ;
%    
%     for m1 = 1:4
%         
%         ind = 2*m1*N/nfold;
%         diff = 1;
%         ar1 = [rx(ind-diff,1)     , ry(ind-diff,1)     , rz(ind-diff,1)];
%         ar2 = [rx(ind+diff,1)     , ry(ind+diff,1)     , rz(ind+diff,1)];
%      mArrow3(ar1,ar2,'stemWidth',0.001,'tipWidth',0.002);
%     end
%     
    
    
end

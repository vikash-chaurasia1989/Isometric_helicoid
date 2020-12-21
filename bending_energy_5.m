clear all
clc

global    lm  N1 p1 count tau  rho f_b  h tx ty tz  u v w tau1 p1 fac

global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p  N   len rx ry rz Lk err   nfold branch sig w1 w2 w3 tau2

format longE
tic


%---- read r and b file and create ruling data ----
col1 = [0,     0.976, 0.968]  ;      % -- color for m = 3
col2 = [0.831, 0.545, 0.247]  ;      % -- color for m = 5
col3 = [0.341, 0.505, 0.819]  ;      % -- color for m = 7
col4 = [0.705, 0.701, 0.070]  ;      % -- color for m = 9
col5 = [0.301, 0.811, 0.498]  ;      % -- color for m = 11

clr = linspecer(15);

col6  = [0.98, 0.325, 0];%[0.1705, 0.2525,  0.4095]   ;    % -- color for m = 13
col7  = col2/2;% [2.941176470588235e-01     5.447058823529412e-01     7.494117647058823e-01]  ;    % -- color for m = 15 
col8  = [3.717647058823529e-01     7.176470588235294e-01     3.611764705882353e-01]/2 ;    % -- color for m = 17
col9  = [ 1.000000000000000e+00     5.482352941176470e-01     9.999999999999998e-02]  ;    % -- color for m = 19
col10 = col1/1.5;%[8.650000000000000e-01     8.109999999999999e-01     4.330000000000000e-01]   ;    % -- color for m = 21
col11 = [0.533, 0.671, 0.525];%col2/1.5;%[6.858823529411765e-01     4.035294117647059e-01     2.411764705882353e-01]   ;    % -- color for m = 23
col12 = [9.717647058823530e-01     5.552941176470587e-01     7.741176470588236e-01]   ;    % -- color for m = 25
col13 = [0.549, 0.035, 0.486];%[0.196, 0.522, 0.788];%col3/1.5;%[9.717647058823530e-01     5.552941176470587e-01     7.741176470588236e-01]/2 ;    % -- color for m = 27
col14 = [0.996, 0.773, 0.243];%col5/1.5;%[     4.350351799515555e-01     7.800722978266937e-01     6.461779788838612e-01]  ;  % -- color for m = 29                                                                        % -- color for m = 29
col15 = [0.11, 0.18, 1];%[0.988, 0.427, 0.412];%col6/1.5;%[     4.350351799515555e-01     7.800722978266937e-01     6.461779788838612e-01]/3  ;  % -- color for m = 31                                                                      % -- color for m = 29



col = {col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12, col13,col14,col15} ;
%col ={clr(1,:),clr(2,:),clr(3,:),clr(4,:),clr(5,:),clr(6,:),clr(7,:),clr(8,:),clr(9,:),clr(10,:),clr(11,:),clr(12,:),clr(13,:),clr(14,:),clr(15,:)};

 
%col = {col1,col2,col3,col4,col5,clr(1,:),clr(2,:),clr(3,:),clr(4,:),clr(5,:),clr(6,:)} ;
colw = [41,44,52]/255;






%==========================================================================
N = 105;
N1 = N-1;
h = 1/N;

nfold=3;
sig = .01;

t1 = [0.1 0.2 0.3 .4] ;

b_ar = [1,2,3,4,5,6,8,9,11];
b_ar = 8;
% = [1,4,6];
%b_ar = [1,4,8];
figure(1)

for num_ar =  1:length(b_ar)
    branch =  b_ar(num_ar);
    path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];
  %  path = ['/Users/rtodres/Documents/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];

    strtau = ['tau_branch_' num2str(branch) '.txt'];
    strLk = ['Linking_branch_' num2str(branch) '.txt'];
    tau1 = load([path strtau]);
   % if(branch<10)
%    Lk = load([path strLk]);
%     else
%         Lk = 5*length(tau1);
%     end
%     if(branch==13)
%         Lk = 5*ones(1,length(tau1));
%     else
%         Lk = load([path strLk]);
%     end
    
      Lk = load([path strLk]);
   % t = t1(num_ar);
    
    for p1  = 1:length(tau1)%262:length(tau1)%1:length(tau1)%3:length(tau1)%48:length(tau1)
        sv = 1;
        tau =   tau1(p1);
        
        
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
            
      
      %  path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];
     %   path = ['/Users/rtodres/Documents/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];

        str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '.txt'];
        
        temp = load([path str0]);
        
        x = temp(1:3*N1,1);
        
        
        bx = zeros(N+1,1);
        by = zeros(N+1,1);
        bz = zeros(N+1,1);
        
        
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
         
        E(branch,p1) = energy_b(bx,by,bz);
        
    %    icol =  round((Lk(p1)-1)/2);
        
        %       figure(1)
        %       plotbrowser on
        % %     plot(tau1(p1)/(2*pi),E(branch,p1),'.','Markersize',25,'color',col{icol})
        %       plot(tau1(p1)/(2*pi),E(branch,p1),'--' ,'color',col{icol})
        %
        %  plot3(tau1(p1)/(2*pi),t, E(branch,p1),'.','Markersize',25,'color',col{icol})
        
        hold on
        
        
        
    end
    
       ext = 53;  % new points added to branch 4 
    if(branch==1)
        ind = 1:54;
        Lk1 = 3;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        ind = 54:192;
        Lk1 = 7;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        ind = 192:319;
        Lk1 = 11;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        ind = 320;
        Lk1 = 13;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        ind = 319:385;%721;%length(tau1);
        Lk1 = 15;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        ind = 385:721;%length(tau1);
        Lk1 = 15;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        
        ind = 721:723;%length(tau1);
        Lk1 = 17;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        ind = 723:726;%length(tau1);
        Lk1 = 19;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        ind = 726:729;%length(tau1);
        Lk1 = 21;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        ind = 729:length(tau1);
        Lk1 = 23;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
    elseif(branch==2)
        ind = 1:48;
        Lk1 = 5;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        ind = 49;
        Lk1 = 7;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        ind = 49:223;
        Lk1 = 11;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
     elseif(branch==3)
        ind = 1:46;
        Lk1 = 7;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        ind = 46:49;
        Lk1 = 11;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
%         ind = 49:313;
%         Lk1 = 15;
%         icol =  round((Lk1-1)/2);
%         plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
%         hold on 
        ind = 49:494;
        Lk1 = 15;   
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        ind = 495:585;
        Lk1 = 23;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        ind = 586:696;
        Lk1 = 31;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        

    elseif(branch==4) 
        ext = 53;
        ind = 1:118;
        Lk1 = Lk(1) ;  % 5
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        ind = 118:248;
        Lk1 = Lk(ind(end)) ;  % 9 
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        ind = 249;
        Lk1 = Lk(ind(end)) ; %11
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        ind = 249:375;
        Lk1 = Lk(ind(end)) ; % 13 
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        ind = 375:680;
        Lk1 = Lk(ind(end)) ;  % 17
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
%         ind = (376+ext):(617+ext);
%         Lk1 = Lk(ind(end)) ;  % 17
%         icol =  round((Lk1-1)/2);
%         plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
%         hold on 
        
        ind = 680:682;
        Lk1 = Lk(ind(end)) ;   % 19
        icol = round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        ind = 682:length(tau1);
        Lk1 = Lk(ind(end)) ;   % 27
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
    elseif(branch==5)
        ind = 1:58;
        Lk1 = Lk(1) ;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        ind = 58:63;
        Lk1 = Lk(ind(end)) ;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        ind = 63:128;
        Lk1 = Lk(ind(end)) ;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        ind = 128:235;
        Lk1 = Lk(ind(end)) ;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        ind = 235:242;
        Lk1 = Lk(ind(end)) ;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        ind = 242:360;
        Lk1 = Lk(ind(end)) ;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
         
        ind = 360:371;
        Lk1 = Lk(ind(end)) ;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        ind = 371:386;
        Lk1 = Lk(ind(end)) ;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        ind = 386:622;
        Lk1 = Lk(ind(end)) ;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        ind = 622:636;
        Lk1 = Lk(ind(end)) ;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        
    elseif(branch==6) 
        
       E(6,1:126) = E(4,193:318);
        
        ind = 1:52;
        Lk1 = Lk(1) ;  % 9
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 

        ind = 52:60;
        Lk1 = Lk(ind(end)) ; %11
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        ind = 60:274;
        Lk1 = Lk(ind(end)) ; % 13
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        ind = 274:303;
        Lk1 = Lk(ind(end)) ; % 17
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        ind = 303:397;
        Lk1 = Lk(ind(end)) ;  % 21
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        
        ind = 397:length(tau1);
        Lk1 = Lk(ind(end)) ;  % 29
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        
        
    elseif(branch==8)
        ind = 1:61;     
        Lk1 = Lk(ind(end)) ;   % 15;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        ind = 61:62;
        Lk1 = Lk(ind(end)) ;   % 17
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        ind = 62:185;
        Lk1 = Lk(ind(end)) ;   % 19
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        ind = 185:186;
        Lk1 = Lk(ind(end)) ;  % 21
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        
        ind = 186:308;
        Lk1 = Lk(ind(end)) ;  % 23
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        ind = 308:309;
        Lk1 = Lk(ind(end)) ;  % 27
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        ind = 309:length(tau1);
        Lk1 = Lk(ind(end)) ;  % 27
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
    elseif (branch==9)
        
        ind = 1:40;     
        Lk1 = Lk(ind(end)) ;   % 9;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        ind = 40:50;     
        Lk1 = Lk(ind(end)) ;   % 15;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        ind = 50:652;     
        Lk1 = Lk(ind(end)) ;   % 19;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        ind = 652:length(tau1);     
        Lk1 = Lk(ind(end)) ;   % 29;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
    elseif (branch==11)
        
        ind = 1:60;     
        Lk1 = Lk(ind(end)) ;   % 17;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        ind = 60:length(tau1);     
        Lk1 = Lk(ind(end)) ;   % 21;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
    elseif (branch==12)
        
        ind = 1:184;     
        Lk1 = Lk(ind(end)) ;   % 19;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
        
        ind = 184:length(tau1);     
        Lk1 = Lk(ind(end)) ;   % 29;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
    elseif (branch==13)
        
        ind = 1:191;     
        Lk1 =  Lk(ind(end)) ;   % 21;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on     
        
        ind = 191:length(tau1);     
        Lk1 =  Lk(ind(end)) ;   % 29;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on  
     elseif (branch==14)
        
        ind = 1:60;     
        Lk1 =  Lk(ind(end)) ;   % 21;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on     
        
        ind = 60:length(tau1);     
        Lk1 =  Lk(ind(end)) ;   % 25;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on  
       elseif (branch==15)
        
        ind = 1:60;     
        Lk1 =  Lk(ind(end)) ;   % 25;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on     
        
        ind = 60:length(tau1);     
        Lk1 =  Lk(ind(end)) ;   % 29;
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on  
        
    else
        ind = 1:length(tau1);
        Lk1 = 3;%Lk(ind(end)) ;  % 27
        icol =  round((Lk1-1)/2);
        plot(tau1(ind)/(2*pi),E(branch,ind),'-','LineWidth',3,'color',col{icol})
        hold on 
    end
    
    
    
    %plot(tau1/(2*pi),E(branch,1:length(tau1)),'--','LineWidth',.5)
    %plot3(t,tau1/(2*pi),E(branch,1:length(tau1)),'--k','LineWidth',.5)
     
    xlim([1.28,85/2/pi])
    box on
    set(gca,'FontSize',25,'LineWidth',.5)
   % set(gcf,'color', colw );
   % set(gca,'XColor','w');
    %
   % set(gca,'color',colw);
   % set(gca,'XColor','w');
    %set(gca,'YColor','w');%hold on
    %=== labels ===
    xlabel('$n$','FontSize',25,'Interpreter','latex');
    ylabel('$F_b$','FontSize',25,'Interpreter','latex');
   % title('Bending energy for $\sigma=0.01$','FontSize',25,'Interpreter','latex','color','w');
    title('Bending energy for $\sigma=0.01$','FontSize',25,'Interpreter','latex');

    
    grid on
    hold on
    
    %=== saving the torsion values for a given branch
    
end

% %== for legend ==
 Lk1 = [3,5,7,9,11,13,15,17,19,21,23,25,27,29,31];

 
fig2 = figure(2)
for i = 1:length(Lk1)


   plot(2,3  , '.','Markersize',45,'color',col{i},'Visible', 'on')
    hold on
end
 Legend=cell(length(Lk1),1);
 for iter=1:length(Lk1)
   Legend{iter}= strcat('Lk=', num2str(Lk1(iter)) );
 end
 leg = legend(Legend,'FontSize',30,'Interpreter','latex' )
 ax = findobj(gcf,'type','axes'); %Retrieve the axes to be copied

fig1 = figure(1);
copyobj([leg,ax],fig1);



%====  saving energy data 
   %
% 
%    str = 'Energy_branches.txt';
%    branch = 1;
%    path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];
%  
%    fileID = fopen([path str],'w');
%    
%    strform = [];
%    
%    for i = 1:max(b_ar)
%        strform = ['%30.16E ' strform];
%    end
%    
%    strform = [strform ' \r\n'];
%    
%    fprintf(fileID,strform,E   );
%    fclose(fileID);
%     
    
    
%======================= Boltzmann distribution ===========================
%

% E = E/55;
% 
% P1 =   exp(-E(1,76:261))./(exp(-E(1,76:261)) +  exp(-E(2,38:223)) + exp(-E(3,3:188)));
% P2 =   exp(-E(2,38:223))./(exp(-E(1,76:261)) +  exp(-E(2,38:223)) + exp(-E(3,3:188)));
% P3 =   exp(-E(3,3:188))./(exp(-E(1,76:261)) +  exp(-E(2,38:223)) + exp(-E(3,3:188)));
% 
%==========================================================================
% figure(2)
% %
% plot(tau1(3:188),P1);
% hold on
% plot(tau1(3:188),P2);
% hold on
% plot(tau1(3:188),P3);
% 
% set(gca,'FontSize',25,'LineWidth',.5)
% hold on
% title('Probability')
% 
% box on
% grid on
% hold on

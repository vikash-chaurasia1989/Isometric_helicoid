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

clr = linspecer(7);

col6  = [9.047058823529412e-01     1.917647058823529e-01     1.988235294117647e-01]  ;    % -- color for m = 13
col7  = col2/2;% [2.941176470588235e-01     5.447058823529412e-01     7.494117647058823e-01]  ;    % -- color for m = 15 
col8  = [3.717647058823529e-01     7.176470588235294e-01     3.611764705882353e-01]/2  ;    % -- color for m = 17
col9  = [ 1.000000000000000e+00     5.482352941176470e-01     9.999999999999998e-02] ;    % -- color for m = 19
col10 = [8.650000000000000e-01     8.109999999999999e-01     4.330000000000000e-01]  ;    % -- color for m = 21
col11 = [6.858823529411765e-01     4.035294117647059e-01     2.411764705882353e-01]  ;    % -- color for m = 23
col12 = [9.717647058823530e-01     5.552941176470587e-01     7.741176470588236e-01]  ;    % -- color for m = 27


col = {col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12, col12/2} ;

%col = {col1,col2,col3,col4,col5,clr(1,:),clr(2,:),clr(3,:),clr(4,:),clr(5,:),clr(6,:)} ;
colw = [41,44,52]/255;

clg = [134,136,138]/255;




%==========================================================================
 
sig = .01;

t1 = [0.1 0.2 0.3 .4] ;

b_ar = [1,2,3,4,5,6];
 b_ar  = [1,4,6];

figure(1)

for num_ar = 1:2
    branch = b_ar(num_ar);
    path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];
    path = ['/Users/rtodres/Documents/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];

    strtau = ['tau2_branch_' num2str(branch) '.txt'];
    strLk = ['Linking_branch_' num2str(branch) '.txt'];
    strst = ['stability_branch_' num2str(branch) '.txt'];
    tau1 = load([path strtau]);
    Lk   = load([path strLk]);
    st   = load([path strst]);
   % t = t1(num_ar);
    
    for p1  = 1:length(tau1)%262:length(tau1)%1:length(tau1)%3:length(tau1)%48:length(tau1)
        sv = 1;
        tau =   tau1(p1);
        
        
        if (branch==1||branch==2||branch==3)
            N=105;
        end
        
        if(branch==4)
            N = 120;
        end
        if(branch==5)
            N = 140;
        end
        
        if(branch== 6)
            N = 90;
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
        
  
        
        if(branch==6&&p1<127)
            E(6,p1) = E(4,128+p1);
        end
        
        if(st(p1)==1)
            plot(tau1(p1)/(2*pi), E(branch,p1),'.','Markersize',10,'color',clg)%col{branch})
            hold on 
        else 
              plot(tau1(p1)/(2*pi), E(branch,p1),'.','Markersize',10,'color','r')
              hold on 
              
        end
        
    end
    
    
    
        
end

 
    set(gca,'FontSize',25,'LineWidth',.5)
    hold on
    title('Bending energy')
    
    box on
    grid on
    hold on
    
    
    
    %== plotting the global minimum energy ====
    
    
tau_int = [1.510000000000000e+01,1.750000000000000e+01, 2.410000000000000e+01,3.060000000000000e+01,3.700000000000000e+01,4.320000000000000e+01 ];



branch = 1; 
    path = ['/Users/rtodres/Documents/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];

    strtau = ['tau2_branch_' num2str(branch) '.txt'];
    strLk = ['Linking_branch_' num2str(branch) '.txt'];
    strst = ['stability_branch_' num2str(branch) '.txt'];
    tau1 = load([path strtau]);
        

branch = 4; 
    path = ['/Users/rtodres/Documents/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];

    strtau = ['tau2_branch_' num2str(branch) '.txt'];
    strLk = ['Linking_branch_' num2str(branch) '.txt'];
    strst = ['stability_branch_' num2str(branch) '.txt'];
    tau4 = load([path strtau]);   
    
    
    
    
branch = 1; 

[ind,temp] = find(tau1<=tau_int(1));

plot(tau1(ind)/2/pi, E(branch,ind),'o--','Markersize',5,'color','k')
hold on 
x1 = tau1(ind(end))/2/pi, y1= E(branch,ind(end));


 

branch = 4; 

[ind,temp] = find(tau_int(1)<=tau4 & tau4<=tau_int(2));

plot(tau4(ind)/2/pi, E(branch,ind),'o--','Markersize',5,'color','k')
hold on 

 

 

x2 = tau4(ind(1))/2/pi, y2= E(branch,ind(1));

plot([x1,x2],[y1,y2],'--', 'Linewidth',5)

branch = 1; 

[ind,temp] = find(tau_int(2)<=tau1 & tau1<=tau_int(3));

plot(tau1(ind)/2/pi, E(branch,ind),'o--','Markersize',5,'color','k')
hold on 


branch = 4; 

[ind,temp] = find(tau_int(3)<=tau4 & tau4 <=tau_int(4));

plot(tau4(ind)/2/pi, E(branch,ind),'o--','Markersize',5,'color','k')
hold on 

branch = 1; 

[ind,temp] = find(tau_int(4)<=tau1 & tau1<=tau_int(5));

plot(tau1(ind)/2/pi, E(branch,ind),'o--','Markersize',5,'color','k')
hold on 



branch = 4; 

[ind,temp] = find(tau_int(5)<=tau4 & tau4<=tau_int(6));

plot(tau4(ind)/2/pi, E(branch,ind),'o--','Markersize',5,'color','k')
hold on 


branch = 1; 

[ind,temp] = find(tau_int(6)<=tau1 & tau1<=tau1(end));

plot(tau1(ind)/2/pi, E(branch,ind),'o--','Markersize',5,'color','k')
hold on 


x1 = tau1(ind(end))/2/pi, y1= E(branch,ind(end));


branch = 4; 

[ind,temp] = find(tau1(end)<=tau4 & tau4<=tau4(end));

plot(tau4(ind)/2/pi, E(branch,ind),'o--','Markersize',5,'color','k')
hold on 


x2 = tau4(ind(1))/2/pi, y2= E(branch,ind(1));

plot([x1,x2],[y1,y2],'--', 'Linewidth',5)
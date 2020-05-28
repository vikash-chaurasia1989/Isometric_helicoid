clear all
clc

global    lm  N1 p1 count tau  rho f_b  h tx ty tz  u v w tau1 p1 fac

global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p  N   len rx ry rz Lk err   nfold branch sig w1 w2 w3 tau2

format longE
tic


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



 

 
%==========================================================================
N = 105;
N1 = N-1;
h = 1/N;

nfold=3;
sig = .01;

t1 = [0.1 0.2 0.3];

b_ar = [2,1,3];

for num_ar = 1:3
    branch = b_ar(num_ar);
    path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];
    strtau = ['tau_branch_' num2str(branch) '.txt'];
    strLk = ['Linking_branch_' num2str(branch) '.txt'];
    tau1 = load([path strtau]);
    
    Lk = load([path strLk]);
    
    t = t1(num_ar);

for p1  =  1:length(tau1)%262:length(tau1)%1:length(tau1)%3:length(tau1)%48:length(tau1)
    sv = 1;
    tau =   tau1(p1);
    
    
    
    path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];
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
    
    
    
    
    E(branch,p1) = energy_b(bx,by,bz);
    
    icol =  round((Lk(p1)-1)/2);
    
    figure(1)
    plotbrowser on
   % plot(tau1(p1)/(2*pi),E(branch,p1),'.','Markersize',25,'color',col{icol})
      plot3(tau1(p1)/(2*pi),t, E(branch,p1),'.','Markersize',25,'color',col{icol})

    hold on 
    
        
   
end
%plot(tau1/(2*pi),E(branch,1:length(tau1)),'--k','LineWidth',.5)
%plot3(t,tau1/(2*pi),E(branch,1:length(tau1)),'--k','LineWidth',.5)

set(gca,'FontSize',25,'LineWidth',.5)
hold on
title('Bending energy')

box on
grid on
hold on

%=== saving the torsion values for a given branch

end
 
% %== for legend ==
% fig2 = figure(2)
% for i = 1:7
%     
%    
%    plot(1, 1, '.','Markersize',45,'color',col{i},'Visible', 'on')
%     hold on 
% end
%  Legend=cell(7,1);
%  Lk1 = [3,5,7,9,11,13,15];
%  for iter=1:7
%    Legend{iter}= strcat('Lk=', num2str(Lk1(iter)) );
%  end
%  leg = legend(Legend,'FontSize',40,'Interpreter','latex' )
%  ax = findobj(gcf,'type','axes'); %Retrieve the axes to be copied
%  
% fig1 = figure(1);
% copyobj([leg,ax],fig1);
%  
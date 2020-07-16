clear all
clc

global    lm  N1 p1 count tau  rho f_b  h tx ty tz  u v w tau1 p1 fac

global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p  N   len rx ry rz Lk err   nfold branch sig w1 w2 w3 tau2

format longE
tic


%---- read r and b file and create ruling data ----
  col1 = [105, 105, 105]/255  ;      % -- color for mode n = 2
  col2 = [112,128,144]/255  ;      % -- color for mode n = 3
  col3 = [47,79,79]/255  ;      % -- color for mode n = 4
   
  col = {col1,col2,col3} ;
  col = {[0,0,0],[0,0,1],[1,0,0]};
colw = [41,44,52]/255;



 

 
%==========================================================================
N = 105;
N1 = N-1;
h = 1/N;

nfold=3;
sig = .01;

figure(1)
sig1 = [0.01 0.02 0.03];
for branch = 1:3
    path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];
    strtau = ['tau_branch_' num2str(branch) '.txt'];
    strLk = ['Linking_branch_' num2str(branch) '.txt'];
    tau1 = load([path strtau]);
    
    Lk = load([path strLk]);

for num_sig = 1:length(sig1)
    sig = sig1(num_sig);
    
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
    
    
    
    
    E(branch,p1,num_sig) = energy_b(bx,by,bz);
    
   
    
    
    
%     plotbrowser on
%     plot(tau1(p1)/(2*pi),E(branch,p1,num_sig),'.','Markersize',25,'color',col{num_sig})
%     hold on 
%     
        
   
end
subplot(1,3,branch)
%plot(tau1/(2*pi),E(branch,1:length(tau1),num_sig),'.','Markersize',25,'color',col{num_sig})
plot(tau1/(2*pi),E(branch,1:length(tau1),num_sig),'color',col{num_sig})

set(gca,'FontSize',25,'LineWidth',.5)
hold on
pbaspect(subplot(1,3,branch),[1.61 1 1])
hold on
hold on
%title('Bending energy')

box on
grid on
hold on

%=== saving the torsion values for a given branch

end

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
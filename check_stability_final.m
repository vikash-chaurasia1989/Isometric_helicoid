%==== This file checks stability
% In this file, we check stability of the data in data_branch1, 2, 3. The
% solutions are obtained for N=105 and using fun_jacobian6 where only end
% points are fixed and the binormals are allowed to rotate

clear all
clc

global    lm  N1 p1 count tau  rho f_b  h tx ty tz  u v w tau1 p1

global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p  N   len rx ry rz Lk err   nfold branch sig w1 w2 w3 tolB

format longE
tic
 
N = 105;
N1=N-1;
h = 1/N;

for branch = 3:3
    path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];
    strtau = ['tau_branch_' num2str(branch) '.txt'];
    tau1 = load([path strtau]);
    %==========================================================================
    p2 = 147;
    st = zeros(1,length(tau1));
    %  st = zeros(1,p2);
    
    
    %==========================================================================
    for p1  = 1:length(tau1)
        sv = 0;
        tau =   tau1(p1);
        
        
        path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];
        str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '.txt'];
        
        x= load([path str0]);
        
       
        
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
        
        
        %======= calling fun_stability2 ====
        
        tolB = 10^-6;%   % tolerance for deciding if the eigen value is zero or negative 
                         % if min((eig))<-tolB then the solution is
                         % unstable
                         % 
        st(p1) = fun_stability6(x);
        
        
        
    end
    plot(st,'-o','LineWidth',.5)
    hold on
    
    strst = ['stability_branch_' num2str(branch) '.txt'];
    
    if(sv==1)
        fileID = fopen([path strst],'w');
        fprintf(fileID,'%30.16E   \r\n',st );
        fclose(fileID);
    end
end
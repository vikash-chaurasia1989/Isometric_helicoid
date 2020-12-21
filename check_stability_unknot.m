clear all
clc

global N  N1 h tau bx by bz id B_temp p1 B


id = eye(3,3);





branch =4;
path   = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];
%path = ['/Users/rtodres/Documents/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];

strtau = ['tau2_branch_' num2str(branch) '.txt'];
tau1 = load([path strtau]);

sv=1;


for p1 =  1:64%length(tau1)
    
    if (branch==1||branch==2||branch==3)
        N=105;
    end
    
    if(branch==4||branch==9)
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
    
    
    
    bx = zeros(N+1,1);
    by = zeros(N+1,1);
    bz = zeros(N+1,1);
    
    % tau1 = 8.8:.1:48;
    
    bx(1,1) = 1;
    by(1,1) = 0;
    bz(1,1) = 0;
    
    bx(N+1,1) = -1;
    by(N+1,1) =  0;
    bz(N+1,1) =  0;
    
    
    
    
    
    tau = tau1(p1);
    
    str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '.txt'];
    
    x = load([path str0]);
    
    st(p1) = fun_stability6(x);
    
    eigB(p1) = B;
end


if(sv==1)
    str = ['stability_branch_' num2str(branch) '.txt'];
    
    fileID = fopen([path str],'w');
    fprintf(fileID,'%30.16E   \r\n',st' );
    fclose(fileID);
    
end
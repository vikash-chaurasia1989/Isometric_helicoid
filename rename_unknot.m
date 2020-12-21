%== This function renames the saved data unknot

clear all
clc

N = 120;
branch =4;

tau1 = 8.8:.1:48;


for p1 = 1:length(tau1)
    tau = tau1(p1);
    
    
    
    path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_unknot_5pi/' ;
    
    str0 = ['branch_unknot_5pi_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '.txt'];
    
    
    x = load([path str0]);
    
    %=== renaming ==
    path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];
    
    str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '.txt'];
    
    
    
    fileID = fopen([path str0],'w');
    fprintf(fileID,'%30.16E   \r\n',x );
    fclose(fileID);
end
    strtau = ['tau_branch_' num2str(branch) '.txt'];
    
    fileID = fopen([path strtau],'w');
    fprintf(fileID,'%30.16E   \r\n',tau1' );
    fclose(fileID);


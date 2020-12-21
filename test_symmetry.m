clear all
clc





N = 120;
N1 = N-1;
branch = 4;

tau =   15.1;
path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];
% path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_unknot_5pi/' ;

   path = ['/Users/rtodres/Documents/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];

str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '.txt'];

x = load([path str0]);
bx(2:N,1)     = x(1:3:3*N1-2,1)    ;
by(2:N,1)     = x(2:3:3*N1-1,1)    ;
bz(2:N,1)     = x(3:3:3*N1,1)      ;
bx(1,1) = 1;
by(1,1) = 0;
bz(1,1) = 0;

bx(N+1,1) = -1;
by(N+1,1) =  0;
bz(N+1,1) =  0;


for i = 1:N+1
    [t1,t2] = cart2sph2(bx(i,1),by(i,1),bz(i,1));
    th(i) = t1;
    ph(i) = t2;
end

% ph(62:end) = ph(62:end) + ph(61)-ph(62);
% plot(ph)

%=========
% 
% 
% ph2 = ph;
% 
% ph2(142:(142+70)) = ph(1:71);

th2 = th;
ph2 = ph;



ind = N+1-N/5;

th((N+1)+ (1:ind)) = th(ind:-1:1);

ph(62:85) = 2*ph(62)-ph(62:85);
ph(86:end) = 2*pi+ ph(86:end);




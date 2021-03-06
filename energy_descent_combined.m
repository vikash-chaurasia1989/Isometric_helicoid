clear all
clc

global N N1 h tau sig fac branch   path  nstep strpath 


smop = {'moving','lowess','lowess','sgolay','rlowess','rloess'};
opid  = 3;


col1 = [0,     0.976, 0.968]  ;      % -- color for mode n = 2
col2 = [0.831, 0.545, 0.247]  ;      % -- color for mode n = 3
col3 = [0.341, 0.505, 0.819]  ;      % -- color for mode n = 4
col4 = [0.705, 0.701, 0.070]  ;      % -- color for mode n = 5
col5 = [0.301, 0.811, 0.498]  ;      % -- color for mode n = 6
col6 = [1, 0.929, 0.278];
col = {col1,col2,col3,col4,col5,col6} ;
colw = [41,44,52]/255;



parameters5(); 


%== prefactor for aspect ratio
 

 branch = 2;
 
 num1 = 3255;
 num3 = 2244;
         
for p1= 1:(num1+num3) %1:253%length(t1)
    
    str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*16)) '_step_' num2str(p1-1)   '.txt'];

    if(p1>num1)
        branch =3;
        str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*16)) '_step_' num2str(num1+num3-p1)   '.txt'];
    end
    
    

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
% % 
%     bx = smooth(bx,smop{opid});
%     by = smooth(by,smop{opid});
%     bz = smooth(bz,smop{opid});
% %     
    bx(1,1) = 1;
    by(1,1) = 0;
    bz(1,1) = 0;
    bx(N+1,1) = -1;
    by(N+1,1) =  0;
    bz(N+1,1) =  0;
    
    E(p1) = energy_b(bx,by,bz);
    
    lenE(p1,:) = sqrt(bx.^2+by.^2 + bz.^2)-1;
end

branch = 3;
 path1 = [strpath num2str(branch) '/'];

strtau = ['tau_branch_' num2str(branch) '.txt'];

tau1 = load([path1 strtau]);

[val,ind] = min(abs(tau1-tau));

tau2 = tau1(ind);


str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau2)) '.txt'];


x = load([path1 str0]);

bx(2:N,1)     = x(1:3:3*N1-2,1)    ;
by(2:N,1)     = x(2:3:3*N1-1,1)    ;
bz(2:N,1)     = x(3:3:3*N1,1)      ;

bx(1,1) = 1;
by(1,1) = 0;
bz(1,1) = 0;

bx(N+1,1) = -1;
by(N+1,1) =  0;
bz(N+1,1) =  0;
% 
% %E(p1+1) = energy_b(bx,by,bz);

%
%==== E1 = 1.091317488162173e+00;  ====
%==== E3 = 4.347956445701464e+00
%==== E2 = 2.853737001466909e+00   ====
%


plot(E,'-o')
hold on
plot(E(1),'-or')
hold on
%plot(p1+1,E(end),'-or')

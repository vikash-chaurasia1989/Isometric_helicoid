   function f =   initial_guess()
 
  
  global   N1  N  qd   id   h   rho lm  tau
 
  global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p
  
  %--------Array initialization -------
 
   %  
   N = 72;
   tau2=12.55;
   
 str0 = ['3fold_N' num2str(N) '_tau_' num2str(10*tau2) '_1.txt'];
strb = ['b_' str0];
strr = ['r_' str0];

strlm = ['lm_' str0];
strrho = ['rho_' str0];
struvw = ['uvw_' str0];
strbext = ['bext_' str0];
% % % % %    


       temp = load(strb);
 
           bx = temp(:,1);
           by = temp(:,2);
           bz = temp(:,3);
%
%               N = 72;
%      
%      s = linspace(0,1,N+1);
%      
%      th(:,1) = pi/2- .8*sin(3*pi*s);
%      ph(:,1) =  pi*s;
%      
%      bx(:,1) = sin(th).*cos(ph);
%      by(:,1) = sin(th).*sin(ph);
%      bz(:,1) = cos(th);
% % % %
%     N = length(bx)-1;
%     N1 = N-1;
%     

    
%        N = length(b)-1;
%        
%      bx = b(1:N/2+1,1);
%        by = b(1:N/2+1,2);
%        bz = b(1:N/2+1,3);
%        N = N/2;
     
%      r = (bx.^2+by.^2 + bz.^2).^.5;
%      bx = bx./r;
%      by = by./r;
%      bz = bz./r;
%      
   N = length(bx)-1;
     
     
%   

   
       u = -5.1566319025341907e-01;% 1.054341983570216e-11 ;
        
       v =  -8.00693449731905190e+00;%-5.308099190019070e-12 ;
 
       w =   3.9967030985498e-01;%-0.01762723061619;
       
       
 %  v = w;       
        temp = load(struvw);
       
       u  = temp(1,1);
       v  = temp(2,1);
       w = 0*temp(3,1);

      
    
%   N = length(bx)-1;
 %   N1 = N-1;
% %  
   
   
%   %    %--- refinement --     
% % %  
%       bx1(1:2:2*(N+1),1)   = bx(1:N+1,1);
%       bx1(2:2:2*(N+1)-1,1) =(bx(1:end-1,1)+bx(2:end,1))/2;
%       
%       by1(1:2:2*(N+1),1)   = by(1:N+1,1);
%       by1(2:2:2*(N+1)-1) =(by(1:end-1,1)+by(2:end,1))/2;
%       
%       bz1(1:2:2*(N+1),1)   = bz(1:N+1,1);
%       bz1(2:2:2*(N+1)-1,1) =(bz(1:end-1,1)+bz(2:end,1))/2;
%       
%       bx = bx1;
%       by = by1;
%       bz = bz1;
%       
      
 N = length(bx)-1;

 N1 = N-1;
%  
%       
           

%---- coarsing the data ---
% 
%        N2 = N/2;
%        
%        bx1 = bx(1:2:N+1);
%        by1 = by(1:2:N+1);
%        bz1 = bz(1:2:N+1);
%        
%        bx = bx1;
%        by = by1;
%        bz = bz1;
%        
%        N = length(bx)-1;
%        N1 = N-1;
%        
       %--- plotting initial guess 
       figure(1)
       
       plotbrowser on 
       title('Initial guess')
       plot3(bx,by,bz,'LineWidth',4)
       set(gca,'DefaultAxesFontSize',4,'FontSize',25,'LineWidth',2)
       box on
       grid on
  
  
     
 
   b(1:3:3*N1-2) = bx(2:N,1)  ;
   b(2:3:3*N1-1) = by(2:N,1)  ;
   b(3:3:3*N1)    = bz(2:N,1)  ;
   
 
   h = 1/N ;
%   h = 1/N  ;
%    
%    
%    u = 0;
%    v = 0;
%    w = .1;
   
   rho(1:N-1,1) = .01*ones(N1,1);
%
  rho = load(strrho);
% 
%  rho = [(rho(1:end-1,1)+rho(2:end,1))'/2 rho(1,1)]';
%   rho = (rho(1:end-1,1)+rho(2:end,1))/2;
%             rho(N-1,1) = rho(1,1);
%             rho(N-2,1) = rho(2,1);
%   

    %  lm(1:N+1,1) = -1/(h*tau*(N+1))*ones(N+1,1);
%
       lm = load(strlm);
   
   N  = length(lm)-1;
%   
 %  lm = (lm(1:end-1,1)+lm(2:end,1))/2;
  %          lm(N+1,1) = lm(1,1);
%              
%             
%---- refinement of lagrange multipliers ---
%       lm1(1:2:2*(N+1),1)     = lm(1:N+1,1);
%       lm1(2:2:2*(N+1)-1,1) = (lm(1:end-1,1) + lm(2:end,1))/2;
%       
%       lm = lm1;
%       
%       rho1(1:2:2*(N-1),1)     = rho(1:N-1,1);
%       rho1(2:2:2*(N-1)-1,1) = (rho(1:end-1,1) + rho(2:end,1))/2;
%       
%       rho = rho1;
%      
%      N = length(lm)-1;
%      N1 = N-1;

%--- linear interpolation
bx1 = 2*bx(1,1)-bx(2,1);
by1 = 2*by(1,1)-by(2,1);
bz1 = 2*bz(1,1)-bz(2,1);

bxN2 = 2*bx(N+1,1)-bx(N,1);
byN2 = 2*by(N+1,1)-by(N,1);
bzN2 = 2*bz(N+1,1)-bz(N,1);


%--- higher order interpolation ---
nfit = 4;

%-- At i = 1;
px = polyfit(linspace(1,nfit+1,nfit+1),bx(1:nfit+1,1)',nfit);
bx1 = px(end);

px = polyfit(linspace(1,nfit+1,nfit+1),by(1:nfit+1,1)',nfit);;
by1 = px(end);

px = polyfit(linspace(1,nfit+1,nfit+1),bz(1:nfit+1,1)',nfit);;
bz1 = px(end);

%-- At i = N+1
px = polyfit(linspace(1,nfit+1,nfit+1),bx(N+1:-1:N+1-nfit,1)',nfit);
bxN2 = px(end);

px = polyfit(linspace(1,nfit+1,nfit+1),by(N+1:-1:N+1-nfit,1)',nfit);
byN2 = px(end);

px = polyfit(linspace(1,nfit+1,nfit+1),bz(N+1:-1:N+1-nfit,1)',nfit);
bzN2 = px(end);


%-- extrapolated points from saved data file 
bext = load(strbext);

bx1 = bext(1);
by1 = bext(2);
bz1 = bext(3);

bxN2 = bext(4);
byN2 = bext(5);
bzN2 = bext(6);

%--- for fun_curve7 ---
  f = [b(1:3*N1)  rho' lm' u v w  bx1 by1 bz1 bxN2 byN2 bzN2 ]'  ;
  
  
  %--- fix bz(2) = 0; bz(N) = 0 ---
  
  f(3,:) = [];
  f(3*N1-1,:) = [];
  
  
  %--- boundary points --
  
  bx(1,1) = 1;
  by(1,1) = 0;
  bz(1,1) = 0;
  
  bx(N+1,1) = -1;
  by(N+1,1) =   0;
  bz(N+1,1) =   0;
  
  bz(2,1) = 0;
  bz(N,1) = 0;
 
%--- for fun_curve9 ---
% f = [b(1:3*N1)  rho' lm' u v w  bx1 by1 bz1]'  ;
 
%------------ Array initialization ----------------

  bxp = zeros(N+1,1);
  byp = zeros(N+1,1);
  bzp = zeros(N+1,1);
  
  bx2p = zeros(N+1,1);
  by2p = zeros(N+1,1);
  bz2p = zeros(N+1,1);
 
  bx4p = zeros(N-1,1);
  by4p = zeros(N-1,1);
  bz4p = zeros(N-1,1);
 
   

%   %-- hessian
    qd  = zeros(3*N1+3,3*N1+3);
    id = [1 0 0;0 1 0;0 0 1];
 
  end
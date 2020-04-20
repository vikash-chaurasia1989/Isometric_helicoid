%=== renaming the saved data file 


%===== tau for 3 fold ====
tau1 = [8.1:.1:9.5 10:1:11.9 11.97871 11.98071 11.98571 11.99071 11.99571  12:.1:12.5 12.55 12.552 12.554 12.555 12.556 12.56 12.57 12.6:.1:13.9 14:.5:21];

tau_3 = [8.1:.1:9.5 10:1:11.9 11.97871 11.98071 11.98571 11.99071 11.99571  12:.1:12.5 12.55 12.552 12.554 12.555 12.556 12.56 12.57 12.6:.1:13.2 13.21 13.215 ...
             13.2165 13.21675 13.216865 13.21686525 13.2168653 13.2168654 13.21686545 13.21686546 13.21686547 13.2168654725 13.216865473  13.216865474355 ] ;
tau_5 = [13.216865474 13.2168655 13.21686575 13.216866 13.21687 13.2169 13.217 13.2175 13.22 13.23];
tau_7 = [13.3:.1:13.9 14.5:.5:21];

%tau1 = [8.1:.1:9.5 10:1:11.9 11.97871 11.98071 11.98571 11.99071 11.99571  12:.1:12.5 12.55 12.552 12.554 12.555 12.556 12.56 12.57 12.6:.1:13.2 13.21 13.1:.1:13.9 14:.5:21];
tau_1 = [tau_3 tau_5 tau_7];
N = 72;
h = 1/N;

N1 = N-1;



for i = 1:length(tau1)
    
         str0 = ['3fold_N' num2str(N) '_tau_' num2str(10*tau) '.txt'];

     if(tau>11.9)
          str0 = ['3fold_N' num2str(N) '_tau_' num2str(10*tau) '_1.txt'];
     end
     
     %=== loading the data ===
     
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
           
           
end
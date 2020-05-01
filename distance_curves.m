function d = distance_curves(p,q,tau)

global N 

temp = [p q];

   path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_evolution/';

for i = 1:length(temp)
    
    branch = temp(i);
    
    
  
  if(branch==1)
      str0 = ['3fold_N' num2str(N) '_tau_' num2str(10^10*tau) '.txt'];
  elseif(branch==2)
      str0 =['5pi_knot_N' num2str(N) '_tau_' num2str(10^10*tau) '.txt'];
      %
  else
      str0 = ['7pi_knot_N' num2str(N) '_tau_' num2str(10^10*tau) '.txt'];
  end
  
  
  strb =  [path 'b_' str0] ;
   
  temp2 = load(strb);
 
  bx(:,i)  = temp2(1:N+1,1);
  by(:,i)  = temp2(1:N+1,2);
  bz(:,i)  = temp2(1:N+1,3);
    
end


d = sum(sqrt((bx(:,1)-bx(:,2)).^2 + (by(:,1)-by(:,2)).^2 + (bz(:,1)-bz(:,2)).^2));
end
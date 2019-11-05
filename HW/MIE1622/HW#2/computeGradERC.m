function gval = computeGradERC (x)

global Q
  
  n = size(Q,1) ;  

  if(size(x,1)==1)
     x = x';
  end
  
  % Insert your gradiant computations here
  % You can use finite differences to check the gradient

  % Define the size of the gradient 
  gval = zeros(n,1); 
  testsize = (n-1)*n/2;
  xdiff = zeros(n,testsize); 
  ydiff = zeros(testsize,1); 
  
  % The objective function is defined as following 
  y = x.*(Q*x);

  for i = 1:n
      for j = 1:n
          diff1 = Q(i,:) * x + Q(i,i) * x(i);
          diff2 = Q(i,j) * x(i);
          g = (y(i)-y(j)) * (diff1 - diff2);
          gval(i,:) = gval(i,:) + g;
      end
      gval(i,:) = 4 *  gval(i,:);
  end
end

%   icount = 1; 
%   % Define two matrix, ydiff and xdiff here for Grandient calc
%   for (i = 1:n-1)
%       for (j = i+1:n)
%           ydiff(icount) = y(i)-y(j);
%           icount = icount + 1;
%       end
%   end
%       
%   icount = 1; 
%   for (round1 = 1:n)
%       for (i = 1:n-1)
%           for (j = i+1:n)
%               if (i==round1)
%                   xdiff(round1,icount) = sum(Q(round1,:)*x)+ Q(i,i) * x(i)- Q(i,j)*x(i);
%               else 
%                   xdiff(round1,icount) = Q(i,round1)*x(round1) - Q(j,round1)*x(round1);
%               end 
%               icount = icount + 1;
%           end
%       end
%       icount = 1; 
%   end 
%   
%   gval = 4*xdiff*ydiff;
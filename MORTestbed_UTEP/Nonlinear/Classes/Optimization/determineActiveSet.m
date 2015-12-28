function  [I] = determineActiveSet(x,A,b)

eps = 1e-10;
I = ( abs(A*x - b) < eps );

end


% function  [I] = determineActiveSet(x,l,u)
% 
% eps = 1e-10;
% I = ( abs(x - u) < eps | abs(x - l) < eps );
% 
% end

A = [eye(2), [1;1]];
b = [1;1];

% parameters
% rho = input('value for rho? \n');
% lambda = input('value for lambda? \n');
rho = 1e1;
lambda = 1e-3; 

rng(0)
x = rand(3,1);
z = rand(3,1);
%% optimize
numIter = 1e5;
aleph = 1e-1;
for k  = 1:numIter
   % compute gradient w/r to x  
   gradx = 2*A'*(A*x-b) + rho*(x-z);
   % search along gradient direction 
   x = x-aleph*gradx;
   % minimize along z gradient
   z = shrinkage(rho, lambda, x);
end
disp(x)

%%
function x = shrinkage(rho, lambda, v)

x = abs(v) .* max(0, 1 - lambda ./ (rho*abs(v)));

end
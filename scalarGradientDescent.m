% clc
% fixed parameters
y1 = 1;
y2 = 2;
n1 = 1;     % air
n2 = 1.5;   % glass
L = 4;

% algorithmic properties 
x = 2;          % initial estimate for x
stepSize = 2;   % step size
numIter = 200;   % number of iterations 
for k=1:numIter
    % compute gradient of cost function
    gradCost = n1*x/sqrt(x^2+y1^2) + ...
               n2*(x-L)/sqrt((x-L)^2 + y2^2);
    % walk downhill
    x = x - stepSize * gradCost;
end
% display result with n digits precision
format long
n = 10;
disp('final estimate for x: ')
disp(round(x*10^n)/10^n)
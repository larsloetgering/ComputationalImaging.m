clc
clear
set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
% fixed parameters
y1 = 1;
y2 = 2;
n1 = 1;
n2 = 1.5;
L = 4;
% algorithmic properties 
x = linspace(0,L,50);          % initial estimate for x
costs = zeros(size(x));    % preallocate costs
stepSize = 2;   % step size
numIter = 20;   % number of iterations 
for k=1:length(x)
    % compute cost function
    costs(k) = n1*sqrt(x(k)^2+y1^2) + ...
               n2*sqrt((x(k)-L)^2 + y2^2);
end
% display result with n digits precision
plot(x, costs, 'ko','MarkerFaceColor','k')
xlabel('x'),ylabel('cost'), grid on, axis square
set(0, 'DefaultAxesFontSize', 24);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
set(0, 'DefaultLineLineWidth', 2);
clear 

N = 20;
rng(1)
x = linspace(0,1,N)' + 1e-2*randn(N,1); % data grid
xd = linspace(0,1,5*N)';    % dense grid
y = exp(x);                 % ground truth
data=y + 1/2*randn(N,1);    % noisy observation

figure(1)
plot(x,data,'ko','MarkerFaceColor','k')
axis square

p = 10;
X = zeros(size(x,1),p+1);
Xd = zeros(size(xd,1),p+1);
for k = 0:p
    X(:,k+1) = x.^k;
    Xd(:,k+1) = xd.^k;
end
% least squares
c = inv(X'*X)*(X'*data);
% ridge regression
lambda = 1;
d = inv(X'*X + lambda*eye(p+1,p+1))*(X'*data);

figure(1), clf
hold on 
plot(xd,Xd*c,'-bx')
plot(xd,Xd*d,'-rx')
plot(x,data,'ko','MarkerFaceColor','k')
plot(x,y,'k--')
hold off
axis square
h = legend('least squares','ridge','data','ground truth');
h.Location = 'SoutheastOutside';
axis([0 1 -inf inf])
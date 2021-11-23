set(0, 'DefaultAxesFontSize', 24);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
set(0, 'DefaultLineLineWidth', 2);
clear 

N = 50;                     % number of data points
rng(1)
x = linspace(0,1,N)';       % data grid
xd = linspace(0,1,5*N)';    % dense grid
y = exp(x);                 % ground truth
data = y + 1/2*randn(N,1);  % noisy observation

figure(1)
plot(x,data,'ko','MarkerFaceColor','k')
axis square
%%
p = 10;
Gimmel = toeplitz(-eye(p+1-1,1), [-1,1,zeros(1,p+1-2)]);
%%
X = zeros(size(x,1),p+1);
Xd = zeros(size(xd,1),p+1);
for k = 0:p
    X(:,k+1) = x.^k;
    Xd(:,k+1) = xd.^k;
end
% least squares
c_lsq = (X'*X)\(X'*data);
% ridge regression
lambda = 1;
c_ridge = (X'*X + lambda*eye(p+1,p+1))\(X'*data);
c_tikhonov = (X'*X + lambda * Gimmel'*Gimmel)\(X'*data);

figure(1), clf
hold on 
plot(xd,Xd*c_lsq,'-b')
plot(xd,Xd*c_ridge,'-r')
plot(xd,Xd*c_tikhonov,'-m')
plot(x,data,'ko','MarkerFaceColor','k')
plot(x,y,'k--')
hold off
axis square
h = legend('least squares','ridge','Tikhonov','data','ground truth');
h.Location = 'SoutheastOutside';
axis([0 1 -inf inf])

figure(2),clf
subplot(1,3,1)
bar(0:length(c_lsq)-1,c_lsq,'FaceColor','b')
xlabel('c_{k;lsq}'), axis square, title('LSQ')

subplot(1,3,2)
bar(0:length(c_ridge)-1,c_ridge,'FaceColor','r')
xlabel('c_{k;ridge}'), axis square, title('ridge')

subplot(1,3,3)
bar(0:length(c_ridge)-1,1./gamma(1:length(c_ridge)),'FaceColor','k')
axis square, title('ground truth')

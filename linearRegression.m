set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
addpath(genpath('/home/user/Dropbox/Codes/fracPty/utils/export_fig2018'))
clf
% matlab code: linear regression

rng(666,'twister')
x = (-5:5)';
N = length(x);
d = 2*x + 4*(2*rand(N,1)-1);
% generate outlier
% d(10) = 20;
% lsq
aEstimate = x'*d / (x'*x);

plot(x, d, 'ko', 'MarkerFaceColor', 'k')
hold on
plot(x,aEstimate*x,'k--','LineWidth',2)
hold off, axis square, grid on, xlabel('x')
legend('data', 'fit','location','Northwest')
export_fig('linearRegression1.png')
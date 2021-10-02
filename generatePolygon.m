set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultFigureColor', 'w');
addpath(genpath('utils'))
clear
n = 11;
F = fft(eye(n,n));
figure(1)
plot(F,'k','lineWidth',1)
% im = getimage(gca);
axis square off
export_fig('../testImages/Polygon.jpg','-m2','-transparent','-silent')
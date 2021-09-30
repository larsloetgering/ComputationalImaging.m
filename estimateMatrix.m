set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultFigureColor', 'w');
addpath(genpath('utils'))
clear
close all
A = toeplitz(1:9);

X = rand(9);
B = A * X;

A_sol = B*X'/(X*X');

figure(1)
subplot(1,2,1)
imagesc(A)
axis image off
colorbar
title('ground truth')

subplot(1,2,2)
imagesc(A_sol)
axis image off, colormap jet
colorbar
title('estimated')
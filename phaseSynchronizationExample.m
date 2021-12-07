set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultFigureColor', 'w');
clear
addpath(genpath('utils'))

charge = 1;
N = 512;
x = linspace(-N/2,N/2,N);
[X,Y] = meshgrid(x);

% get aberration
[THETA,R] = cart2pol(X,Y);
[W,~] = zernike(R/(N/2),THETA,-1,3);

rng(0)
noiseAmplitude = 5e-2;
amplitude1 = exp(-(X.^2 + Y.^2) / (N/4)^2) + ...
            noiseAmplitude*randn(N,N); 
beam1 = amplitude1 .*... 
        exp(1i * (X.^2 + Y.^2)/1e3) .* ...
        exp(1i * charge * atan2(Y,X)) .* ...
        exp(1i*20*W);
    
amplitude2 = exp(-(2*X.^2+Y.^2) / (N/4)^2) + ...
            noiseAmplitude*randn(N,N); 
beam2 = amplitude2 .* ...
        exp(1i * (X.^2 + Y.^2)/1e3) .* ...
        exp(1i * charge * atan2(Y,X)) .* ...
        exp(1i*20*W)* exp(1i*pi);
    
exp_iPhi = phaseSynchronization(beam1, beam2);

figure(1)
subplot(1,3,1)
hsvplot(beam1)
title('a')
axis off
 
% figure(2)
subplot(1,3,2)
hsvplot(exp_iPhi*beam1)
title('a{\cdot}exp(i\phi)')
axis off

subplot(1,3,3)
hsvcolorbarplot(beam2, 'colorbar', 'test')
title('b')
axis off
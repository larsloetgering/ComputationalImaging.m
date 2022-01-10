% config
set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
addpath(genpath('utils'))
clear

% read images
g1 = imread('..\testImages\Dog.png');
g2 = imread('..\testImages\Windmill.png');
% exctract only red channel
g1 = single( g1(:,:,1) );
g2 = single( g2(:,:,1) );
% interpolate to same size
g1 = imresize(g1, [256 256],'bilinear');
g2 = imresize(g2, [256 256],'bilinear');
%% compute Fourier transform
G1 = fft2c(g1);
G2 = fft2c(g2);
%% extract modulus and phase
mod1 = abs(G1);
phi1 = angle(G1);
mod2 = abs(G2);
phi2 = angle(G2);
%% exchange modulus and phase
H1 = mod1 .* exp(1i * phi2);
H2 = mod2 .* exp(1i * phi1);
%% compute inverse Fourier transform
h1 = ifft2c(H1);
h2 = ifft2c(H2);
%% show results
% original image 1
figure(1)
imagesc(g1)
axis image off
colormap gray
% original image 2
figure(2)
imagesc(g2)
axis image off
colormap gray
% final image 1
figure(3)
imagesc(abs(h1))
axis image off
colormap gray
% final image 2
figure(4)
imagesc(abs(h2))
axis image off
colormap gray
addpath(genpath('utils'))
clc, clear
set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultFigureColor', 'w');

N = 2^9;
x = linspace(-1,1,N);
[X,Y] = meshgrid(x);
[theta,r] = cart2pol(X,Y);

c = 50;
% source support
A = exp(-c*(X.^2 + Y.^2)) - exp(-2*c*(X.^2 + Y.^2));
% normalization
A = A / max(A(:));
% target support
B = A;

lambda = 1e0;
rng(0)
v = exp(1i*2*pi*rand(N,N));

figure(1)
imagesc(A)
axis image off
colormap(jet)

%% power method
v = v / sqrt(v(:)'*v(:));
numIter = 1e4;
figureUpdate = 1000;

for k = 0:numIter
    v = linOp(v,B,A,lambda);
    v = v / sqrt(v(:)'*v(:));
    
    if mod(k,figureUpdate) == 0
        figure(2)
        hsvplot(v)
        zoom(4*sqrt(c/50))
        
        figure(3)
        hsvplot(fft2c(v))
        zoom(4*sqrt(c/50))
        title(['iteration ', num2str(k)])
        drawnow
    end
end

function r = linOp(v,B,C,lambda)
r = 1./(1+lambda-C) .* ifft2c(B.*fft2c(v));
end

function G = fft2c(g)
G = fftshift(fft2(ifftshift(g)))/numel(g);
end

function g = ifft2c(G)
g = fftshift(ifft2(ifftshift(G)))/numel(G);
end

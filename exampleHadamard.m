I = single(imread('cameraman.tif'));

[M, N] = size(I);
y = linspace(-M/2, M/2-1, M);
x = linspace(-N/2, N/2-1, N);

% normalize
I = I / sum(sum(I));

xs = sum(sum(I .* x));
ys = sum(sum(I .* y'));

figure(1)
imagesc(I)
axis image, colormap gray
%%


I = I/max(I(:));
I2 = 1-I.^2;
% reLU

figure(2)
subplot(1,2,1)
imagesc(I)
axis image

subplot(1,2,2)
imagesc(I2)
axis image

%% mutual information, NGF, SSIM, BM3D ...
% argmax_P(\lambda),O(\lambda)

figure(3)
plot(I(:), I2(:), 'o')

% daten schicken, code benchmarken, schneller, rekonstruktionsqualität?
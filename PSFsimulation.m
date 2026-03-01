%% configuration
fsize = 32;
set(0, 'DefaultAxesFontSize', fsize);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');

% Parameters
lambda = 550e-9; % Wavelength of light in meters (e.g., 550 nm for green light)
NA = 0.95; % Numerical aperture of the optical system
M = 200; % Magnification
f = 200e-3; % Focal length of the lens in meters (e.g., 200 mm)
pixel_size = 1e-2*6.5e-6; % Pixel size of the camera in meters (e.g., 6.5 um)
image_size = 1024; % Size of the image (image_size x image_size pixels)

% Calculate the spatial frequency coordinates
x = linspace(-image_size/2, image_size/2-1, image_size) * pixel_size / M;
y = linspace(-image_size/2, image_size/2-1, image_size) * pixel_size / M;
[X, Y] = meshgrid(x, y);
r = sqrt(X.^2 + Y.^2);

% Calculate the Airy disk pattern
k = 2 * pi / lambda; % Wave number
alpha = NA / f; % Aperture angle
PSF = (2 * besselj(1, k * alpha * r) ./ (k * alpha * r)).^2;
% PSF(r == 0) = 1; % Correct the center value
PSF(isnan(PSF(:))) = max(PSF(:));
% Normalize the PSF
PSF = PSF / max(PSF(:));

% Display the PSF
figure(1);
imagesc(x * 1e6, y * 1e6, nthroot(PSF+1e-3,3));
colormap(flipud(gray));
% colorbar;

dticks =  0.025;
xticks(-0.2:dticks:0.2)
yticks(-0.2:dticks:0.2)
axis image off;
grid on
% xlabel('$x (\mu m)$');
% ylabel('$y (\mu m)$');
% title('2D In-Focus Point Spread Function');

%% product of shifted PSFs
shift = 100;
PSF1 = circshift(PSF, [0 -shift]);
PSF2 = circshift(PSF, [0 +shift]);
productPSF = PSF1 .* PSF2;
productPSF = productPSF / max(productPSF(:));

figure(2); clf
plot(x,(PSF1(end/2,:)),'k','linewidth',2)
hold on
plot(x,(PSF2(end/2,:)),'k','linewidth',2)
plot(x,(productPSF(end/2,:)),'k--','linewidth',2)
xticks('')
yticks([0:0.5:1])
grid on
axis tight square
h = gca;
h.GridLineWidth = 3;
legend('$PSF_i$','$PSF_d$','$\prod_{k=i,d} PSF_k$')




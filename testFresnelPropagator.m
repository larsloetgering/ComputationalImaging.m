set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
addpath('utils')

% source grid (primed)
N = 2^9;
dxp = 2e-6;
L = N*dxp;
x = linspace(-L/2,L/2,N);
[X,Y] = meshgrid(x);
wavelength = 633e-9;

% define aperture
w = 200e-6;
aperture = rect(X/(L/2)).*rect(Y/(L/2));
% aperture = circ(X,Y,w);
% aperture = imfilter(aperture,...
%             fspecial('disk',2),'same');

figure(1)
h = subplot(2,2,1);
imagesc(x,x,aperture)
xlabel('x / m'),ylabel('y / m')
title('aperture')
axis image, colormap gray

%%
z = 2e-2;
[aperture_propagated, Qin, Qout] = ...
fresnelPropagator(aperture, z, wavelength, L);

subplot(2,2,2)
imagesc(sqrt(abs(aperture_propagated)))
axis image
title('propagated')
subplot(2,2,3)
hsvplot(Qin)
title('Q_{in}')

subplot(2,2,4)
hsvplot(Qout)
title('Q_{out}')
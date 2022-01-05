clear
set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
addpath('utils')

N = 2^9;
dx = 5e-6;
L = N*dx;
x = linspace(-L/2,L/2,N);
[X,Y] = meshgrid(x);
wavelength = 633e-9;

% define aperture
w = 200e-6;
% aperture = rect(X/(L/2)).*rect(Y/(L/2));
aperture = circ(X,Y,w);
aperture = imfilter(aperture,...
            fspecial('disk',2),'same');

figure(1)
h = subplot(1,3,1);
imagesc(x,x,aperture)
xlabel('x / m'),ylabel('y / m')
title('aperture')
axis image, colormap gray
%%
M = 50;
z = linspace(0,10e-2,M);
dz = z(2) - z(1);

aperture_diffracted = zeros(N,N,M);
tic
alpha =0.6;
for k = 1:length(z)
    aperture_diffracted(:,:,k) = ...
        exp(-1i*2*pi/wavelength*z(k))*...
        aspw(aperture, z(k),wavelength,L);
    
    subplot(1,3,2)
    hsvplot_axis(z,x,...
        squeeze(aperture_diffracted(:,N/2,:)),...
        'intensityScale', [0,alpha])
    axis square, title('y-z cross section')
    
    subplot(1,3,3)
    hsvplot_axis(x,x,...
        aperture_diffracted(:,:,k),...
        'intensityScale', [0,alpha])
    title(['z = ', num2str(round(z(k)*1e3)),' mm'])
    drawnow
end
% reset colormap of panel 1
h.Colormap = gray;
toc
%% configuration
fsize = 32;
set(0, 'DefaultAxesFontSize', fsize);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');
addpath('utils')

% Parameters
lambda = 500e-9; % emission wavelength
NA = 0.5;        % Numerical aperture of the optical system
dx = 0.125e-6;   % lateral pixel size of the camera in meters
dz = 0.25e-6;    % axial pixel size
N = 1024;        % Size of the image (image_size x image_size pixels)
K = 150;

%% define 3D FFTs

FT3 = @(x) fftshift(fftn(ifftshift(x))) / sqrt(numel( x ) );
iFT3 = @(x) fftshift(ifftn(ifftshift(x))) * sqrt(numel( x ) );

%% get spatial frequencies

fx = (-N/2:N/2-1)/(N*dx);
fy = (-N/2:N/2-1)/(N*dx)';
fperp = sqrt(fx.^2 + fy'.^2);
fz = reshape((-(K-1)/2:(K-1)/2)/(K*dz), [1,1,K]);

[Fx, Fy] = meshgrid(fx,fy);

pupil = circ(Fx,Fy,2*NA/lambda);

figure(1)
imagesc(pupil)
axis image off
colormap gray

%% simulate defocus CSF

z = reshape( (-(K-1)/2:(K-1)/2)*dz, [1,1,K]);
csf = ifft2c( pupil .* exp(1i * 2*pi/lambda * z .* real(sqrt(1-(lambda*fperp).^2))) );

%% get PSF

psf = abs(csf).^2;

figure(2)
imagesc(psf(:,:,21))
axis square
colormap gray
zoom(4)

%% get OTF

% numerical approach
OTF_wf = FT3(psf);
OTF_conf = FT3(psf.^2);

%%
p = 4;
figure(3)
subplot(1,2,1)
imagesc(nthroot(abs(squeeze(OTF_wf(:,end/2+1,:))),p), [1 inf])
axis square off
colormap gray

subplot(1,2,2)
imagesc( imfilter(nthroot(abs(squeeze(OTF_conf(:,end/2+1,:))),p), ones(15)), [100 inf])
axis square off
colormap gray

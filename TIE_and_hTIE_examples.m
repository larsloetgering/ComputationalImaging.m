clear
% appearance of figures
addpath('utils')
set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')

% path to data
folder = '../datasets/';
file = 'Tonsil_1um_TIE_BF_AF4_smallROI.czi';
stack = CZIzstackReader([folder,file]);

%% parameters

pixel_size = 172e-9;
wavelength = 500e-9;
dz = 1e-6;          % axial slice separation
epsilon = 0.01;     % controls merge point for hTIE
zoom_fac = 1;       % zomm factor for figures

%% compute standard tie at short and large axial distance

phi = tie(stack(:,:,5:7), 1e-6, 500e-9, 172e-9, 1e1, 100);
phi_large_reg = tie(stack(:,:,5:7), 1e-6, 500e-9, 172e-9, 1e3, 100);
phi_large_dist = tie(stack(:,:,1:5:11), 4e-6, 850e-9, 172e-9, 1e1, 1000);

figure(11)
imagesc(var_stab(phi), 3*[-1 1])
axis image off
colormap gray
zoom(zoom_fac)

figure(12)
imagesc(var_stab(phi_large_reg), 3*[-1 1])
axis image off
colormap gray
zoom(zoom_fac)

figure(13), clf
imagesc(var_stab(phi_large_dist), 3*[-1 1])
axis image off
colormap gray
zoom(zoom_fac)

%% example hTIE

mergeFreq = sqrt(sqrt(6*epsilon)/(pi*wavelength*dz));

[Ny, Nx] = size(stack, [1,2]);

% Create frequency coordinates
% Frequency spacing
dfx = 1 / (Nx * pixel_size);
dfy = 1 / (Ny * pixel_size);
fx = (-Nx/2 : Nx/2-1) * dfx;
fy = (-Ny/2 : Ny/2-1)' * dfy;

W = exp(-log(2)*(fx.^2 + fy.^2)/(mergeFreq)^2);
idx = find( abs( W(:)-0.5 ) < 0.1);
PHI_low = fft2c(phi_large_dist);
PHI_high = fft2c(phi);
hilo = real(ifft2c( W.* PHI_low + (1-W).*PHI_high ));

figure(14)
imagesc(var_stab(hilo), 3*[-1 1])
colormap gray
axis image off
zoom(zoom_fac)
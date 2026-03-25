%% configuration
clear
fsize = 32;
set(0, 'DefaultAxesFontSize', fsize);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');
lw = 2;
addpath('utils')

%%  function definitions

FT2 =  @(x) fftshift(fft2(ifftshift(x)));
iFT2 =  @(x) fftshift(ifft2(ifftshift(x)));

%% mirror padding illustration

N = 512;
x = (-N/2:N/2-1)/(N/2);
[X,Y] = meshgrid(x);
[theta, rho] = cart2pol(X,Y);

t = ((1/2*sin(0*rho + 7*theta))>0) .* (1-exp(-(X.^2+Y.^2)/(2*(0.05).^2) ));
ckernel = exp(-(X.^2+Y.^2)/(2*(0.1).^2) ) - exp(-(X.^2+Y.^2)/(2*(0.05).^2) );
ckernel = ckernel / max(ckernel(:));
tresult = real(iFT2(FT2(t).*FT2(ckernel)));

figure(4)
imagesc(t,[0,1])
axis image off
colormap gray

figure(5)
imagesc(ckernel,[0,1])
axis image off
colormap gray

figure(6)
imagesc(tresult)
axis image off
colormap gray


%%
tp = mirror_pad(t,100);
ckernelp = mirror_pad(ckernel,100);

tresultp = real(iFT2(FT2(tp).*FT2(ckernelp)));

figure(7)
imagesc(tp,[0,1])
axis image off
colormap gray

figure(8)
imagesc(ckernelp,[0,1])
axis image off
colormap gray

figure(9)
imagesc(crop(tresultp,100))
axis image off
colormap gray

%% illustrate Wiener DCV with generalized cross-validation

% step 1: generate noisy observation
bit_depth = 8;
I = tresultp / max(tresultp(:)) * 2^bit_depth;
I(:) = poissrnd( I(:) );

figure(10)
imagesc(crop(I,100))
axis image off
colormap gray

%% try Wiener filter on noisy observation

lambda = 8.685114e+03;
H = FT2(ckernelp);
wiener_dcv = iFT2( FT2(I) .* conj(H) ./ (abs(H).^2 + lambda) );

figure(11)
imagesc(crop( wiener_dcv, 100))
axis image off
colormap gray

%%

figure(12)
imagesc(crop( g_optimal, 100))
axis image off
colormap gray

%%

[g_optimal, p_optimal, gcv_values, p_values] = ...
    gcv_deconvolution(I, ckernelp,'p_min', 1e1, 'p_max', 1e7,...
    'reg_type', 'gradient');


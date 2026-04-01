clear
fsize = 32;
set(0, 'DefaultAxesFontSize', fsize);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');
lw = 2;

%% parameters for optimization

maxIter = 1e3; % max number of iterations

%% 2D Deconvolution Example

% Create synthetic image
fileName = '../datasets/nuclei_DAPI.tif';
img = imread(fileName);
img = single(rgb2gray(img));
img = img / max(img(:));

[M,N] = size(img);
x = (-N/2:N/2-1)/(N/2);
y = (-M/2:M/2-1)'/(M/2);
[X,Y] = meshgrid(x);
[~, rho] = cart2pol(X,Y);

% Create Gaussian PSF
s = 0.01;
psf = exp(-(X.^2+Y.^2)/(2*(s).^2) ) - exp(-(X.^2+Y.^2)/(2*(1*s/4).^2) );
psf = psf / max(psf(:));

figure(3)
subplot(1,2,1)
imagesc(psf)
axis image off
colormap gray

% Simulate blur and additional Poisson noise
bit_depth = 2^12;
blurred = imfilter(img, psf, 'conv', 'circular');
blurred = blurred / max(blurred(:));
blurred(:) = poissrnd(blurred(:) * bit_depth);
blurred = blurred / max(blurred(:));

subplot(1,2,2)
imagesc(blurred)
axis image off
colormap gray

%% generalized cross validation
% uncoment to estimate regularization parameters

% [g_optimal, p_optimal, gcv_values, p_values] = gcv_deconvolution(blurred, psf,'reg_type','tikhonov');
% [g_optimal, p_optimal, gcv_values, p_values] = gcv_deconvolution(blurred, psf,'reg_type','gradient');

%% Wiener filter (equivalent to L2 data term + Tikhonov prior)
% regularizers = {'tikhonov', 'tv', ...
% 'goods', 'hessian'};
% dataTerms = {'l2', 'l1', 'poisson'};

lr = 1e2;
deconv = Deconvolver(psf, 'Regularizer', 'tikhonov', ...
    'DataTerm', 'l2', 'Lambda', 7.5e-1);

% Perform deconvolution
restored = deconv.deconvolve(blurred, 'Iterations', maxIter, ...
    'LearningRate', lr, 'Verbose', true);

figure(2), clf
subplot(2,3,1)
imagesc(img)
axis image off
title('ground truth')
colormap gray

subplot(2,3,2)
imagesc(blurred)
axis image off
title('raw data')

subplot(2,3,3)
imagesc(real(restored))
axis image off
title('Wiener')

%% Good's roughness

lr = 1e2;
% Goods filter
deconv = Deconvolver(psf, 'Regularizer', 'goods', ...
    'DataTerm', 'l2', 'Lambda', 2.2e0);

% Perform deconvolution
restored = deconv.deconvolve(blurred, 'Iterations', maxIter, ...
    'LearningRate', lr, 'Verbose', true);

figure(2)
subplot(2,3,4)
imagesc(real(restored))
axis image off
title('Good`s')

%% tv

lr = 1e2;
deconv = Deconvolver(psf, 'Regularizer', 'tv', ...
    'DataTerm', 'l2', 'Lambda', 1e-1);

% Perform deconvolution
restored = deconv.deconvolve(blurred, 'Iterations', maxIter, ...
    'LearningRate', lr, 'Verbose', true);

figure(2)
subplot(2,3,5)
imagesc(max(real(restored),0))
axis image off
title('TV')

%% Hessian

lr = 1e2;
deconv = Deconvolver(psf, 'Regularizer', 'hessian', ...
    'DataTerm', 'l2', 'Lambda', 1e-1);

% Perform deconvolution
restored = deconv.deconvolve(blurred, 'Iterations', maxIter, ...
    'LearningRate', lr, 'Verbose', true);

figure(2)
subplot(2,3,6)
imagesc(max(real(restored),0))
axis image off
title('Hessian')
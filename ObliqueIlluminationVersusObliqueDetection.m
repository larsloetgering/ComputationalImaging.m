fsize = 23;
set(0, 'DefaultAxesFontSize', fsize);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');
addpath(genpath('utils'))
% compute transfer functions

N = 1000;
x = (-(N-mod(N,2))/2:(N+mod(N,2))/2-1);
[X,Y] = meshgrid(x);
rho = sqrt(X.^2 + Y.^2);
theta = atan2(Y,X);

source = zeros(N,N,2);
source(:,:,1) = circ(X,Y,N/2-1) .* (X>0);  % oblique illumination
source(:,:,2) = circ(X,Y,N/2-1);  % cicular illumination

pupil = source(:,:,[2,1]);

show_index = 2;
figure(1)
imagesc(source(:,:,show_index))
axis image off
colormap gray

figure(2), clf
imagesc(pupil(:,:,show_index))
axis image off
colormap gray

% compute ATF and PTF
PTF = 0*source;
for k = 1:2
    [PTF_temp, ~] = getTransferFunctions(...
        imfilter(source(:,:,k), fspecial('Gaussian', 8, 1),'same'), ...
        imfilter(pupil(:,:,k), fspecial('Gaussian', 8, 1),'same'));
    PTF(:,:,k) = PTF_temp;
end

figure(3), clf
imagesc(imag(PTF(:,:,show_index)))
axis image off
colormap(slanCM('RdGy'))
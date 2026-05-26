fsize = 23;
set(0, 'DefaultAxesFontSize', fsize);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');
addpath('utils')
% compute transfer functions

N = 1000;
x = (-(N-mod(N,2))/2:(N+mod(N,2))/2-1);
[X,Y] = meshgrid(x);

source = zeros(N,N,3);
source(:,:,1) = circ(X,Y,N/20);    % quasi-coherent
source(:,:,2) = circ(X,Y,N/2-1); % standard bright field case (matched NA bright field)
source(:,:,3) = circ(X,Y,N/2-1) - circ(X,Y,N/2-15); % ring shaped

pupil = circ(X,Y,N/2-1);

show_index = 3;
figure(1)
imagesc(source(:,:,show_index))
axis image off
colormap gray

figure(2)
imagesc(pupil)
axis image off
colormap gray

%% compute ATF
ATF = 0*source;
for k = 1:3
    [~, ATF_temp] = getTransferFunctions(imfilter(source(:,:,k), fspecial('Gaussian', 8, 1),'same'), ...
        imfilter(pupil, fspecial('Gaussian', 8, 1),'same'));
    ATF(:,:,k) = ATF_temp;
end

show_index = 1;
figure(3)
imagesc(real(ATF(:,:,show_index)))
axis image off
colormap(gray)

%% show cross-sections
markerTypes = {'k-','k--','k-*'};

figure(4)
hold on
for k = 1:3
plot(real(ATF(end/2+1,:,k)),markerTypes{k},'lineWidth', 2)
axis square off
end
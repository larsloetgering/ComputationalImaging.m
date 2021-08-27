set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultFigureColor', 'w');
% create coordinate meshgrid
N = 2^7; dx = 100e-9;
L = N * dx;
% x = linspace(-N/2,N/2,N)*dx; z = x;
x = (-N/2:N/2-1)*dx; z = x;
[X,Z] = meshgrid(x);
% determine wavevector
k = 2*pi/633e-9;

% source coordinate
xs = 0;
zs = -1.25*L;
% displacement meshgrid

% example 1: Uin = spherical wave
R = sqrt((X-xs).^2 + (Z-zs).^2);
Uin = exp(1i*k*R)./(1i*k*R);

% example 2: Uin = plane wave
% Uin = exp(1i*k*Z);

figure(1)
hsvplot(Uin)
axis image off
% imagesc(x,x,abs(Uin)), axis image, colormap jet
title('U_{in}')

%% determine random source points

cmap = jet(50);
aleph = linspace(0,1,10);
cmap = [aleph' * cmap(1,:);cmap];
rng(4)
idx = 1:10:N^2; numScatterers = length(idx); idx = min(idx+randi([1 4],size(idx)),N^2);
aleph = 5e-2*ones(numScatterers,1);
source = zeros(N,N);
source(idx) = aleph;
R = sqrt(X.^2 + Z.^2);
source = source .* (R < (N/4)*dx);
idx = find(source>0); numScatterers = length(idx);

figure(2)
imagesc(x,z,source)
axis image
colormap gray

%% normalization factor

G = mean(aleph)*exp(1i * k * R) ./ (1i * k * R + eps);

% method 1: Foldy Lax
M = zeros(N^2,N^2);
for loop = 1:length(idx)
    R = sqrt((X-X(idx(loop))).^2 + (Z-Z(idx(loop))).^2);
    % 2D Green's function: Hankel
    G = aleph(loop) * 1i*pi*besselh(0,k * R);
    G(idx(loop)) = 0;  
    M(:,idx(loop)) = G(:);
end

tic
Ufoldy = reshape( sparse(eye(N^2) - M) \ Uin(:), [N,N]);
toc
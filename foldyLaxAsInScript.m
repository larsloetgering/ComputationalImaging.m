set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultFigureColor', 'w');
%addpath(genpath('C:\Users\Lars Loetgering\Dropbox\Codes\fracPty\utils\export_fig2018'))

% create coordinate meshgrid
N = 2^6; dx = 100e-9; L = N * dx;
x = (-N/2:N/2-1)*dx; z = x;
[X,Z] = meshgrid(x);

% determine wavevector
k = 2*pi/633e-9;

% source coordinate
xs = 0;zs = -1.25*L;

% spherical wave 
R = sqrt((X-xs).^2 + (Z-zs).^2);
Uin = exp(1i*k*R)./(1i*k*R);

%% determine random point scatter distribution

rng(0,'twister')
idx = 1:10:N^2; 
numScatterers = length(idx); 
idx = min(idx+randi([1 4],size(idx)),N^2);
aleph = 5e-2;
pointScatterDist = zeros(N,N);
pointScatterDist(idx) = aleph;
pointScatterDist = pointScatterDist .* (X.^2+Z.^2 < (L/4)^2);
idx = find(pointScatterDist(:)>0);

%% compute scattering matrix

G = mean(aleph)*exp(1i * k * R) ./ (1i * k * R + eps);

% method 1: Foldy Lax
M = zeros(N^2,N^2);
for loop = 1:length(idx)
    R = sqrt((X-X(idx(loop))).^2 + (Z-Z(idx(loop))).^2);
    % 2D Green's function: Hankel
    G = aleph * 1i*pi*besselh(0,k * R);
    G(idx(loop)) = 0;  
    M(:,idx(loop)) = G(:);
end

Ufoldy = reshape( sparse(eye(N^2) - M) \ Uin(:), [N,N]);

% method 2: 1st Born approximation
Uborn = Uin(:) + M*Uin(:);
Uborn = reshape( Uborn, [N,N]);

maxAbs = max(max(abs(Ufoldy(:)),abs(Uborn(:))));

%% show results
% modify colormap
cmap = jet(50);
c = linspace(0,1,10);
cmap = [c' * cmap(1,:);cmap];

figure(3)
subplot(1,4,1)
imagesc(x,x,abs(Uin),[0 maxAbs]), axis image, colormap(cmap)
title({'U_{in}',''}), colorbar

subplot(1,4,2)
imagesc(x,x,pointScatterDist), axis image, colormap(gray)
h = gca; title({'point scatter distr.',''}), colorbar

subplot(1,4,3)
imagesc(x,x,abs(Ufoldy),[0 maxAbs]), axis image, colormap(cmap)
title({'U_{foldy}',''}), colorbar

subplot(1,4,4)
imagesc(x,x,abs(Uborn),[0 maxAbs]), axis image, colormap(cmap)
title({'U_{born}',''}), colorbar

colormap(h, gray)

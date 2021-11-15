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
xs = 0;zs = -1.1*L;
% xs = 0;zs = L/2;

% spherical wave 
R = sqrt((X-xs).^2 + (Z-zs).^2);
Uin = exp(1i*k*R)./(1i*k*R);
Uin(isnan(Uin)) = 0;

% plane wave
% Uin = exp(1i*k*Z);

%% determine random point scatter distribution

rng(0,'twister')
% idx = 1:10:N^2;
idx = 1:10:N^2;
numScatterers = length(idx); 
idx = min(idx+1*randi([1 4],size(idx)),N^2);
aleph = 10e-2;
pointScatterDist = zeros(N,N);
pointScatterDist(idx) = aleph;
pointScatterDist = pointScatterDist .* (X.^2+Z.^2 < (L/4)^2);
% pointScatterDist = (X.^2+4*Z.^2 < (L/4)^2);
idx = find(pointScatterDist(:)>0);

%% compute scattering matrix
G = exp(1i * k * R) ./ (1i * k * R + eps);
% method 1: Foldy Lax
M = zeros(N^2,N^2);
for loop = 1:length(idx)
    R = sqrt((X-X(idx(loop))).^2 + (Z-Z(idx(loop))).^2);
    % 2D Green's function: Hankel
    G = 1i*pi*besselh(0,k * R);
    G(idx(loop)) = 0;  
    M(:,idx(loop)) = G(:);
end
M = aleph * M;
Ufoldy = reshape( sparse(eye(N^2) - M) \ Uin(:), [N,N]);

% method 2: 1st Born approximation
Uborn = Uin(:) + M*Uin(:);
Uborn = reshape( Uborn, [N,N]);

maxAbs = max(max(abs(Ufoldy(:)),abs(Uborn(:))));

%% higher order Born
% 
k = 5;
UbornHigher = Uin(:);
% gamma = pointScatterDist/1;
gamma = 75e-2; % if not convergent, reduce this parameter
Mborn = sparse( gamma*M + (1-gamma)*eye(N^2) );
for loop = 1:k
    UbornHigher(:) = gamma * Uin(:) + Mborn * UbornHigher(:);
    
    if mod(loop,k/10) == 0
    figure(4)
    imagesc(x,x,abs(reshape( UbornHigher, [N,N])),[0 maxAbs]), axis image, colormap(jet)
    title({'U_{born}',['iteration: ', num2str(loop)]}), colorbar
    end
end
UbornHigher = reshape( UbornHigher, [N,N]);

lambda_max = abs(eigs(M,1));
disp(['largest eigenvalue of M: ', num2str(lambda_max)])
disp(['maximum gamma: ', num2str(1/lambda_max)])

%% show results
% modify colormap
cmap = jet(50);
c = linspace(0,1,10);
cmap = [c' * cmap(1,:);cmap];

figure(3)
subplot(1,5,1)
imagesc(x,x,abs(Uin),[0 maxAbs]), axis image, colormap(cmap)
title({'U_{in}',''}), colorbar

subplot(1,5,2)
imagesc(x,x,pointScatterDist), axis image, colormap(gray)
h = gca; title({'point scatter distr.',''}), colorbar

subplot(1,5,3)
imagesc(x,x,abs(Ufoldy),[0 maxAbs]), axis image, colormap(cmap)
title({'U_{foldy}',''}), colorbar

subplot(1,5,4)
imagesc(x,x,abs(Uborn),[0 maxAbs]), axis image, colormap(cmap)
title({'U_{born}',''}), colorbar

subplot(1,5,5)
imagesc(x,x,abs(UbornHigher),[0 maxAbs]), axis image, colormap(jet)
title({'U_{born,higher}',''}), colorbar
colormap(h, gray)
set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultFigureColor', 'w');
% create coordinate meshgrid
N = 2^6; dx = 100e-9;
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
% imagesc(x,x,abs(Uin)), axis image, colormap jet
title('U_{in}')

%% determine random source points

cmap = jet(50);
aleph = linspace(0,1,10);
cmap = [aleph' * cmap(1,:);cmap];
rng(4)
% idx = [N^2/2+N/2,N^2/2+N/2-10,N^2/2+N/2+10*N]; aleph = ones(3,1);
% numScatterers = 5; idx = randi([1,N^2],numScatterers,1);aleph = 1e-6*rand(numScatterers,1);

% numScatterers = 40; idx = randi([1,N^2],numScatterers,1);aleph = 2e-8*ones(numScatterers,1);
% numScatterers = 10; idx = randi([1,N^2],numScatterers,1);aleph = 1e-1*ones(numScatterers,1);
% numScatterers = 40; idx = N^2/2+N/4 + randi([1,N/2],numScatterers,1);aleph = 20e-8*rand(numScatterers,1);
% numScatterers = 5; idx = N^2/2+N/4 + randi([1,20],numScatterers,1);aleph = 20e-8*ones(numScatterers,1);
% create image of source distribution
% idx = 1:100:N^2; numScatterers = length(idx); aleph = 1e1*ones(numScatterers,1);
% idx = 1:100:N^2; numScatterers = length(idx);
idx = 1:10:N^2; numScatterers = length(idx); idx = min(idx+randi([1 5],size(idx)),N^2);
% idx = 1:10:N^2; numScatterers = length(idx);
aleph = 10e-1*ones(numScatterers,1);
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
normFactor = sqrt(G(:)' * G(:));
fillValue = mean(aleph) * (G(end/2,end/2+1) + G(end/2+2,end/2+1) + G(end/2+1,end/2) + G(end/2+1,end/2+2))/4; 

% method 1: Foldy Lax
M = zeros(N^2,N^2);
for loop = 1:length(idx)
    R = sqrt((X-X(idx(loop))).^2 + (Z-Z(idx(loop))).^2);
    G = aleph(loop) * exp(1i * k * R) ./ (1i * k * R + eps);
    G(idx(loop)) = fillValue;  
%     G(idx(loop)) = 0;  
%     G = G / normFactor;
    M(:,idx(loop)) = G(:);
end

tic
Ufoldy = reshape( (eye(N^2) - M) \ Uin(:), [N,N]);
toc
%%
% 

% %
% figure(3)
% subplot(2,2,1)
% imagesc(x,x,abs(Ufoldy)), axis image, colormap(cmap)
% title('U_{out}')
% subplot(2,2,2)
% hsvplot(x,z,Ufoldy-Uin)
% title('U_{out}-U_{in}')
% 
% % Born approximation
% 
% Uborn = Uin(:) + M*Uin(:);
% Uborn = reshape( Uborn, [N,N]);
% 
% subplot(2,2,3)
% imagesc(x,x,abs(Uborn)), axis image, colormap(cmap)
% title('U_{born}')
% 
% subplot(2,2,4)
% hsvplot(x,z,phaseSynchronization(Ufoldy-Uin,Uborn-Uin))
% title('U_{born}-U_{in}')

%%

% Born approximation
Uborn = Uin(:) + M*Uin(:);
Uborn = reshape( Uborn, [N,N]);

maxAbs = max(max(abs(Ufoldy(:)),abs(Uborn(:))));

figure(3)
subplot(1,4,1)
imagesc(x,x,abs(Uin),[0 maxAbs]), axis image, colormap(cmap)
title('U_{in}')
colorbar

subplot(1,4,2)
imagesc(x,x,abs(Ufoldy),[0 maxAbs]), axis image, colormap(cmap)
title('U_{foldy}')
colorbar

subplot(1,4,3)
imagesc(x,x,abs(Uborn),[0 maxAbs]), axis image, colormap(cmap)
title('U_{born}')
colorbar

subplot(1,4,4)
imagesc(x,x,abs(Ufoldy-Uborn),[0 maxAbs]), axis image, colormap(cmap)
title('\Delta U')
colorbar

%% compute first two eigenfunctions

% [V,D] = eigs(eye(N^2,N^2)+M,2);
% 
% figure(5)
% subplot(1,2,1)
% imagesc(x,x,abs(reshape(V(:,1),[N,N]))), axis image, colormap jet
% 
% subplot(1,2,2)
% imagesc(x,x,abs(reshape(V(:,2),[N,N]))), axis image, colormap jet

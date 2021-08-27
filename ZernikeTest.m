set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
%

N = 2^8;
x = linspace(-1,1,N);
dx = x(2) - x(1);
[X,Y] = meshgrid(x);
[THETA, R] = cart2pol(X,Y);

%% calculate Zernike for m,n
m = 2;
n = 2;
Z22 = zernike(R,THETA,m,n) .* circ(X,Y,2);

m = 1;
n = 2; % n >= m
Z12 = zernike(R,THETA,m,n) .* circ(X,Y,2);
% check orthonormality
sum(Z22(:) .* Z22(:)) * dx.^2
sum(Z12(:) .* Z22(:)) * dx.^2

% renormalize
Z22 = Z22 / max(abs(Z22(:)));
Z12 = Z12 / max(abs(Z12(:)));
figure(1)
subplot(1,2,1)
imagesc(x,x,Z12)
xticks([-1 0 1]), yticks([-1 0 1])
axis image
h = colorbar;
h.Ticks = [min(Z12(:)), max(Z12(:))];
h.TickLabels = {'min', 'max'};

subplot(1,2,2)
imagesc(x,x,Z22)
xticks([-1 0 1]), yticks([-1 0 1])
axis image
cmap = [[1 0 0];...
        [1 1 1];...
        [0 0 1]];
cmap = imresize(cmap, [100, 3],'bilinear');
cmap(end/2,:) = [1,1,1];
cmap(end/2 + 1,:) = [1,1,1];
colormap(cmap)
h = colorbar;
h.Ticks = [min(Z22(:)), max(Z22(:))];
h.TickLabels = {'min', 'max'};
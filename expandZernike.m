set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
addpath(genpath('/home/user/Dropbox/Codes/fracPty/utils/export_fig2018'))
% matlab code: expandZernike.m

N = 2^6;
x = linspace(-1,1,N);
dx = x(2) - x(1);
[X,Y] = meshgrid(x);
% convert Cartesian to polar coordinates
[THETA, R] = cart2pol(X,Y);

%% calculate Zernike for m,n
cmap = [[1 0 0];...
        [1 1 1];...
        [0 0 1]];
cmap = imresize(cmap, [100, 3],'bilinear');
% function to be expanded
f =  cos(THETA - 2*pi * (X.^2 + Y.^2)) .* (R<=1);
Z = [];
fExpansion = 0;
for n = 0:10
    for m = -n:2:n
        disp(m)
        disp(n)
        [temp, J] = zernike(R,THETA,m,n);
        c = sum(sum(temp .* f)) * dx^2;
        fExpansion = fExpansion + c * temp;
        
        figure(1)
        subplot(1,3,1)
        imagesc(x,x,temp)
        xticks([-1 0 1]), yticks([-1 0 1])
        axis image off
        colormap(cmap)
        title({['n = ',num2str(n)], ['m = ',num2str(m)],...
               ['J = ',num2str(J)]})
        
        subplot(1,3,2)
        imagesc(f)
        axis image off
        title('f')
        
        subplot(1,3,3)
        imagesc(fExpansion)
        axis image off
        title('approx.')
    end
end

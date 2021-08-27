set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
addpath(genpath('/home/user/Dropbox/Codes/fracPty/utils/export_fig2018'))
% matlab code: showZernike.m

N = 2^8;
x = linspace(-1,1,N);
dx = x(2) - x(1);
[X,Y] = meshgrid(x);
[THETA, R] = cart2pol(X,Y);

%% calculate Zernike for m,n
cmap = [[1 0 0];...
        [1 1 1];...
        [0 0 1]];
cmap = imresize(cmap, [100, 3],'bilinear');
f =  (exp(1i*THETA).*exp(-1i * 2*pi * (X.^2 + Y.^2)) + 1).^2;
f = real(f-2) .* (R<=1);
f = (R<=1) .*f+ (1-(R<=1)) * (max(f(:).*(R(:)<=1))+min(f(:).*(R(:)<=1)))/2;

figure(1)
imagesc(f)
axis image off
colormap(cmap)
%%
Z = [];
fExpansion = 0;
for n = 0:17
    for m = -n:2:n
        
        [temp, J] = zernike(R,THETA,m,n);
        c = sum(sum(temp .* f)) *dx^2; % / sum(sum(temp .* temp))
        fExpansion = fExpansion + c * temp;
%         Z = cat(3,Z, temp .* circ(X,Y,2));
        
        figure(1)
        subplot(1,3,1)
        imagesc(x,x,temp)
        xticks([-1 0 1]), yticks([-1 0 1])
        axis image off
        colormap(cmap)
        title({['n = ',num2str(n)], ['m = ',num2str(m)], ['J = ',num2str(J)]})
        
        subplot(1,3,2)
        imagesc(f.* (R<=1))
        axis image off
        
        subplot(1,3,3)
        imagesc(fExpansion.* (R<=1))
        axis image off
        
        figure(4)
        imagesc(fExpansion.* (R<=1) + (1-(R<=1)) * (max(fExpansion(:).*(R(:)<=1))+min(fExpansion(:).*(R(:)<=1)))/2)
        axis image off
        colormap(cmap)
        if ismember(J, [10 20 50 100 200])
            export_fig([num2str(J),'.png'], '-transparent')
        end
    end
end

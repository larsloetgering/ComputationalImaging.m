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
cmap = cat(1,...
       cat(2,linspace(1,1,50),linspace(1,0,50)),... %R
       cat(2,linspace(0,1,50),linspace(1,0,50)),... %G
       cat(2,linspace(0,1,50),linspace(1,1,50)))';  %B

Z = [];
for n = 0:5
    for m = -n:2:n
        [temp, J] = zernike(R,THETA,m,n);
        Z = cat(3,Z, temp);
        figure%(1)
        imagesc(x,x,Z(:,:,end))  % scaled image display
        xticks([-1 0 1]), yticks([-1 0 1])
        axis image off
        colormap(cmap)        
        title(['n = ',num2str(n),...
               ', m = ',num2str(m), ', J = ',num2str(J)])
        pause(1)
    end
end

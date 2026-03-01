set(0, 'DefaultAxesFontSize', 64);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
set(0, 'DefaultLineLineWidth', 2);

addpath(genpath('utils'))
%% compute physical coordinate grid

M = 2048;
N = M;
wavelength = 500e-9;
dx = wavelength / 2;
L = N*dx;
x = (-N/2:(N/2-1))*dx;
[X,Y] = meshgrid(x);
[Theta,R] = cart2pol(Y,X);

%% generate source

% point souce (unmatched)
% source = circ( X,Y+(floor(N/8)-0)*dx,20*dx);

% point souce (matched)
% source = circ( X,Y+(floor(N/4)-0)*dx,20*dx);

% upper half plane
% source = circ( X,Y+0*(floor(N/4))*dx,N/2*dx) .* (Y<0);

% half ring
% source = circ( X,Y+0*(floor(N/4))*dx,N/2*dx) .* (Y<0) - ...
%          circ( X,Y+0*(floor(N/4))*dx,(N/2-20)*dx) .* (Y<0);

% azimuthal cosine
source = (1+cos(Theta+pi))/2;
source = circ( X,Y+0*(floor(N/4))*dx,N/2*dx) .* source - ...
         circ( X,Y+0*(floor(N/4))*dx,(N/2-20)*dx) .* source;

% generate transfer function
pupil = circ(X,Y,(N/2+1)*dx);
% pupil_grad = circ(X,Y,(N/2+10)*dx) - circ(X,Y,(N/2-10)*dx);

% get transfer function
[PTF, ~] = getTransferFunctions(source, pupil);
% [PTF2, ~] = getTransferFunctions(j*rot90(source,1), pupil);

figure(2), clf
subplot(1,3,1)
hsvxplot( PTF/max(abs(PTF(:))),...
    'intensityScale', [0 1])
% hsvxplot( PTF/max(abs(PTF(:))),...
%     'intensityScale', [0 1],'colorbar','test')
h = gca; h.FontSize = 80;
PTF(isnan(PTF(:))) = 0;
subplot(1,3,2)
% figure(99), clf
% hsvxplot( PTF + PTF2,...
%     'intensityScale', [0 1])
hsvxplot( PTF + exp(j*pi/2.01)*rot90(PTF) )
h = gca; h.FontSize = 80;

subplot(1,3,3)
% imagesc( gray2rgb(source, [1 1 1]) + gray2rgb(pupil_grad, [1 1 1]))
imagesc(source, [0 1])
% colormap gray
h = gca;
h.Colormap = gray;
axis image off

%%
source = zeros(M,N,3);
PTF = zeros(M,N,3);
for k = 1:3
source(:,:,k) = (1+cos(Theta+pi + (k-1)*2*pi/3))/2;
source(:,:,k) = circ( X,Y+0*(floor(N/4))*dx,N/2*dx) .* source(:,:,k) - ...
           circ( X,Y+0*(floor(N/4))*dx,(N/2-20)*dx) .* source(:,:,k);
[temp, ~] = getTransferFunctions(source(:,:,k), pupil);
PTF(:,:,k) = temp;
end

%%
n = 3;
figure(2)
subplot(1,2,1)
imagesc(source(:,:,n))
axis image off
colormap gray

subplot(1,2,2)
hsvxplot(PTF(:,:,n) .* exp(1j*(n-1)*2*pi/3))
%%

Hg = sum( exp(j*2*pi/3*reshape(0:2,[1,1,3]) ) .* PTF ,3);
figure(3)
hsvxplot(Hg,'colorbar','test')
axis image off
title('\Sigma')


% 
% %% simulate fully coherent point source a la FPM forward model
% % s = -sind(30); % 0.5
% % s = 0.5; % 0.5
% s = 0.45;
% % s = 0;
% % s = 0.51;
% illu = exp(1j * 2*pi/wavelength * (s*X)  );
% % object_imaged = abs(ifft2c(fft2c( object .* illu ) .* pupil)).^2;
% k = 0; % mirror padding
% object_imaged = abs(ifft2c(fft2c( mirror_pad(object .* illu, k) ) .* fft2c(mirror_pad(ifft2c( pupil ), k ) )) ).^2;
% 
% figure(3)
% imagesc((object_imaged), [0.5 1.5])
% % imagesc(sqrt(object_imaged), [0 0.4])
% axis image
% colormap gray
% xticks('')
% yticks('')
% zoom(1.8)
% 
% %%
% 
% % figure(3)
% % imagesc((object_imaged1 - object_imaged2)/2, [-1.5 1.5])
% % % imagesc(sqrt(object_imaged), [0 0.4])
% % axis image
% % colormap gray
% % xticks('')
% % yticks('')
% % zoom(1.8)
% 
% %%
% function r = crossGrating(object,X,Y,xs,ys,w,p,type)
% switch type
%     case 'phase'
%         support1 = rect( (X+xs) / (3*w) ) .* rect( (Y+ys) / w );
%         support2 = rect( (X+xs) / (3*w) ) .* rect( (Y+0.25*ys) / w );
%         r = support1 .* exp(1j*(cos(2*pi*X/p))) + ...
%             support2 .* exp(1j*(cos(2*pi*Y/p))) + ...
%             (1-support1-support2) .* object;
%         % r = support1 .* exp(1j*2*pi*X/p) + ...
%         %     support2 .* exp(1j*2*pi*Y/p) + ...
%         %     (1-support1-support2).*object;
%     case 'amplitude'
%         support1 = rect( (X+xs) / (3*w) ) .* rect( (Y+ys) / w );
%         support2 = rect( (X+xs) / (3*w) ) .* rect( (Y+0.25*ys) / w );
%         r = support1 .* ((cos(2*pi*X/p)+1)/2) + ...
%             support2 .* ((cos(2*pi*Y/p)+1)/2) + ...
%             (1 - support1 - support2).*object;
% end
% end
% 
% function r = starPattern(X,Y,afreq)
% 
% [theta, rho] = cart2pol(X,Y);
% 
% r = (1+cos(afreq*theta))/2;
% 
% end
% 
% 
% 
% 

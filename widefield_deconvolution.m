%% configuration
fsize = 32;
set(0, 'DefaultAxesFontSize', fsize);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');
addpath('utils')

% Parameters
lambda = 464e-9; % emission wavelength
NA = 0.8;        % Numerical aperture of the optical system
dx = 0.343e-6;     % lateral pixel size of the camera in meters
dz = 2e-6;       % axial pixel size

%% read data

file_name = '../datasets/widefield_fluorescence_DCV_1024.czi';
stack = CZIzstackReader(file_name);

pad_px = 20;
stack = mirror_pad(stack, pad_px);

N = size(stack, 1); % Size of the image (image_size x image_size pixels)
K = size(stack,3);

%% define 3D FFTs

FT3 = @(x) fftshift(fftn(ifftshift(x))) / sqrt(numel( x ) );
iFT3 = @(x) fftshift(ifftn(ifftshift(x))) * sqrt(numel( x ) );

%% get spatial frequencies

fx = (-N/2:N/2-1)/(N*dx);
fy = (-N/2:N/2-1)/(N*dx)';
fperp = sqrt(fx.^2 + fy'.^2);
fz = reshape((-(K-1)/2:(K-1)/2)/(K*dz), [1,1,K]);
[Fx, Fy] = meshgrid(fx,fy);

% compute pupil
pupil = circ(Fx,Fy,2*NA/lambda);

figure(1)
imagesc(pupil)
axis image off

%% compute 3D CSF, PSF, OTF

z = reshape( (-(K-1)/2:(K-1)/2)*dz, [1,1,K]);
csf = ifft2c( pupil .* exp(1i * 2*pi/lambda * z .* real(sqrt(1-(lambda*fperp).^2))) );
psf = abs(csf).^2;
OTF = FT3(psf);

figure(2)
imagesc(psf(:,:,21))
axis square
colormap gray
zoom(4)

figure(3)
subplot(1,2,1)
imagesc(abs(OTF(:,:,21)))
axis image off

subplot(1,2,2)
imagesc(abs(squeeze(OTF(:,end/2+1,:))))
axis square off

%% Wiener filter

tic
regularization = 1e5;
wiener_dcv = abs( iFT3( FT3(stack) .* conj(OTF) ./ (abs(OTF).^2 + regularization) ) );
toc

%% crop reconstructed data

wiener_dcv = crop(wiener_dcv, pad_px);

%% visualize Wiener filter result

plane = 21;

figure(4)
subplot(1,2,1)
imagesc(stack(:,:,plane)', [0 0.75*max(stack(:,:,plane),[],[1,2])])
axis image off
colormap gray
title('Wide field')

subplot(1,2,2)
imagesc(wiener_dcv(:,:,plane)', [0 0.75*max(wiener_dcv(:,:,plane),[],[1,2])])
axis image off
colormap gray
title('Wiener')


%% Richardsen Lucy DCV

% initialization
rl_dcv = stack;
stack = stack / max(stack(:));
num_iter = 200;
regularization = 1e-3;
denominator = iFT3(FT3(ones(size(stack))) .* conj(OTF));
residual = [];
figureUpdate = 25;

% switches
useGPU = true;

if useGPU
    denominator = gpuArray(single(denominator));
    OTF = gpuArray(single(OTF));
    stack = gpuArray(single(stack));
    rl_dcv= gpuArray(single(rl_dcv));
end

tic
for loop = 1:num_iter
    
    % RL DCV
    y = iFT3(FT3(rl_dcv) .* OTF);
    numerator = real( iFT3(FT3((stack + aleph)./(y + aleph)) .* conj(OTF)) );
    rl_dcv = (numerator./denominator) .* rl_dcv;

    % compute residual
    residual = [residual, real(gather(-sum( stack(:).*log(y(:)) - y(:) )))];

    if mod(pad_px,figureUpdate)==0
        
        figure(5)
        subplot(1,2,1)
        imagesc(gather(stack(:,:,plane))', [0 0.75*max(stack(:,:,plane),[],[1,2])])
        axis image off
        colormap gray
        title('wide field')
        % colorbar

        subplot(1,2,2)
        imagesc(gather(rl_dcv(:,:,plane))', [0 0.3*max(rl_dcv(:,:,plane),[],[1,2])])
        axis image off
        colormap gray
        % colorbar
        title(['RL iterations: ', num2str(loop)])

        if loop >= 1
            figure(6)
            loglog(1:loop,residual(1:loop)/residual(1),'-ok','linewidth',2,'MarkerFaceColor','k')
            xlabel('iteration'), ylabel('residual')
            axis square
            grid on
        end
        drawnow
    end
end
toc

%% crop reconstructed data

rl_dcv = crop(rl_dcv,pad_px);

%% max projection

figure(100)
subplot(1,3,1)
imagesc(max(stack,[],3)', [0 1*max(max(stack,[],3), [],[1,2])])
axis image off
colormap gray

subplot(1,3,2)
imagesc(max(wiener_dcv,[],3)', [0 0.7*max(max(wiener_dcv,[],3), [],[1,2])])
axis image off
colormap gray

subplot(1,3,3)
imagesc(max(rl_dcv,[],3)', [0 0.45*max(max(rl_dcv,[],3), [],[1,2])])
axis image off
colormap gray

%% axial cross-section

figure(102)
subplot(1,3,1)
imagesc(squeeze(stack(:,end/2,:))', [0 0.7*max(max(stack(:,:,plane),[],3), [],[1,2])])
axis square off
colormap gray

subplot(1,3,2)
imagesc(squeeze(wiener_dcv(:,end/2,:))', [0 0.45*max(max(wiener_dcv(:,:,plane),[],3), [],[1,2])])
axis square off
colormap gray

subplot(1,3,3)
imagesc(squeeze(rl_dcv(:,end/2,:))', [0 0.2*max(max(rl_dcv(:,:,plane),[],3), [],[1,2])])
axis square off
colormap gray
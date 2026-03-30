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
colormap gray

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


%% Richardson-Lucy DCV with Anderson Acceleration

% initialization
rl_dcv = stack;
stack = stack / max(stack(:));
num_iter = 200;
regularization = 1e-3;
denominator = iFT3(FT3(ones(size(stack))) .* conj(OTF));
residual = [];
figureUpdate = 25;

% Anderson Acceleration parameters
aa_depth = 5;              % number of previous iterates to store (m)
aa_start = 1;              % iteration to start AA
aa_beta = 1e-8;            % regularization for least squares
use_AA = true;             % enable/disable AA

% AA storage
aa_F = [];                 % residual history: F_k = g(x_k) - x_k
aa_X = [];                 % iterate history
aa_iter = 0;

% switches
useGPU = true;

if useGPU
    denominator = gpuArray(single(denominator));
    OTF = gpuArray(single(OTF));
    stack = gpuArray(single(stack));
    rl_dcv = gpuArray(single(rl_dcv));
end

tic
for loop = 1:num_iter
    
    % Standard RL update (fixed-point map g(x))
    y = iFT3(FT3(rl_dcv) .* OTF);
    numerator = real(iFT3(FT3((stack + regularization)./(y + regularization)) .* conj(OTF)));
    rl_dcv_new = (numerator./denominator) .* rl_dcv;
    
    % Anderson Acceleration
    if use_AA && loop >= aa_start
        % Compute residual: F_k = g(x_k) - x_k
        F_k = rl_dcv_new(:) - rl_dcv(:);
        
        % Store current iterate and residual
        if aa_iter == 0
            aa_F = F_k;
            aa_X = rl_dcv(:);
        else
            aa_F = [aa_F, F_k];
            aa_X = [aa_X, rl_dcv(:)];
        end
        
        % Limit depth
        m_k = min(aa_iter, aa_depth);
        if size(aa_F, 2) > aa_depth
            aa_F = aa_F(:, end-aa_depth+1:end);
            aa_X = aa_X(:, end-aa_depth+1:end);
        end
        
        % Compute differences
        if m_k > 0
            dF = diff(aa_F, 1, 2);  % F_{k} - F_{k-1}, ..., F_{k} - F_{k-m_k}
            
            % Solve least squares: min ||dF * gamma - F_k||^2 + beta*||gamma||^2
            % Using normal equations with Tikhonov regularization
            M = dF' * dF + aa_beta * eye(m_k);
            rhs = dF' * F_k;
            gamma = M \ rhs;
            
            % Compute accelerated iterate
            dX = diff(aa_X, 1, 2);
            x_aa = rl_dcv_new(:) - dF * gamma + dX * gamma;
            
            % Reshape and apply non-negativity constraint
            rl_dcv = reshape(x_aa, size(rl_dcv));
            rl_dcv = max(rl_dcv, 0);  % ensure non-negativity
        else
            rl_dcv = rl_dcv_new;
        end
        
        aa_iter = aa_iter + 1;
    else
        rl_dcv = rl_dcv_new;
    end
    
    % Compute residual (Poisson log-likelihood)
    y = iFT3(FT3(rl_dcv) .* OTF);  % recompute with updated rl_dcv
    residual = [residual, real(gather(-sum(stack(:).*log(y(:)+eps) - y(:))))];
    
    % Visualization
    if mod(loop, figureUpdate) == 0
        
        figure(5)
        subplot(1,2,1)
        imagesc(gather(stack(:,:,plane))', [0 0.75*max(stack(:,:,plane),[],[1,2])])
        axis image off
        colormap gray
        title('Wide Field')
        
        subplot(1,2,2)
        imagesc(gather(rl_dcv(:,:,plane))', [0 0.3*max(rl_dcv(:,:,plane),[],[1,2])])
        axis image off
        colormap gray
        if use_AA
            title(['RL+AA (m=', num2str(aa_depth), '), iter: ', num2str(loop)])
        else
            title(['RL iterations: ', num2str(loop)])
        end
        
        if loop >= 1
            figure(6)
            loglog(1:loop, residual(1:loop)/residual(1), '-ok', 'linewidth', 2, 'MarkerFaceColor', 'k')
            xlabel('Iteration'), ylabel('Normalized Residual')
            axis square
            grid on
            if use_AA
                title(['Convergence (AA depth m=', num2str(aa_depth), ')'])
            end
        end
        drawnow
    end
end
toc

fprintf('Final residual reduction: %.2e\n', residual(end)/residual(1));

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
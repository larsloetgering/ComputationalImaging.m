function [g_optimal, p_optimal, gcv_values, p_values] = gcv_deconvolution(I, h, varargin)
% GCV_DECONVOLUTION Deconvolution using Generalized Cross-Validation
%
% Inputs:
%   I       - Observed blurred image (2D array)
%   h       - Point spread function/convolution kernel (2D array)
%   varargin - Optional parameters as name-value pairs:
%              'p_min'    - Minimum regularization parameter (default: 1e-6)
%              'p_max'    - Maximum regularization parameter (default: 1e2)
%              'n_p'      - Number of p values to test (default: 50)
%              'reg_type' - Regularization type: 'tikhonov' or 'gradient' (default: 'tikhonov')
%              'verbose'  - Display progress (default: true)
%
% Outputs:
%   g_optimal  - Reconstructed image at optimal p
%   p_optimal  - Optimal regularization parameter
%   gcv_values - GCV function values for all tested p
%   p_values   - Tested regularization parameters

%% Parse input arguments
p = inputParser;
addRequired(p, 'I');
addRequired(p, 'h');
addParameter(p, 'p_min', 1e-6);
addParameter(p, 'p_max', 1e2);
addParameter(p, 'n_p', 50);
addParameter(p, 'reg_type', 'tikhonov');
addParameter(p, 'verbose', true);
parse(p, I, h, varargin{:});

p_min = p.Results.p_min;
p_max = p.Results.p_max;
n_p = p.Results.n_p;
reg_type = p.Results.reg_type;
verbose = p.Results.verbose;

%% Setup
[M, N] = size(I);
n_pixels = M * N;

% Ensure h has same size as I (zero-pad if necessary)
[Mh, Nh] = size(h);
h_padded = zeros(M, N);
h_padded(1:Mh, 1:Nh) = h;

% Shift kernel to center (important for proper FFT convolution)
h_padded = ifftshift(h_padded);

% Compute FFT of kernel and image
H_fft = fft2(h_padded);
I_fft = fft2(I);

%% Setup regularization operator in Fourier domain
switch lower(reg_type)
    case 'tikhonov'
        % Identity operator: L = I, so L^T L = I
        % In Fourier domain: |F{L}|^2 = 1
        L_fft_squared = ones(M, N);
        
    case 'gradient'
        % Gradient operator (Laplacian): approximates smoothness
        % In Fourier domain, use discrete Laplacian
        [kx, ky] = meshgrid(0:N-1, 0:M-1);
        kx = kx - floor(N/2);
        ky = ky - floor(M/2);
        
        % Laplacian in Fourier domain: |k|^2
        L_fft_squared = (2*pi*kx/N).^2 + (2*pi*ky/M).^2;
        L_fft_squared = ifftshift(L_fft_squared);
        
    otherwise
        error('Unknown regularization type: %s', reg_type);
end

%% Generate candidate p values (logarithmically spaced)
p_values = logspace(log10(p_min), log10(p_max), n_p);
gcv_values = zeros(n_p, 1);

%% Compute GCV for each p value
if verbose
    fprintf('Computing GCV for %d regularization parameters...\n', n_p);
    tic;
end

for i = 1:n_p
    p_current = p_values(i);
    
    % Solve in Fourier domain: G(p) = (H^T H + p L^T L)^{-1} H^T I
    % In Fourier space: G_fft = (|H|^2 + p |L|^2)^{-1} * conj(H) * I_fft
    denominator = abs(H_fft).^2 + p_current * L_fft_squared;
    
    % Avoid division by zero
    denominator(denominator < eps) = eps;
    
    % Compute solution in Fourier domain
    G_fft = (conj(H_fft) .* I_fft) ./ denominator;
    
    % Transform back to spatial domain
    G = real(ifft2(G_fft));
    
    % Compute residual: r = I - H * G
    HG_fft = H_fft .* G_fft;
    HG = real(ifft2(HG_fft));
    residual = I - HG;
    
    % Numerator: ||I - H*G||^2
    numerator = sum(residual(:).^2);
    
    % Denominator: [tr(I - A(p))]^2
    % tr(A(p)) = sum of diagonal elements of influence matrix
    % In Fourier domain: tr(A(p)) = sum of |H|^2 / (|H|^2 + p |L|^2)
    influence_trace = sum(abs(H_fft(:)).^2 ./ denominator(:));
    
    % Effective degrees of freedom
    dof = n_pixels - influence_trace;
    
    % GCV formula
    if abs(dof) < eps
        gcv_values(i) = inf;
    else
        gcv_values(i) = (n_pixels * numerator) / (dof^2);
    end
    
    if verbose && mod(i, max(1, floor(n_p/10))) == 0
        fprintf('  Progress: %d/%d (%.1f%%)\n', i, n_p, 100*i/n_p);
    end
end

if verbose
    fprintf('Completed in %.2f seconds.\n', toc);
end

%% Find optimal p
[~, idx_optimal] = min(gcv_values);
p_optimal = p_values(idx_optimal);

if verbose
    fprintf('Optimal regularization parameter: p = %.6e\n', p_optimal);
end

%% Compute optimal solution
denominator_opt = abs(H_fft).^2 + p_optimal * L_fft_squared;
denominator_opt(denominator_opt < eps) = eps;
G_optimal_fft = (conj(H_fft) .* I_fft) ./ denominator_opt;
g_optimal = real(ifft2(G_optimal_fft));

%% Display results
if verbose
    figure('Name', 'GCV Analysis Results');
    
    % Plot 1: GCV curve
    subplot(2,3,1);
    semilogx(p_values, gcv_values, 'b-', 'LineWidth', 1.5);
    hold on;
    semilogx(p_optimal, gcv_values(idx_optimal), 'ro', ...
             'MarkerSize', 10, 'MarkerFaceColor', 'r');
    xlabel('Regularization Parameter $$p$$', 'Interpreter', 'latex');
    ylabel('GCV Function', 'Interpreter', 'latex');
    title('GCV Curve', 'Interpreter', 'latex');
    grid on;
    legend('GCV', 'Optimal $$p$$', 'Interpreter', 'latex', 'Location', 'best');
    
    % Plot 2: Observed image
    subplot(2,3,2);
    imagesc(I);
    colormap(gray);
    colorbar;
    axis image;
    title('Observed Image $$I$$', 'Interpreter', 'latex');
    
    % Plot 3: PSF/Kernel
    subplot(2,3,3);
    imagesc(h);
    colormap(gray);
    colorbar;
    axis image;
    title('Convolution Kernel $$h$$', 'Interpreter', 'latex');
    
    % Plot 4: Reconstructed image
    subplot(2,3,4);
    imagesc(g_optimal);
    colormap(gray);
    colorbar;
    axis image;
    title(sprintf('Reconstructed $$G$$ ($$p$$ = %.2e)', p_optimal), ...
          'Interpreter', 'latex');
    
    % Plot 5: Residual
    residual_opt = I - real(ifft2(H_fft .* G_optimal_fft));
    subplot(2,3,5);
    imagesc(residual_opt);
    colormap(gray);
    colorbar;
    axis image;
    title('Residual $$I - h \star G$$', 'Interpreter', 'latex');
    
    % Plot 6: Residual statistics
    subplot(2,3,6);
    histogram(residual_opt(:), 50);
    xlabel('Residual Value', 'Interpreter', 'latex');
    ylabel('Count', 'Interpreter', 'latex');
    title('Residual Distribution', 'Interpreter', 'latex');
    grid on;
end

end

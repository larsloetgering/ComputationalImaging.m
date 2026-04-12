function phi = tie(I, dz, wavelength, pixel_size, reg_param, pad_px)
% Solves the Transport of Intensity Equation, neglecting intensity gradients
% in illumination
% 
% EXAMPLE use:
% phi = tie(I, dz, wavelength, pixel_size, regLaplace, pad_px)
%
% INPUTS:
%   I           - focus stack of three adjacent, equispaced intensities 
%   dz          - Defocus distance (in meters), equidistant from focus
%   wavelength  - Wavelength of illumination (in meters)
%   pixel_size  - Pixel size in object plane (in meters)
%   pad_px      - number of pixels to be padded    
%
% OUTPUT:
%   phi         - Recovered phase (in radians)
    
    % 2D intensity image before focal plane (defocused by -dz)
    I_minus = mirror_pad(I(:,:,1),pad_px); 
    % 2D intensity image at focal plane
    I_focus = mirror_pad(I(:,:,2),pad_px); 
    % 2D intensity image after focal plane (defocused by +dz)
    I_plus = mirror_pad(I(:,:,3),pad_px);
    
    %% Get image dimensions
    [Ny, Nx] = size(I_focus);
    
    %% Create frequency coordinates
    % Frequency spacing
    dfx = 1 / (Nx * pixel_size);
    dfy = 1 / (Ny * pixel_size);
    
    % Frequency vectors (centered)
    fx = ifftshift((-floor(Nx/2) : ceil(Nx/2)-1) * dfx);
    fy = ifftshift((-floor(Ny/2) : ceil(Ny/2)-1) * dfy);
    
    % 2D frequency grids
    [Fx, Fy] = meshgrid(fx, fy);
    
    % Squared transverse frequency
    F_perp_sq = Fx.^2 + Fy.^2;
    
    %% Compute inverse Laplacian operator with regularization

    inv_laplacian = F_perp_sq ./( F_perp_sq.^2 + (reg_param/pixel_size).^2 );
    
    %% Step 1: Compute longitudinal intensity derivative (dI/dz)

    % Using central difference approximation
    dI_dz = pi / (wavelength*dz) * (I_plus - I_minus) / mean(I_focus(:));
    % note: by assumption, the object is transparent and we have a
    % homogeneous illumination; hence mean(I) in denominator,
    % which avoids risk of hot/dead pixels

    %% Step 2: Apply inverse Laplacian to dI/dz
    dI_dz_fft = fft2(dI_dz);
    phi_fft = dI_dz_fft .* inv_laplacian;
    phi = crop(real(ifft2(phi_fft)),pad_px);
    % zero mean constraint (absolute phase not uniquely determined)
    phi = phi-mean(phi(:)); 
end

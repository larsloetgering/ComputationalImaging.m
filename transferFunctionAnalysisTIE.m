clear
fsize = 32;
set(0, 'DefaultAxesFontSize', fsize);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');

%% parameters

pixel_size = 172e-9;
wavelength = 500e-9;
N = 1024;
NAi = 0.27;
NAo = 0.8;
Delta = wavelength / (1 - sqrt(1-NAo^2));
dz = Delta ;          % axial slice separation

% get Fourier space coordinates

dfx = 1 / (N * pixel_size);
dfy = 1 / (N * pixel_size);
fx = (-N/2 : N/2-1) * dfx;
fy = (-N/2 : N/2-1)' * dfy;
[Fx,Fy] = meshgrid(fx,fy);
fmax_coh = NAo/wavelength;

% simulate point source illumination

source = circ(Fx,Fy,2*NAi/wavelength);
pupil  = circ(Fx,Fy,2*NAo/wavelength) .* ...
        exp(1i*2*pi/wavelength*dz*real(sqrt(1-(wavelength*Fx).^2-(wavelength*Fy).^2)));

%

[PTF, ATF] = getTransferFunctions(source, pupil);

figure(1)
imagesc(real(PTF), [-1 1])
axis image off
colormap(greenWhiteBlue(256))

figure(2)
plot(fx/fmax_coh,real(PTF(end/2,:)),'k-','lineWidth',2)
axis([min(fx)/fmax_coh, max(fx)/fmax_coh, -1 1])
xlabel('$f_{\perp} [ NA / \lambda ]$')
axis square
grid on
h = gca;
h.GridLineWidth = 2;

figure(3)
imagesc(source)
axis image
colormap(flipud(gray))


function cmap = greenWhiteBlue(n)
% greenWhiteBlue creates a divergent colormap: green -> white -> magenta
%
%   cmap = greenWhiteBlue(n) returns an n-by-3 colormap matrix.
%   Default n = 256 if not specified.

    if nargin < 1
        n = 256;
    end
    
    % Define anchor colors (RGB values normalized to [0, 1])
    green   = [0, 0.7, 0];      % Green
    white   = [1, 1, 1];        % White
    blue = [0, 0, 0.7];    % Magenta
    
    % Number of points for each half
    half = floor(n / 2);
    
    % Interpolate green -> white
    r1 = linspace(green(1), white(1), half)';
    g1 = linspace(green(2), white(2), half)';
    b1 = linspace(green(3), white(3), half)';
    
    % Interpolate white -> magenta
    r2 = linspace(white(1), blue(1), n - half)';
    g2 = linspace(white(2), blue(2), n - half)';
    b2 = linspace(white(3), blue(3), n - half)';
    
    % Combine both halves
    cmap = [r1, g1, b1; r2, g2, b2];
end

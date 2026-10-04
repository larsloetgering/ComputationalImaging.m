% config
fsize = 23;
set(0, 'DefaultAxesFontSize', fsize);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');
addpath('utils')
%% Ewald sphere (2D) visualization for arbitrary illumination + finite-NA detection
% Didactic tool: draws the spherical caps (arcs in k-space) that a specimen
% samples when illuminated from an arbitrary direction and detected through an
% objective of a given numerical aperture. Only (kx, kz) are considered; ky is
% ignored, so the "Ewald sphere" is a 2D circle and each measurement is an arc.
%
% Physics (Fourier diffraction theorem, elastic scattering):
%   incident wavevector : k_in  = k * [sin(alpha), cos(alpha)]      (|k_in| = k)
%   detected wavevector : k_out = k * [sin(theta), cos(theta)],  theta in +/- asin(NA)
%   measured object freq: K     = k_out - k_in
% For a fixed illumination this traces a circular arc of radius k centered at
% -k_in, spanning the objective acceptance cone about the +z detection axis.
%
% Special cases obtainable with the parameters below:
%   a) OCT   : illumAngles_deg = -180, sweep wavelengths_nm  -> caps stack near (0, 2k)
%   b) ODT   : illumAngles_deg in +/- asind(NA), single wavelength -> rotated caps near origin
%   c) Dark field : |illumAngles_deg| > asind(NA) (illumination NA > objective NA)

clear; clc; close all;

%% ---------------- User parameters ----------------

% illumAngles_deg = -45:1:45;transparency = 0.15;      % (1) illumination directions [deg], row vector in [-180, 180] (tODT)
% illumAngles_deg = [(-90:-45),(45:90)]; transparency = 0.15;      % (1) illumination directions [deg], row vector in [-180, 180] (t dark)
% illumAngles_deg = [(-180:1:-135),(135:1:179)]; transparency = 0.05;      % (1) illumination directions [deg], row vector in [-180, 180] (r dark)
% illumAngles_deg = [(-135:1:-90),(90:1:135)]; transparency = 0.2;      % (1) illumination directions [deg], row vector in [-180, 180] (r dark)
% illumAngles_deg = [(-90:1:90)]; transparency = 0.15;      % combine eEPI, td, bf
% illumAngles_deg = [(-5:1:5)]; transparency = 0.05;      % combine eEPI, td, bf
% illumination_deg = -
% dtheta=10;illumAngles_deg = -(180-dtheta):dtheta:180; transparency = 0.5;       % (1) illumination directions [deg], row vector in [-180, 180] (the diffraction "torus")
dtheta=10;illumAngles_deg = -(130):dtheta:130; transparency = 1;       % (1) illumination directions [deg], row vector in [-180, 180] (the diffraction "torus")

% illumAngles_deg = [(-129:1:-90),(90:1:129)];      % (1) illumination directions [deg], row vector in [-180, 180] (r dark)

% illumAngles_deg = 0;      % (1) illumination directions [deg], row vector in [-180, 180] (r dark)
NA              = 0.7;              % (2) numerical aperture of the DETECTION objective
% wavelengths_nm  = 450:100:650;   % (3) illumination wavelengths [nm], row vector
% wavelengths_nm  = 450:5:650;         % OCT
wavelengths_nm = 450;
% transparency    = 0.02;             % curve opacity: 0 = fully transparent, 1 = opaque

% --- Example presets (uncomment one) ---
% a) OCT:        illumAngles_deg = -180;                 wavelengths_nm = 700:20:900;
% b) ODT:        illumAngles_deg = linspace(-asind(NA), asind(NA), 15);  wavelengths_nm = 532;
% c) Dark field: illumAngles_deg = [-120 -90 90 120];    wavelengths_nm = 532;   % |angle| > asind(NA)

%% ---------------- Precompute ----------------
theta_NA = asin(NA);                     % objective acceptance half-angle [rad]
thetaCap = linspace(-theta_NA, theta_NA, 200);   % sampling along each cap

figure('Color', 'w'); hold on; axis equal; grid on; box on;

legHandles = gobjects(1, numel(wavelengths_nm));
legLabels  = cell(1, numel(wavelengths_nm));

%% ---------------- Nested loop: illumination angles x wavelengths ----------------
for iw = 1:numel(wavelengths_nm)
    for ia = 1:numel(illumAngles_deg)
        alpha = deg2rad(illumAngles_deg(ia));

    
        lam_nm = wavelengths_nm(iw);
        k = 2*pi / (lam_nm * 1e-3);      % wavenumber [rad/um]  (lambda converted nm -> um)

        % incident wavevector and Ewald-cap center (-k_in)
        k_in = k * [sin(alpha), cos(alpha)];
        c    = -k_in;

        % arc of measured spatial frequencies K = k_out - k_in
        Kx = c(1) + k * sin(thetaCap);
        Kz = c(2) + k * cos(thetaCap);

        % color = wavelength (numColors = 1 -> pure spectral RGB)
        rgb = wavelength2rgb(lam_nm);
        rgb = min(max(double(rgb), 0), 1);      % clamp to valid [0,1]

        p = plot(Kx, Kz, 'LineWidth', 2);
        % p.Color = [rgb, transparency];          % 4th element = per-line transparency
        p.Color = [[0 0 0], transparency];
        % keep one handle per wavelength for the legend
        if ~isgraphics(legHandles(iw))
            legHandles(iw) = p;
            legLabels{iw}  = sprintf('$%g$ nm', lam_nm);
        end
    end
end

%% ---------------- Limiting circle (backscatter bound, bluest wavelength) ----------------
% Same radius/center as the -180 deg cap at the shortest (bluest) wavelength,
% i.e. same curvature and tangent to that cap, just completed to a full circle.
grid on
lambda_blue_nm = min(wavelengths_nm);
k_blue         = 2*2*pi / (lambda_blue_nm * 1e-3);
alpha_ref      = deg2rad(-180);
c_ref          = -k_blue * [sin(alpha_ref), cos(alpha_ref)];
tFull          = linspace(0, 2*pi, 361);
Kx_lim = c_ref(1) + k_blue * sin(tFull);
Kz_lim = 0.01*c_ref(2) + k_blue * cos(tFull);
axis(k_blue*[-1 1 -0.15 1])
xticklabels([])
yticklabels([])
h = gca;
h.GridLineWidth  = 2;
% xticks(''), yticks('')
% hLim = plot(Kx_lim, Kz_lim, 'k--', 'LineWidth', 2);

% ---------------- Cosmetics ----------------
plot(0, 0, 'k+', 'MarkerSize', 20, 'LineWidth', 2);   % origin (DC / unscattered)
% xlabel('$k_x$ [rad/$\mu$m]');
% ylabel('$k_z$ [rad/$\mu$m]');
% axis off
% title(sprintf('2D Ewald sphere coverage ($NA = %.2f$, $\\theta_{NA} = \\pm%.1f^{\\circ}$)', ...
%       NA, rad2deg(theta_NA)));
% legend([legHandles, hLim], [legLabels, {sprintf('Limiting circle ($%g$ nm, $-180^{\\circ}$)', lambda_blue_nm)}], ...
%        'Location', 'bestoutside');
% set(gca, 'FontSize', 11);

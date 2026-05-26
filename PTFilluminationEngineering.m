fsize = 23;
set(0, 'DefaultAxesFontSize', fsize);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');
addpath(genpath('utils'))
% compute transfer functions

N = 1000;
x = (-(N-mod(N,2))/2:(N+mod(N,2))/2-1);
[X,Y] = meshgrid(x);

source = zeros(N,N,3);
source(:,:,1) = circ(X,Y,N/4) .* (X>0);   % unmatched NA
source(:,:,2) = circ(X,Y,N/2-1) .* (X>0); % matched NA
source(:,:,3) = (circ(X,Y,N/2-1) - circ(X,Y,N/2-15)) .* (X>0); % annular matched

pupil = circ(X,Y,N/2-1);

show_index = 3;
figure(1)
imagesc(source(:,:,show_index))
axis image off
colormap gray

figure(2), clf
imagesc(pupil)
axis image off
colormap gray

%% compute ATF and PTF
PTF = 0*source;
for k = 1:3
    [PTF_temp, ~] = getTransferFunctions(imfilter(source(:,:,k), fspecial('Gaussian', 8, 1),'same'), ...
        imfilter(pupil, fspecial('Gaussian', 8, 1),'same'));
    PTF(:,:,k) = PTF_temp;
end

show_index = 3;
figure(3), clf
imagesc(imag(PTF(:,:,show_index)))
axis image off
% colormap(slanCM('bwr'))
colormap(slanCM('RdGy'))
% colormap(diverging_cmap(100, 'white'))

%% show cross sections
markerTypes = {'k-','k--','k-'}

figure(4), clf
hold on
for k = 1:3
    plot(1:N,imag(PTF(end/2+1,:,k)),markerTypes{k},'lineWidth', 2)
    if k == 3
        plot(1:40:N,imag(PTF(end/2+1,1:40:end,k)),'ko','lineWidth', 2,'MarkerFaceColor','k')
    end
    axis square off
end

function cmap = diverging_cmap(n, variant)
%DIVERGING_CMAP Diverging colormap with optional protan-oriented variant.
%
%   CMAP = DIVERGING_CMAP(N) returns an N-by-3 colormap that transitions
%   from a cool color, through a neutral point, to a warm color.
%
%   CMAP = DIVERGING_CMAP(N, VARIANT) lets you choose the endpoint hues:
%
%       VARIANT = 'default'  : dark blue -> gray -> light orange
%       VARIANT = 'protan'   : dark blue -> white -> vivid light yellow
%                             (designed to be more robust for protan
%                             color-vision deficiencies; center is white)
%
%   Inputs
%   ------
%   N       : Number of colors (default 256, minimum 2).
%   VARIANT : (optional) 'default' or 'protan'.
%
%   Example
%   -------
%       imagesc(peaks(200));
%       axis image
%       colormap(diverging_cmap(256, 'protan'));  % protan-oriented, white center
%       colorbar

    if nargin < 1 || isempty(n)
        n = 256;
    end
    if nargin < 2 || isempty(variant)
        variant = 'default';
    end

    % Ensure integer >= 2
    n = max(2, round(n));

    % Normalize and interpret variant
    variant = lower(string(variant));
    switch variant
        case {"default", "blue-orange"}
            % ---- DEFAULT: dark blue -> gray -> light orange ----
            %
            % Center (neutral gray)
            L0 = 70;    % L* of central gray
            C0 = 0;     % C* = 0 at the center (neutral)
            % Left (blue)
            L1 = 40;    % darker than center
            C1 = 55;    % chroma
            h1 = 260;   % hue in degrees (blue)
            % Right (orange)
            L2 = 90;    % lighter than center
            C2 = 55;    % chroma
            h2 = 40;    % hue in degrees (orange)

        case {"protan", "protan-safe", "blue-yellow"}
            % ---- PROTAN-ORIENTED: dark blue -> white -> vivid light yellow ----
            %
            % Center (white)
            % CIELAB(D65) white: L* = 100, a* = 0, b* = 0  =>  C* = 0
            L0 = 100;   % central white
            C0 = 0;     % neutral (no chroma)

            % Left endpoint (dark, saturated blue)
            L1 = 20;    % darker -> more contrast vs white
            C1 = 60;    % stronger chroma for a more vivid blue
            h1 = 250;   % slightly purple-blue

            % Right endpoint (light, more saturated yellow)
            L2 = 90;    % light yellow, noticeably darker than white
            C2 = 55;    % stronger chroma for more saturated yellow
            h2 = 90;    % yellow (stays on blue–yellow axis, not orange)

        otherwise
            error('Unknown VARIANT "%s". Use ''default'' or ''protan''.', variant);
    end

    % ---- 1. Define diverging path in CIE L*C*h° ----
    %
    % Parameter t runs from -1 (left) to +1 (right)
    t = linspace(-1, 1, n).';

    L = zeros(n, 1);
    C = zeros(n, 1);
    h = zeros(n, 1);

    % Left side: t < 0 (cool -> neutral)
    maskL = t < 0;
    sL = -t(maskL);  % maps t in [-1,0) to s in (0,1]
    L(maskL) = L0 + (L1 - L0) .* sL;
    C(maskL) = C0 + (C1 - C0) .* sL;
    h(maskL) = h1;

    % Right side: t > 0 (neutral -> warm)
    maskR = t > 0;
    sR = t(maskR);   % maps t in (0,1] to s in (0,1]
    L(maskR) = L0 + (L2 - L0) .* sR;
    C(maskR) = C0 + (C2 - C0) .* sR;
    h(maskR) = h2;

    % Center: t = 0 (exactly neutral, if present)
    mask0 = (t == 0);
    L(mask0) = L0;
    C(mask0) = C0;
    % h is irrelevant at the center because C* = 0

    % ---- 2. Convert L*C*h° -> L*a*b* ----
    a = C .* cosd(h);
    b = C .* sind(h);
    Lab = [L, a, b];

    % ---- 3. Convert L*a*b* -> sRGB ----
    if exist('lab2rgb', 'file') == 2
        % Use Image Processing Toolbox if available
        cmap = lab2rgb(Lab, 'OutputType', 'double');
    else
        % Fallback: built-in implementation
        cmap = lab2rgb_local(Lab);
    end

    % Clip to [0,1] for numerical safety
    cmap = min(max(cmap, 0), 1);
end

% -------------------------------------------------------------------------
% Local helper: CIELAB (D65) -> sRGB conversion
% -------------------------------------------------------------------------
function rgb = lab2rgb_local(Lab)
%LAB2RGB_LOCAL Approximate CIELAB(D65) to sRGB conversion.
%
%   RGB = LAB2RGB_LOCAL(LAB) converts an N-by-3 matrix of CIELAB values
%   (rows are [L*, a*, b*] with L* in [0,100]) to an N-by-3
%   matrix of sRGB values in [0,1].

    L = Lab(:,1);
    a = Lab(:,2);
    b = Lab(:,3);

    % Reference white (D65)
    Xn = 95.047;
    Yn = 100.000;
    Zn = 108.883;

    % Lab -> XYZ
    fy = (L + 16) / 116;
    fx = fy + a / 500;
    fz = fy - b / 200;

    xyz = [fx, fy, fz];

    delta = 6/29;
    mask = xyz > delta;
    xyz(mask)  = xyz(mask)  .^ 3;
    xyz(~mask) = (xyz(~mask) - 16/116) / 7.787;

    X = Xn * xyz(:,1);
    Y = Yn * xyz(:,2);
    Z = Zn * xyz(:,3);

    % XYZ -> linear sRGB
    % Matrix for D65-adapted sRGB (IEC 61966-2-1)
    M = [ 3.2404542, -1.5371385, -0.4985314; ...
         -0.9692660,  1.8760108,  0.0415560; ...
          0.0556434, -0.2040259,  1.0572252];

    XYZ = [X, Y, Z] / 100;          % scale to reference white
    rgb_lin = XYZ * M.';            % linear sRGB

    % Clamp negatives before gamma correction
    rgb_lin = max(rgb_lin, 0);

    % Linear sRGB -> gamma-encoded sRGB
    threshold = 0.0031308;
    rgb = zeros(size(rgb_lin));
    mask = rgb_lin <= threshold;
    rgb(mask)  = 12.92 * rgb_lin(mask);
    rgb(~mask) = 1.055 * (rgb_lin(~mask).^(1/2.4)) - 0.055;

    % Final clamp to [0,1]
    rgb = min(max(rgb, 0), 1);
end

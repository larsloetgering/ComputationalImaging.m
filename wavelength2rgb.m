function rgb = wavelength2rgb(wavelength_nm)
%WAVELENGTH2RGB Approximate RGB color for a visible wavelength.
%   rgb = wavelength2rgb(630) returns the RGB triplet for 630 nm light.
%   Input wavelength is in nanometers and must be in the range 350..800.
%   Output is an RGB triplet in the range [0, 1].
%
%   This implementation does not depend on wave2rgb.m.

    if nargin ~= 1
        error('wavelength2rgb:InvalidInput', 'wavelength2rgb requires exactly one input value.');
    end

    if ~isnumeric(wavelength_nm) || isempty(wavelength_nm) || any(~isfinite(wavelength_nm(:)))
        error('wavelength2rgb:InvalidInput', 'Input wavelength must be a finite numeric scalar or vector.');
    end

    wl = double(wavelength_nm(:));

    if any(wl < 350) || any(wl > 800)
        error('wavelength2rgb:OutOfRange', 'Wavelength must be between 350 and 800 nm.');
    end

    n = numel(wl);
    rgb = zeros(n, 3);

    for k = 1:n
        x = wl(k);

        if x < 380
            r = 0;
            g = 0;
            b = 1;
            factor = (x - 350) / (380 - 350);
        elseif x >= 380 && x < 440
            r = -(x - 440) / (440 - 380);
            g = 0;
            b = 1;
            factor = 1;
        elseif x >= 440 && x < 490
            r = 0;
            g = (x - 440) / (490 - 440);
            b = 1;
            factor = 1;
        elseif x >= 490 && x < 510
            r = 0;
            g = 1;
            b = -(x - 510) / (510 - 490);
            factor = 1;
        elseif x >= 510 && x < 580
            r = (x - 510) / (580 - 510);
            g = 1;
            b = 0;
            factor = 1;
        elseif x >= 580 && x < 645
            r = 1;
            g = -(x - 645) / (645 - 580);
            b = 0;
            factor = 1;
        elseif x >= 645 && x < 700
            r = 1;
            g = 0;
            b = 0;
            factor = 1;
        elseif x >= 700 && x <= 800
            r = 1;
            g = 0;
            b = 0;
            factor = 1 - (x - 700) / (800 - 700);
        end

        % Reduce intensity near the ends of the visible range.
        if x < 380
            factor = 0.3 + 0.7 * factor;
        elseif x > 700
            factor = max(0, factor);
        end

        rgb(k, :) = [r, g, b] * factor;
        rgb(k, :) = max(0, min(1, rgb(k, :)));
    end

    if isscalar(wavelength_nm)
        rgb = rgb(1, :);
    end
end

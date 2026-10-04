function cmap = rainbowCMAP(k)
%RAINBOWCMAP  Rainbow-plus-white colormap with k colors.
%
%   cmap = rainbowCMAP(k) returns a k-by-3 RGB colormap.
%
%   Special small cases:
%       k = 1 : white
%       k = 2 : green, red
%       k = 3 : blue, green, red
%       k = 4 : blue, green, red, white
%       k = 5 : blue, green, yellow, red, white
%   For k > 5 the colormap is obtained by interpolating along the rainbow
%   plus white:  violet, blue, turquoise, green, yellow, red, white.

    % ---- base colors (RGB, 0..1) ----
    white     = [1 1 1];
    red       = [1 0 0];
    yellow    = [1 1 0];
    green     = [0 1 0];
    turquoise = [0 1 1];
    blue      = [0 0 1];
    violet    = [0.5 0 1];

    switch k
        case 1
            cmap = white;
        case 2
            cmap = [green; red];
        case 3
            cmap = [blue; green; red];
        case 4
            cmap = [blue; green; yellow; red];
        case 5
            cmap = [blue; green; yellow; red; white];
        otherwise
            % k > 5 : interpolate along the full rainbow + white
            anchors = [violet; blue; turquoise; green; yellow; red; white];
            nAnchor = size(anchors, 1);

            % positions of the anchor colors and of the k output samples,
            % both on [0, 1]
            xAnchor = linspace(0, 1, nAnchor);
            xQuery  = linspace(0, 1, k);

            cmap = zeros(k, 3);
            for c = 1:3
                cmap(:, c) = interp1(xAnchor, anchors(:, c), xQuery, 'linear');
            end

            % guard against tiny numerical overshoot from interpolation
            cmap = min(max(cmap, 0), 1);
    end
end

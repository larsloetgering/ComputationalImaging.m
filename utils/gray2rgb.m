function r = gray2rgb( x, cmap, varargin )

p = inputParser;
% p.addParameter('intensityScale', [ ])
p.addParameter('sigma', [ ])
p.parse( varargin{:} );

if ~isempty(p.Results.sigma)
    x = x-mean(x(:));
    x = x/std(x(:));
    x(x > p.Results.sigma) = p.Results.sigma; 
    x(x < -p.Results.sigma) = -p.Results.sigma;
end

% x = x / ( max(x(:)) + eps );
% if ~isempty(p.Results.intensityScale)
%     x = not_negative(x - p.Results.intensityScale(1));
%     x = x/(p.Results.intensityScale(2) - p.Results.intensityScale(1));
%     x(x > p.Results.intensityScale(2))  = p.Results.intensityScale(2);
% end

cmap = linspace(0,1,255)' * cmap;
x = x - min(x(:));
x = x / max(max(x, [], 1), [], 2);
r = ind2rgb(uint8(round(x*255)), cmap);
end

function r = not_negative(x)
r = (abs(x) + x)/2;
end
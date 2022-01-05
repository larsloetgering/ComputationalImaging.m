function varargout = hsvplot_axis(x,y,u, varargin)
% hsvplot(u) generates hue-brightness plot of two dimensional input u 
% last change: 3rd January 2022

p = inputParser;
p.addParameter('intensityScale', [ ])
p.addParameter('pixelSize', [1,1])
p.parse( varargin{:} );

% fetch data from gpu (if needed)
u = gather(u);

% normalize birghtness (value) to range [0,1]
r = abs(u);
r = r / ( max(r(:)) + eps );

% readjust intensity maximum
if ~isempty(p.Results.intensityScale)
    r = posit(r - p.Results.intensityScale(1));
    r = r/(p.Results.intensityScale(2) - p.Results.intensityScale(1));
end

% normalize angle 
phi = angle( u );
phi = ( phi + pi )/( 2 * pi );

% normalization of phase saturation
B = zeros(size(r, 1), size(r, 2), 3, 'like', r);         % Declare RGB array
B(:,:,1) = phi;
B(:,:,2) = 1;
B(:,:,3) = r;
A = hsv2rgb(B);

imagesc(x,y,A); axis image 

set(gcf, 'Color', 'w');

switch nargout
    case 1
        varargout{1} = A;
    otherwise
end
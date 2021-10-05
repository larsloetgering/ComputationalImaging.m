function varargout = hsvplot(u, varargin)
% hsvplot(u) generates hue-brightness plot of two dimensional input u 
% last change: 3rd March 2018

p = inputParser;
p.addParameter('intensityScale', [ ])
p.addParameter('pixelsize', [1,1])
p.parse( varargin{:} );

% fetch data from gpu (if needed)
u = gather(u);

% set grid
N = size(A,2); x = (1:N) * p.Results.pixelsize(1);
M = size(A,1); y = (1:M) * p.Results.pixelsize(2);

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

if nargin == 1
    imagesc(A); axis image 
elseif nargin == 2
    imagesc(x,x,A); axis image 
elseif nargin == 3
    imagesc(x,y,A); axis image 
end
set(gcf, 'Color', 'w');

switch nargout
    case 1
        varargout{1} = A;
    otherwise
end
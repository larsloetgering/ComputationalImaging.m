function varargout = hsvplot(varargin)
% hsvplot(u) generates hue-brightness plot of two dimensional input u 
% alternative usage:
% hsvplot(x,u)
% hsvplot(x,y,u)
% last change: 3rd March 2018
% author: Lars Loetgering

if nargin == 1
    u = varargin{1};
elseif nargin == 2
    x = varargin{1};
    u = varargin{2};
elseif nargin == 3
    x = varargin{1};
    y = varargin{2};
    u = varargin{3};
end

u = gather(u);

% normalize birghtness (value) to range [0,1]
r = abs(u);
r = r / ( max(r(:)) + eps );

% normalize angle 
phi = angle( u );
phi = ( phi + pi )/( 2 * pi );

% normalization of phase saturation
B = zeros(size(r, 1), size(r, 2), 3, 'like', r); % Declare hsv array
B(:,:,1) = phi; % hue
B(:,:,2) = 1;   % saturation
B(:,:,3) = r;   % value
A = hsv2rgb(B); % hsv > rgb

if nargin == 1
    imagesc(A); axis image off
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
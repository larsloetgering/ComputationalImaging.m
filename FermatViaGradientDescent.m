addpath(genpath('utils'))

aleph = 1.5;
beth = 1e-2;
L = 2*pi/(sqrt(2*aleph/beth));
N = 100;
dx = L/N;
x = (-N/2:N/2-1)*dx;
% spatial frequency grid
df = 1/L;
f = (-N/2:N/2-1)*df;

numIter = 1e5;

% initialize y
% y = sin(2*pi*x/L);
y = rand(size(x));
stepSize = 1*1e-2;

for k = 1:numIter
    
%     yx = ifftc( fftc(y) .* (-1i*2*pi*f) );
%     yxx = ifftc( fftc(y) .*(-1i*2*pi*f).^2 );
    yx = [0,diff(y)];
    yxx = [0,0,diff(y,2)];
    
    denom = sqrt(1+yx.^2);
    yGrad = -2*beth*y.*sqrt(1+yx.^2) + ...
        2*beth*y.*yx.^2 ./sqrt(1+yx.^2) - ...
        (aleph - beth * y.^2).*(yxx.*(1+yx.^2) - yx.^2)./...
        (1+yx.^2).^(3/2);
    yGrad = yGrad ./ denom.^2;
%     yGrad = -2*beth*y + ...
%              2*beth*y.*yx./(1+yx.^2) - ...
%             (aleph - beth * y.^2).*(yxx.*(1+yx.^2) - yx.^2)./...
%                                     (1+yx.^2).^(2); 
                                
    y = y - stepSize * yGrad;
    y(1) = 0;
    y(end) = 0;
%     y = y/max(y);
    
    if mod(k,numIter/10)==0
    figure(1)
    plot(x,real(y))
    axis square
    title(num2str(k))
    drawnow
    end
end
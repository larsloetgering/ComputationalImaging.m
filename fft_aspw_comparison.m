clear
set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
addpath(genpath('utils'))

% N = 2^5;
N = 13;
dx = 50e-6;
L = N*dx;
x = linspace(-L/2,L/2,N);
wavelength = 633e-9;

% define aperture
w = L/2;
aperture = rect(x/w);

figure(1)
h = subplot(1,3,1);
plot(x,aperture,'k','lineWidth',2)
xlabel('x / m')
title('aperture')
axis square

z = 4e-3;
fx = (-N/2:N/2-1)/L;
f_max = L/(wavelength*sqrt(L^2+4*z^2));   
W = rect(fx/(2*f_max));
H = exp( 1i * 2*pi/wavelength * z * sqrt(complex(1 - (fx*wavelength).^2) )) .* W;
% H = exp( -1i * 2*pi/wavelength * z * (fx*wavelength).^2/2) .* W;

aperture_prop = ifftc( H.* fftc(aperture) );
subplot(1,3,2)
plot(x,abs(aperture_prop),'k','lineWidth',2)
xlabel('x / m')
title({'aperture','propagated'})
axis square

F = fftc(eye(N));
% ASP = ifftc( H.*fftc( eye(N) ) );
ASP = H.*fftc( eye(N) );

c = linspace(1,1,N);

% cmap = [hsv(size(F,2)/2); flipud(hsv(size(F,2)/2))].*gray(size(F,2));
cmap = hsv(size(F,2));
% cmap(1,:) = [0 0 0];
% cmap = hsv(size(F,2)/2);cmap = [cmap;flipud(cmap)];
figure(2), clf
subplot(1,2,1)
for k = 1:size(F,2)
    subplot(1,2,1)
    hold on
    plot(c(k)*real(F(:,k)),c(k)*imag(F(:,k)),'o-','color',cmap(k,:),...
            'MarkerFaceColor',cmap(k,:))
    axis square off
    hold off
    title(num2str(k))
%     axis(0.04*[-2 2 -2 2])
    
    subplot(1,2,2)
    hold on
    plot(c(k)*real(ASP(:,k)),c(k)*imag(ASP(:,k)),'o-','color',cmap(k,:),...
            'MarkerFaceColor',cmap(k,:))
    axis square off
    hold off
%     axis(0.04*[-2 2 -2 2])
    
    drawnow
    
    pause(0.1)
end

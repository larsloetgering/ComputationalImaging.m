fsize = 23;
set(0, 'DefaultAxesFontSize', fsize);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');
addpath(genpath('utils'))
clear

% preallocate source patterns
% N = 1001;
N = 257;
wavelength = 500e-9;
dx = wavelength/4;
df = 1/(N*dx);
f = (-(N-1)/2:(N-1)/2)*df;
% f = (-N/2:N/2-1)*df;
[Fx,Fy] = meshgrid(f);
NAi = linspace(0.2,0.8,4);

numz = 2*N+1;
dz = wavelength / (4*(1-sqrt(1-NAi(end)^2) ));
z = (-(numz-1)/2:(numz-1)/2)*dz/4;

numSources = length(NAi);
source = zeros(N,N,numSources);
for sourceLoop = 1:numSources
    % half pupil illumination
    % source(:,:,sourceLoop) = (circ(Fx,Fy,2*NAi(sourceLoop)/wavelength)).* (Fx>0);  % oblique illumination
    % azimuthal ring
    source(:,:,sourceLoop) = (1+cos(atan2(Fy,Fx)))/2 .* ...
                    (circ(Fx,Fy,2*NAi(sourceLoop)/wavelength) - ...
                             circ(Fx,Fy,4/5*2*NAi(sourceLoop)/wavelength));
    
end
pupil = circ(Fx,Fy,2*NAi(end)/wavelength);
%%
show_index = 4;
figure(1)
imagesc(source(:,:,show_index))
axis image off
colormap gray

figure(2), clf
imagesc(pupil)
axis image off
colormap gray

%% compute ATF and PTF

PTF = zeros(N,N,numz,numSources);
kz = 2*pi/wavelength*real(sqrt(1-(wavelength*Fx).^2-(wavelength*Fy).^2));
for sourceLoop = 1:numSources
    for zLoop = 1:numz
        % mixed representation
        propagator = exp(1i*z(zLoop)*kz);
        % [PTF_temp, ~] = getTransferFunctions(...
        %     source(:,:,sourceLoop), ...
        %     pupil.*propagator);
        [PTF_temp, ~] = getTransferFunctions(...
            imfilter(source(:,:,sourceLoop),fspecial('Gaussian',8,1)), ...
            imfilter(pupil.*propagator,fspecial('Gaussian',8,1)));
        PTF(:,:,zLoop,sourceLoop) = PTF_temp;
    end
end

figure(3), clf
imagesc(imag(PTF(:,:,(numz+1)/2,1)))
axis image off
colormap(slanCM('RdGy'))

%%

PTF3d = imag(fftshift(fft(ifftshift(PTF,3),[],3),3));
%%
PTF3d_crossSection = squeeze(PTF3d((N+1)/2,:,:,4));

figure(4), clf
% imagesc(log(abs(PTF3d_crossSection)+1e1).*sign(PTF3d_crossSection))
% imagesc(log(abs(PTF3d_crossSection)+1e0).*sign(PTF3d_crossSection))
% imagesc(nthroot(abs(PTF3d_crossSection).^2+1e-16,4).*sign(PTF3d_crossSection))
imagesc(imfilter(PTF3d_crossSection,fspecial('Gaussian',4,1)))
axis square off
% colormap(slanCM('fusion'))
% cmap = slanCM('seismic');
cmap = slanCM('RdGy');
colormap((cmap).^4)
% colormap(slanCM('PuOr'))
zoom(1.25)
% colormap('jet')
% colormap gray
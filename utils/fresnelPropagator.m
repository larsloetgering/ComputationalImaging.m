function [r, Qin, Qout] = fresnelPropagator(u, z, lambda, L)
% Fresnel propagator
% u = fresnelPropagator(u,z,lambda,L)
% u: field distribution at z = 0 
% (u is assumed to be square, i.e. N x N)
% z: propagation distance
% lambda: wavelength
% L: total size [m] of source field
% assume square array input, 
% even number of pixel per dimension

% get source (primed) coordinates
N = size(u,2);
dxp = L / N;
xp = (-N/2:N/2-1)*dxp;
[Xp,Yp] = meshgrid(xp);
% get target (unprimed) coordinates
dx = lambda * z / L;
x = (-N/2:N/2-1)*dx;
[X,Y] = meshgrid(x);
% get wavenumber
k = 2*pi/lambda;
% step 1: quadratic phase inside integral
Qin = exp(1i * k/(2*z) * (Xp.^2 + Yp.^2));
% step 2: Fourier transform
r = fft2c(u .* Qin);
% step 3: multiply quadratic phase 
% outside integral
Qout = exp(1i * k/(2*z) * (X.^2 + Y.^2));
r = Qout.*r;
end
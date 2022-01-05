function [u, H] = aspw(u, z, lambda, L)
% ASPW wave propagation
% u = aspw(u,z,lambda,L)
% u: field distribution at z = 0 
% (u is assumed to be square, i.e. N x N)
% z: propagation distance
% lambda: wavelength
% L: total size [m] of 
% following: Matsushima et al., 
% "Band-Limited Angular Spectrum Method for
% Numerical Simulation of Free-Space
% Propagation in Far and Near Fields", 
% Optics Express, 17 (22), pp. 19662-19673 
% (2009), https://doi.org/10.1364/OE.17.019662

k = 2*pi/lambda;
N = size(u,1);
[Fx, Fy] = meshgrid((-N/2:N/2-1)/L);
f_max = L/(lambda*sqrt(L^2+16*z^2));   
W = circ(Fx, Fy, 2*f_max); % low pass filter
H = exp(1i*k*z*sqrt(1-(Fx*lambda).^2-(Fy*lambda).^2));
u = ifft2c(bsxfun(@times, fft2c(u), H.*W));

end

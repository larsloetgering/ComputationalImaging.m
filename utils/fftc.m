function G = fftc(g)
% G = fft2(g)
% fftc performs centered 1-dimensional Fourier transformation
% fftc is normalized (i.e. norm(g) = norm(G) ), 
% i.e. it preserves the L2-norm

G = fftshift(fft(ifftshift(g))) / sqrt(numel( g ) );
% note that numel is only executed along a single slice, so that it is
% independent of the number of slices! (important for correct normalization 
% under partially coherent conditions)
end

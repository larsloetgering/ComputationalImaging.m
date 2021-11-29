function g = ifftc(G)
% G = ifftc(g)
% ifftc performs 2-dimensional inverse Fourier transformation
% ifftc is normalized (i.e. norm(g) = norm(G) ), 
% i.e. it preserves the L2-norm

% g = fftshift(ifft2(ifftshift(G))) * sqrt(numel(G));
g = fftshift(ifft(ifftshift(G))) * sqrt( numel(G) );
% note that numel is only executed along a single slice, so that it is
% independent of the number of slices! (important for correct normalization 
% under partially coherent conditions)
end

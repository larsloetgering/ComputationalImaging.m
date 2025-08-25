function r = mirror_pad(x, k)
% introduces semi-circular boundaries 
% (mitigates circular convolution artefacts)
% k: number of pixels padded on each side
r = zeros(size(x,1) + 2*k, size(x,2) + 2*k, size(x,3));
for loop = 1:size(x, 3)
    y = x(:,:,loop);
    z = [fliplr( y(:,1:k) ), y, fliplr( y(:,end-k+1:end) )];
    r(:,:,loop) = [flipud( z(1:k,:) ); z; flipud( z(end-k+1:end,:) )];
end
end
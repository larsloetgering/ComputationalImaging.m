function r = crop(x, px)
% think of this as the adjoint of mirror_pad and zero_pad
r = x(px+1:end-px, px+1:end-px,:,:);

end
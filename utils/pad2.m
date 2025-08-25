function r = pad2( psi )
% input: psi, 2D array with even number of pixels per dimension
[M,N] = size(psi);
r = circshift( padarray( psi, [M N],'post'), [M N]/2);
end
function [W,H,residual] = L2_NMF(A,r,numIter)
% NMF to minimize norm(A-W*H,'fro') wrt W,H
% input:
% A: MxN matrix (images arranged as row vectors)
% r: rank
% numIter: number of iterations
% output:
% W: interpretable output images
% H: weights

% initialize by absolute value of each 
% truncated SVD factors
[U,S,V] = svd(A,'econ');
W = (abs(U(:,1:r)) + 1)/2;
H = (abs(S(1:r,1:r)*V(:,1:r)') + 1);
delta = 1;      % prevent division by zero

for k = 1:numIter
    H = H .* (W'*A + delta)./((W'*W)*H + delta);
    W = W .* ((A*H' + delta)./(W*(H*H') + delta));
end

residual = norm(A - W*H,'fro') / norm(A, 'fro');
end
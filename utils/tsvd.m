function [U, S, V] = tsvd(A, r)
% truncated svd
[U, S, V] = svd(A,'econ');
U = U(:,1:r);
S = S(1:r,1:r);
V = V(:,1:r);

end
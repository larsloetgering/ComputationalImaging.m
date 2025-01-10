function [Y, Lambda] = DMD(X, r)
% dynamic mode decomposition of a time series
% DMD model: X(:,2:end) = A*X(:,1:end-1) 
% input:
% X: input data (entire time series M x N)
% r: rank
% output:
% Y: high-dimensional eigenmodes of A
% Lambda: oscillation/decay rates

% step 1: truncated SVD
[U,S,V] = svd(X(:,1:end-1),'econ');
Ur = U(:,1:r);
Sr = S(1:r,1:r);
Vr = V(:,1:r);
% Step 2: compute projected matrix
Aprojected = Ur' * X(:,2:end) * (Vr/Sr); 
% Step 3: eigendecomposition
[Z,Lambda] = eig(Aprojected);  
% step 4: high-dimensional eigenmodes
Y = X(:,2:end)*(Vr/Sr)*Z;
end
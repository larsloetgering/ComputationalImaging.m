% gradient descent with matrix valued derivatives
clc

set(0, 'DefaultAxesFontSize', 24);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')

rng(0)
m = 5;
n = 5;
X = rand(m,n);
B = cat(1, X, X);
% B = flipud(X);

Alsq = B*X'/(X*X');
% gradient descent
numIter = 1e4;
A = rand(size(B));
% A = Alsq;
aleph = 5e-2;
for k = 1:numIter
    A = A - aleph * (A*X-B)*X';
end

figure(1)
subplot(1,2,1)
imagesc(A)
axis image off, colormap gray
title('gradient descent')

subplot(1,2,2)
imagesc(Alsq)
axis image off, colormap gray
title('least square [via normal Eq.]')


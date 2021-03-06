% gradient descent with matrix valued derivatives
clc
addpath(genpath('utils'))
set(0, 'DefaultAxesFontSize', 24);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')

rng(0)
m = 5;
n = 10;
X = rand(m,n);
B = cat(1, X, flipud(X));

Alsq = B*X'/(X*X');
% gradient descent
numIter = 1e3;
cost = zeros(numIter,1);
A = rand(size(B,1),size(X,1));
% A = Alsq;
stepSize = 1e-1;
R = A*X - B;
cost(1) = norm(R,'fro');
for loop = 2:numIter
    % gradient descent step
    A_test = A - stepSize * R*X';
    % test residual
    R = A_test*X - B;
    % test cost
    cost_test = norm(R,'fro');
    % backtracking
    if cost_test <= cost(loop-1)
        cost(loop) = cost_test;
        A = A_test;
    else
        stepSize = 0.9*stepSize;
        cost(loop) = cost(loop-1);
        disp(loop)
    end
end

figure(1)
subplot(1,3,1)
imagesc(A)
axis image off, colormap winter
title('gradient descent')

subplot(1,3,2)
imagesc(Alsq)
axis image off, colormap winter
title('lsq [via normal Eq.]')

subplot(1,3,3)
loglog(cost,'k','linewidth',2)
axis square, title('cost')
xlabel('iteration'), grid on


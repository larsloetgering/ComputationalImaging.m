%% configuration
clear
fsize = 32;
set(0, 'DefaultAxesFontSize', fsize);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');

%% fixed parameters of loss function

A = [-1/2 1];
b = 1;

%% initialization
clc
% rng(0);
W = diag(rand(2,1));
x = rand(2,1);
numIter = 5e1;
delta = 1e-7;
lambda = 1;
mu = 1e-7;

for loop = 1:numIter
    % step 1: update W
    xold = x;

    % step 2: update x
    x = inv(A'*A + lambda*W'*W + mu*eye(2)) * (A'*b + mu*x);
    
    % step 3: update W
    W = diag( 1 ./ max( sqrt( abs(xold) ), delta ) );

    % step 4: decrease lambda
    lambda = lambda * 0.9;
end
disp(x)



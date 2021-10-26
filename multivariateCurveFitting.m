clc
clear, close all
set(0, 'DefaultAxesFontSize', 24);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
rng(0)
N = 3e1;
noise = 0.1*randn(1,N);
x = linspace(0,1,N);
a = -1; b = 1;
d = a*x + b + noise;
figure(1)
plot(x,d,'ko','MarkerFaceColor','k')
axis square
xlabel('x'), ylabel('data'), grid on
%% choose initial parameters
a = -2;
b = 0;
%% gradient descent
numIter = 100;
stepSize = 1e-3;
figureUpdate = 10;
cost = zeros(numIter,1);
for loop = 1:numIter
    % compute residual
    r = 2*(a*x+b-d);
    % compute gradient
    aGrad = sum(r .* x);
    bGrad = sum(r);
    % search along negative gradient
    a = a - stepSize * aGrad;
    b = b - stepSize * bGrad;
    % get cost
    cost(loop) = sum( r.^2 );
    % show current estimate
    if mod(loop,figureUpdate)==0
        figure(2)
        subplot(1,2,1)
        plot(x,d,'ko','MarkerFaceColor','k')
        hold on
        plot(x,a*x+b,'bo-','MarkerFaceColor','b','linewidth',2)
        legend('data', 'fit')
        axis square, xlabel('x'), ylabel('data'), grid on
        hold off
        
        subplot(1,2,2)
        loglog(1:loop, cost(1:loop), 'ko','MarkerFaceColor','k')
        axis square, grid on
        xlabel('iteration'), ylabel('cost')
        drawnow
    end
end


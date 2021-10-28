clc
clear, close all
set(0, 'DefaultAxesFontSize', 24);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
rng(0)
N = 3e1;
noise = 0.15*randn(1,N);
x = linspace(0,1,N);
a = -1; b = 1;
d = a*x + b + noise;
figure(1)
plot(x,d,'ko','MarkerFaceColor','k')
axis square
xlabel('x'), ylabel('data'), grid on
%% choose initial parameters
a = a+1;
b = b+1;
a_bt = a;
b_bt = b;
%% gradient descent
numIter = 100;
stepSize = 1e-3;
stepSize_bt = 1;
figureUpdate = 100;
cost = zeros(numIter,1);
cost_bt = zeros(numIter,1);
% compute residuals
r = 2*(a*x-b-d);
r_bt = 2*(a_bt*x-b_bt-d);
% compute initial cost
cost(1) = sum( r.^2 );
cost_bt(1) = sum( r_bt.^2 );
%%
for loop = 2:numIter
    % compute gradients
    aGrad = sum(r .* x);
    bGrad = sum(r);
    a_btGrad = sum(r_bt .* x);
    b_btGrad = sum(r_bt);
    % search along negative gradient
    a = a - stepSize * aGrad;
    b = b - stepSize * bGrad;
    a_test = a_bt - stepSize_bt * a_btGrad;
    b_test = b_bt - stepSize_bt * b_btGrad;
    % compute new residuals
    r = 2*(a*x+b-d);
    r_bt_test = 2*(a_test*x+b_test-d);
    % get cost
    cost(loop) = sum( r.^2 );
    cost_bt_test = sum( r_bt_test.^2 );
    if cost_bt_test < cost_bt(loop-1)
        cost_bt(loop) = cost_bt_test;
        a_bt = a_test;
        b_bt = b_test;
        r_bt = r_bt_test;
    else
        stepSize_bt = 2/3*stepSize_bt;
        cost_bt(loop) = cost_bt(loop-1);
    end
    % show current estimate
    if mod(loop,figureUpdate)==0
        figure(2)
        subplot(1,2,1)
        plot(x,d,'ko','MarkerFaceColor','k')
        hold on
        plot(x,a*x+b,'bo-','MarkerFaceColor','b','linewidth',2)
        plot(x,a_bt*x+b_bt,'go-','MarkerFaceColor','g','linewidth',2)
        legend('data', 'gradient descent','backtracking')
        axis square, xlabel('x'), ylabel('data'), grid on
        hold off, title(['iteration: ',num2str(loop)])
        
        subplot(1,2,2)
        loglog(1:loop, cost(1:loop), 'bo','MarkerFaceColor','b')
        hold on
        loglog(1:loop, cost_bt(1:loop), 'go','MarkerFaceColor','g')
        axis square, grid on
        xlabel('iteration'), ylabel('cost')
        hold off
        drawnow
    end
end

disp(cost(end))
disp(cost_bt(end))
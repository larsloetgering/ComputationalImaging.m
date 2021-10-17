clc
clear
set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
rng(0)
N = 5e2;
x = randn(N,1);
K = 20;
[d,t] = histcounts(x, K);
t = t(1:K) + 0.5;
figure(1)
plot(t,d,'ko','MarkerFaceColor','k')
xlabel('x value'), ylabel('counts'), grid on
%% choose initial parameters
A = max(d);
B = sqrt(sum(t.^2.*d/sum(d)));
C = 0;
%% gradient descent
numIter = 100;
stepSize = 1e-4;
figureUpdate = 10;
cost = zeros(numIter,1);
for loop = 1:numIter
    % compute gradient
    r = 2*(log(A)-B*(t-C).^2-log(d));
    Agrad = sum(r)/A;
    Bgrad = -sum(r.*(t-C).^2);
    Cgrad = -2*sum(r.*B.*(C-t));
    % search along negative gradient
    A = A - stepSize * Agrad;
    B = B - stepSize * Bgrad;
    C = C - stepSize * Cgrad;
    % get cost
    cost(loop) = sum( (log(A) - B*(t-C).^2 -log(d)).^2 );
    % show current estimate
    if mod(loop,figureUpdate)==0
        figure(2)
        subplot(1,2,1)
        plot(t,d,'ko','MarkerFaceColor','k')
        hold on
        plot(t,A*exp(-B*(t-C).^2),'bo','MarkerFaceColor','b')
        legend('data', 'fit')
        axis square, xlabel('t'), ylabel('counts'), grid on
        hold off, title(['iteartion: ',num2str(loop)])
        
        subplot(1,2,2)
        loglog(1:loop, cost(1:loop), 'ko','MarkerFaceColor','k')
        axis square, grid on
        xlabel('iteartion'), ylabel('cost')
        drawnow
    end
end


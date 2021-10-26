clc
clear
set(0, 'DefaultAxesFontSize', 24);
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
A = 0.99*max(d);
B = sqrt(sum(t.^2.*d/sum(d)));
C = 0;
%% gradient descent
numIter = 100;
stepSize = 1e-5;
for loop = 1:numIter
    % compute gradient
    Agrad = 2*sum( (A*exp(-B*(t-C).^2) - d) .* exp(-B*(t-C).^2) );
    Bgrad = 2*sum( (A*exp(-B*(t-C).^2) - d) .* A.* exp(-B*(t-C).^2).*(t-C).^2 .* (-1) );
    Cgrad = 2*sum( (A*exp(-B*(t-C).^2) - d) .* A.* exp(-B*(t-C).^2).*B.*(t-C) );
    % search along negative gradient
    A = A - stepSize * Agrad;
    B = B - stepSize * Bgrad;
    C = C - stepSize * Cgrad;
    % show current estimate
    figure(2)
    plot(t,d,'ko','MarkerFaceColor','k')
    hold on
    plot(t,A*exp(-B*(t-C).^2),'bo','MarkerFaceColor','b')
    legend('data', 'fit')
    hold off, title(['iteartion: ',num2str(loop)]), drawnow
end


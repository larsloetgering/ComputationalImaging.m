%% configuration
fsize = 32;
set(0, 'DefaultAxesFontSize', fsize);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');

%% 

xx = linspace(0,2,100);
f1 = xx;
f2 = xx + 1;
f3 = -0.5*xx + 1;
f4 = -0.5*xx + 2;

figure(1), clf
plot(xx,f1,'k-','lineWidth',2)
hold on
plot(xx,f2,'-','lineWidth',2, 'color', [0.7 0 0])
plot(xx,f3,'-','lineWidth',2, 'color', [0 0.7 0])
plot(xx,f4,'-','lineWidth',2, 'color', [0 0 0.7])

axis square
grid on
xlabel('$x_1$'), ylabel('$x_2$')

%%

A1 = [-1 1];
A2 = [-1 1];
A3 = [0.5 1];
A4 = [0.5 1];
b1 = 0;
b2 = 1;
b3 = 1;
b4 = 2;

rsol = inv(A1'*A1 + A2'*A2 + A3'*A3+A4'*A4)*...
          (A1'*b1+A2'*b2+A3'*b3+A4'*b4);

%% initialization

rng(1)
a = 2.5*rand(2,1);

%% AD algorithm

x = a;
numIter = 100;
lambda = 1e-4;
plot(x(1),x(2),'kd','lineWidth',2,'MarkerSize',7)
plot(rsol(1), rsol(2), 'co', 'MarkerFaceColor', [0 1 1])

for loop = 1:numIter
    xold = x;
    % step 1: b1=0, so a bit redundant (in favor of readibility and context)
    x = inv(A1'*A1 + lambda*eye(2))*(A1'*b1 + lambda*x);
    plot([xold(1) x(1)],[xold(2) x(2)],'k-','lineWidth',2)
    plot(x(1),x(2),'ko','lineWidth',2,'MarkerSize',7,...
        'MarkerFaceColor','k')
    % step 2: 
    xold = x;
    x = inv(A2'*A2 + lambda*eye(2))*(A2'*b2 +lambda*x);
    
    plot([xold(1) x(1)],[xold(2) x(2)],'k-','lineWidth',2)
    plot(x(1),x(2),'ko','lineWidth',2,...
        'MarkerSize',7,'MarkerFaceColor','k')

    x = inv(A3'*A3 + lambda*eye(2))*(A3'*b3 +lambda*x);
    
    plot([xold(1) x(1)],[xold(2) x(2)],'k-','lineWidth',2)
    plot(x(1),x(2),'ko','lineWidth',2,...
        'MarkerSize',7,'MarkerFaceColor','k')

    x = inv(A4'*A4 + lambda*eye(2))*(A4'*b4 +lambda*x);
    
    plot([xold(1) x(1)],[xold(2) x(2)],'k-','lineWidth',2)
    plot(x(1),x(2),'ko','lineWidth',2,...
        'MarkerSize',7,'MarkerFaceColor','k')
end
axis([0 2 0 2])
% h = legend('$x_2 = x_1$','$x_2 = 0.5x_1 + 0.5$');
% h.Location = 'northwest';



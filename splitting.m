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
f2 = 0.5*xx+0.5;

figure(1), clf
plot(xx,f1,'b-','lineWidth',2)
hold on
plot(xx,f2,'-','lineWidth',2, 'color', [0 0.7 0])
axis square
grid on
xlabel('$x_1$'), ylabel('$x_2$')

%%

A1 = [-1 1];
A2 = [-0.5 1];
b1 = 0;
b2 = 0.5;

%% splitting algorithm

x = [0.25;1];
y = x;
numIter = 100;
mu = 1;
gamma = 1.1;
plot(x(1),x(2),'kx','lineWidth',2,'MarkerSize',7)
plot(y(1),y(2),'mx','lineWidth',2,'MarkerSize',7)
for loop = 1:numIter
    xold = x;
    yold = y;
    % step 1: b1=0, so a bit redundant (in favor of readibility and context)
    x = inv(A1'*A1 + mu*eye(2))*(A1'*b1 + mu*y);

    % step 2: 
    y = inv(A2'*A2 + mu*eye(2))*(A2'*b2 + mu*x);
    
    plot([xold(1) x(1)],[xold(2) x(2)],'k-','lineWidth',2)
    plot([yold(1) y(1)],[yold(2) y(2)],'m-','lineWidth',2)
    
    plot(x(1),x(2),'ko','lineWidth',2,'MarkerSize',7,'MarkerFaceColor','k')
    plot(y(1),y(2),'mo','lineWidth',2,'MarkerSize',7,'MarkerFaceColor','m')
    
    
end

% h = legend('$x_2 = x_1$','$x_2 = 0.5x_1 + 0.5$');
% h.Location = 'northwest';



fsize = 32;
set(0, 'DefaultAxesFontSize', fsize);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');
clear

%% penalty method with augmented Lagrangian 
% (set rho = 0 to retrieve simple solver)

t = linspace(0,pi,100);
x1 = linspace(-4,4,100);
y1 = 1/2*x1 + 2; A = [-1/2 1]; b = 2;

figure(1), clf

hold on
cmap = gray(40);
cmap(1:10,:) = [];
for k = 1:size(cmap,1)
    c = 0.2;
    xc = c*k*cos(t);
    yc = c*k*sin(t);
    plot(xc, yc,'-','color',cmap(k,:),'lineWidth',2)
end
plot(x1,y1,'k--','lineWidth',2)
grid on
axis equal
fill([-4 0 0 -4], [0 0 4 4], [0.5 0.5 0.5],'FaceAlpha',0.3)
axis([-4 4 0 4])
xlabel('$x_1$')
ylabel('$x_2$')

% solve problem via Lagrange multipliers

numIter = 200;
rng(0)
x = 3*rand(2,1);
z = 3*rand(2,1);
lambda = rand(1,1);
stepSize = 1/10;
mu = rand(2,1);
rho = 0;
gamma = 1.01;

for k = 1:numIter
    figure(1)
    hold on
    plot(x(1),x(2),'kx','lineWidth',2,'MarkerSize',10)

    xGrad = (2*x + A'*lambda - mu);
    
    % step 1: x gradient (argmin)
    x = x - stepSize * xGrad;

    % step 2: Lagrange multiplier update (argmax)
    lambdaGrad = A*x - b;
    lambda = lambda + stepSize*lambdaGrad;
    
    % step 3: mu update (argmax)
    muGrad = -x;
    mu = mu + stepSize * muGrad;
    mu = max(mu, 0);
end

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

%% splitting algorithm

% a = [0.5; 2];
% b = [1; 2];
% c = [1.5;2];
% d = [2;2];

rng(5)
a = rand(2,1);
b = rand(2,1);
c = rand(2,1);
d = rand(2,1);
numIter = 10;
mu = 1e-1;
gamma = 1.5;

plot(a(1),a(2),'kd','lineWidth',2,'MarkerSize',7,'MarkerFaceColor','k')
plot(b(1),b(2),'rd','lineWidth',2,'MarkerSize',7,'MarkerFaceColor','r')
plot(c(1),c(2),'gd','lineWidth',2,'MarkerSize',7,'MarkerFaceColor','g')
plot(d(1),d(2),'bd','lineWidth',2,'MarkerSize',7,'MarkerFaceColor','b')

plot(rsol(1), rsol(2), 'co', 'MarkerFaceColor', [0 1 1])
for loop = 1:numIter

    aold = a;
    % step 1: b1=0, so a bit redundant (in favor of readibility and context)
    a = inv(A1'*A1 + mu*eye(2))*(A1'*b1 + mu*b);

    bold = b;
    % step 2: 
    b = inv(A2'*A2 + mu*eye(2))*(A2'*b2 + mu*c);
    
    cold = c;
    % step 3: 
    c = inv(A3'*A3 + mu*eye(2))*(A3'*b3 + mu*d);

    dold = d;
    % step 4:
    d = inv(A4'*A4 + mu*eye(2))*(A4'*b4 + mu*a);
    
    mu = mu*gamma;

    plot([aold(1) a(1)],[aold(2) a(2)],'k-','lineWidth',2)
    plot([bold(1) b(1)],[bold(2) b(2)],'r-','lineWidth',2)
    plot([cold(1) c(1)],[cold(2) c(2)],'g-','lineWidth',2)
    plot([dold(1) d(1)],[dold(2) d(2)],'b-','lineWidth',2)

    plot(a(1),a(2),'ko','lineWidth',2,'MarkerSize',7,'MarkerFaceColor','k')
    plot(b(1),b(2),'ro','lineWidth',2,'MarkerSize',7,'MarkerFaceColor','r')
    plot(c(1),c(2),'go','lineWidth',2,'MarkerSize',7,'MarkerFaceColor','g')
    plot(d(1),d(2),'bo','lineWidth',2,'MarkerSize',7,'MarkerFaceColor','b')

    % plot(x(1),x(2),'ko','lineWidth',2,'MarkerSize',7)
    % plot(y(1),y(2),'mo','lineWidth',2,'MarkerSize',7)
    % plot(z(1),z(2),'yo','lineWidth',2,'MarkerSize',7)
    
end
axis([0 2 0 2])
% h = legend('$x_2 = x_1$','$x_2 = 0.5x_1 + 0.5$');
% h.Location = 'northwest';



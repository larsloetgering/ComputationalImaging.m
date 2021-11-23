set(0, 'DefaultAxesFontSize', 24);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
set(0, 'DefaultLineLineWidth', 2);
clear, close all

rng(0)  % seed random number generator
N = 30; % number of data points (in total)
M = round(N/5); % number of test points for cross-validation

x = linspace(0,1,N)';       % data grid
xd = linspace(0,1,5*N)';    % dense grid
indices = randperm(N);      % randomized indices (for test data)
z = x(indices(1:M));        % cross-validation grid  
x(indices(1:M)) = [ ];      % reduced data grid
y = exp(x);                 % ground truth
data= y + 1/2*randn(N-M,1);   % noisy observation
dataz = exp(z) + 1/2*randn(M,1); % noisy observation for cross-validation

figure(1)
plot(x,data,'ko','MarkerFaceColor','k')
axis square

p = 10;
X = zeros(size(x,1), p+1);
Xd = zeros(size(xd,1),p+1);
Z = zeros(size(z,1), p+1);
for k = 0:p
    X(:,k+1) = x.^k;
    Xd(:,k+1) = xd.^k;
    Z(:,k+1) = z.^k;
end
% least squares
c_lsq = (X'*X)\(X'*data);
figure(1)
hold on
plot(xd,Xd*c_lsq,'-b')
hold off

% ridge regression
lambda = logspace(-2,2,20);
C_ridge = zeros(length(c_lsq),length(lambda));
for k = 1:length(lambda)
    C_ridge(:,k) = (X'*X + lambda(k)*eye(p+1,p+1))\(X'*data);
end

cmap = 0.9*hot(p+1);
figure(2), clf
subplot(1,2,1)
for k = 1:size(C_ridge,1)
    semilogx(lambda, C_ridge(k,:),'color',cmap(k,:))
    hold on
end
hold off, grid on, axis square
hold on
GT = 1./(gamma(linspace(0,p,p+1)+1))' * ones(1,length(lambda));
for k = 2:4
    semilogx(lambda, GT(k,:),'--','color',cmap(k,:))
    hold on
end
xlabel('\lambda')
legend('c_0','c_1','c_2','c_3','c_4','c_5','c_6','c_7','c_8','c_9','c_{10}',...
       'gt_1','gt_2','gt_3')
legend('location','North','Numcolumns',11)

% compute ground truth
cost_GT = zeros(1,length(lambda));
cost_CV = zeros(1,length(lambda));
for k = 1:size(C_ridge,2)
    cost_GT(k) = norm( C_ridge(:,k) - GT(:,k),2)^2;
    cost_CV(k) = norm( Z*C_ridge(:,k) - dataz,2)^2;
end

cost_GT = cost_GT/max(cost_GT);
cost_CV = cost_CV/max(cost_CV);

color1 = 0.5*[1 0 1];
color2 = 0.5*[0 1 0];
h = subplot(1,2,2);
yyaxis left
semilogx(lambda, cost_GT, '-o','color',color1,'MarkerFaceColor',color1)
h.YColor = color1;
ylabel('ground truth')
hold on
yyaxis right
semilogx(lambda, cost_CV, '-o','color',color2,'MarkerFaceColor',color2)
h.YColor = color2;
ylabel('cross-validation')
grid on, axis square, grid on, xlabel('\lambda')
legend('ground truth cost','cross-validation')
legend('location','North')

[~,minIndex] = min(cost_CV);
optimalLambda = lambda(minIndex);
disp(['optimal lambda: ', num2str(optimalLambda)])
figure(1)
hold on 
plot(xd, Xd*C_ridge(:,minIndex),'-r')
plot(x,y,'k--')
hold off
legend('data','LSQ','ridge','ground truth')
legend('location','SouthEast')
xlabel('x')
grid on
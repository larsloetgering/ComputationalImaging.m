set(0, 'DefaultAxesFontSize', 24);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
set(0, 'DefaultLineLineWidth', 2);
clear 

N = 50;                     % number of data points
rng(1)
x = linspace(0,1,N)';       % data grid
xd = linspace(0,1,5*N)';    % dense grid
y = exp(x);                 % ground truth
data = y + 1e-1*randn(N,1);  % noisy observation

figure(1)
plot(x,data,'ko','MarkerFaceColor','k')
axis square
%%
p = 10;
X = zeros(size(x,1),p+1);
Xd = zeros(size(xd,1),p+1);
for k = 0:p
    X(:,k+1) = x.^k;
    Xd(:,k+1) = xd.^k;
end
% least squares
c_lsq = (X'*X)\(X'*data);
% ridge regression
lambda = 1e0;
c_ridge = (lambda * X'*X + eye(p+1,p+1))\(X'*data);

% duality
c_aug = rand(size(c_ridge));
lambda_aug = rand(size(X,1),1);
numIter = 40;
for k = 1:numIter
    % compute c search step
    c_aug = -1/2*X'*lambda_aug;
    % improve lambda_aug 
    lambda_aug = lambda_aug - 1/100*(X*c_aug-data);
end

figure(2),clf
subplot(2,2,1)
bar(0:length(c_lsq)-1,c_lsq,'FaceColor','b')
xlabel('c_{k;lsq}'), axis square, title('LSQ')

subplot(2,2,2)
bar(0:length(c_ridge)-1,c_ridge,'FaceColor','r')
xlabel('c_{k;ridge}'), axis square, title('ridge')

subplot(2,2,3)
bar(0:length(c_aug)-1,c_aug,'FaceColor','g')
xlabel('c_{k;aug}'), axis square, title('augmented')

subplot(2,2,4)
bar(0:length(c_ridge)-1,1./gamma(1:length(c_ridge)),'FaceColor','k')
axis square, title('ground truth')
xlabel('c_{k}')

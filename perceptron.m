set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultFigureColor', 'w');
addpath(genpath('utils'))
clear, close all

K = 50;
rng(1)

% example 1 (no covariance)
mu1 = [5;0];
mu2 = [15;5];
sigma1 = [3,0;0,3];
sigma2 = [3,0;0,3];

% example 2 (with covariance)
% mu1 = [5;0];
% mu2 = [5;5];
% p = 0.9;
% sigma1 = 10*[1,p;p,1];
% sigma2 = 10*[1,p;p,1];

% training data
X1 = mvnrnd(mu1,sigma1,K)';
X2 = mvnrnd(mu2,sigma2,K)';
% validation data
X1val = mvnrnd(mu1,sigma1,K)';
X2val = mvnrnd(mu2,sigma2,K)';

%% convert to training vector

Xtrain = [X1, X2];
Xtrain = [ones(1,2*K); Xtrain];
ytrain = [-1*ones(1,50), 1*ones(1,50)];

%% initialize decision boundary

w = rand(size(Xtrain,1), 1);

% subgradient perceptron learning algorithm

numIter = 1000;
stepSize = 1;
figureUpdate = 10;
for k = 1:numIter
    % sgPLA
    idx = randi(2*K);
    
    subGrad = -(ytrain(idx) - sign(w'*Xtrain(:,idx)))*Xtrain(:,idx);
    w = w - stepSize * subGrad;
    
    if mod(k,figureUpdate) == 0
        showResults(Xtrain,X1,X2,w,idx)
        title(['iteration: ',num2str(k)])
%         pause(0.5)
    end
end



function showResults(Xtrain,X1,X2,w,idx)
%% generate decision boundary

s1 = linspace(min([X1(1,:),X2(1,:)]),max([X1(1,:),X2(1,:)]),100);
s2 =  -w(2)/w(3)*s1 - w(1)/w(3);

%% show results

figure(1), clf
scatter(X1(1,:),X1(2,:),'mo','filled'), hold on
scatter(X2(1,:),X2(2,:),'go','filled')
plot(s1,s2,'--k','LineWidth',2)
plot(Xtrain(2,idx),Xtrain(3,idx),'ok','LineWidth',2)
minAx1 = min([X1(1,:),X2(1,:)]);
minAx2 = min([X1(2,:),X2(2,:)]);
maxAx1 = max([X1(1,:),X2(1,:)]);
maxAx2 = max([X1(2,:),X2(2,:)]);
axis([minAx1, maxAx1, minAx2, maxAx2])
xlabel('feature 1'), ylabel('feature 2')
grid on, axis square
hold off
end
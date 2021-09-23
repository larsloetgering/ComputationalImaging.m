set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultFigureColor', 'w');
addpath(genpath('utils'))
clear, close all
K = 50;
rng(1)
mu1 = [5;0];
mu2 = [6;5];
% example 1
sigma1 = [3,0;0,3];
sigma2 = [3,0;0,3];
% example 2
% sigma1 = [3,3;3,4];
% sigma2 = [3,3;3,4];
% training data
X1 = mvnrnd(mu1,sigma1,K)';
X2 = mvnrnd(mu2,sigma2,K)';
% validation data
X1val = mvnrnd(mu1,sigma1,K)';
X2val = mvnrnd(mu2,sigma2,K)';

m1 = mean(X1,2);
m2 = mean(X2,2);
% between class separation matrix
Sb = (m2 - m1)*(m2 - m1)';
% within class covariance
Sw = (X1 - m1)*(X1 - m1)' + ...
     (X2 - m2)*(X2 - m2)';
% separation surface normal
w = Sw\(m2 - m1);
wlsq = m2 - m1;
% global mean
m = 1/2*(m1 + m2);
% threshold
c = -w'*m;
clsq = (m1 - m2)'*m;
% generate 
s1 = linspace(min([X1(1,:),X2(1,:)]),max([X1(1,:),X2(1,:)]),100);
s2 = -1/w(2)*(c+w(1)*s1);
slsq = -1/wlsq(2)*(clsq+wlsq(1)*s1);

figure(1), clf
scatter(X1(1,:),X1(2,:),'mo','filled'), hold on
scatter(X2(1,:),X2(2,:),'go','filled')
plot(s1,s2,'--k','LineWidth',2)
plot(s1,slsq,'-k','LineWidth',2)
minAx1 = min([X1(1,:),X2(1,:)]);
minAx2 = min([X1(2,:),X2(2,:)]);
maxAx1 = max([X1(1,:),X2(1,:)]);
maxAx2 = max([X1(2,:),X2(2,:)]);
axis([minAx1, maxAx1, minAx2, maxAx2])
xlabel('feature 1'), ylabel('feature 2')
grid on, axis square
hold off

figure(2), clf
scatter(X1val(1,:),X1val(2,:),'md','filled'), hold on
scatter(X2val(1,:),X2val(2,:),'gd','filled')
plot(s1,s2,'--k','LineWidth',2)
plot(s1,slsq,'-k','LineWidth',2)
axis([minAx1, maxAx1, minAx2, maxAx2])
xlabel('feature 1'), ylabel('feature 2')
grid on, axis square
hold off

figure(3)
subplot(1,2,1)
bar( 1:K, sign(w'*X1val+c),'FaceColor','m','EdgeColor','w' ), hold on
bar( K+1:2*K, sign(w'*X2val+c),'FaceColor','g','EdgeColor','w' ), hold on
title('LDA validation')
axis square, hold off
axis([0 2*K -1 1])
xticks([1 50 100]), yticks([-1 1])
ylabel('class')

subplot(1,2,2)
bar( 1:K, sign(wlsq'*X1val+clsq),'FaceColor','m','EdgeColor','w' ), hold on
bar( K+1:2*K, sign(wlsq'*X2val+clsq),'FaceColor','g','EdgeColor','w' ), hold on
title('LSQ validation')
axis square, hold off
axis([0 2*K -1 1])
xticks([1 50 100]), yticks([-1 1])
ylabel('class')
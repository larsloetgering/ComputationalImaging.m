set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultFigureColor', 'w');
% clear, close all
K = 50;
rng(1)
mu1 = [10;0];
mu2 = [20;10];
% example 1
% M1 = [2,0;0,2];
% M2 = [2,0;0,2];
% example 2
M1 = [2,2;1,2];
M2 = [2,2;1,2];
% training data
X1 = mu1+M1*randn(2,K);
X2 = mu2+M2*randn(2,K);
% validation data
X1val = mu1+M1*randn(2,K);
X2val = mu2+M2*randn(2,K);

figure(1), clf
scatter(X1(1,:),X1(2,:),'mo','filled'), hold on
scatter(X2(1,:),X2(2,:),'go','filled')

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

plot(s1,s2,'--k','LineWidth',2)
plot(s1,slsq,'-k','LineWidth',2)
minAx1 = min([X1(1,:),X2(1,:)]);
minAx2 = min([X1(2,:),X2(2,:)]);
maxAx1 = max([X1(1,:),X2(1,:)]);
maxAx2 = max([X1(2,:),X2(2,:)]);
axis([minAx1, maxAx1, minAx2, maxAx2])
grid on, axis square
hold off

figure(2), clf
scatter(X1val(1,:),X1val(2,:),'md'), hold on
scatter(X2val(1,:),X2val(2,:),'gd')
plot(s1,s2,'--k','LineWidth',2)
plot(s1,slsq,'-k','LineWidth',2)
axis([minAx1, maxAx1, minAx2, maxAx2])
grid on, axis square
hold off

figure(3)
subplot(1,2,1)
bar( 1:K, sign(w'*X1val+c),'FaceColor','m','EdgeColor','w' ), hold on
bar( K+1:2*K, sign(w'*X2val+c),'FaceColor','g','EdgeColor','w' ), hold on
title('LDA validation')
axis square, hold off
axis([0 2*K -1 1])

subplot(1,2,2)
bar( 1:K, sign(wlsq'*X1val+clsq),'FaceColor','m','EdgeColor','w' ), hold on
bar( K+1:2*K, sign(wlsq'*X2val+clsq),'FaceColor','g','EdgeColor','w' ), hold on
title('LSQ validation')
axis square, hold off
axis([0 2*K -1 1])
set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultFigureColor', 'w');
addpath(genpath('utils'))
clear, close all
K = 50;
rng(1)

% example 1 (no covariance)
% mu1 = [5;0];
% mu2 = [15;5];
% sigma1 = [3,0;0,3];
% sigma2 = [3,0;0,3];
% example 2 (with covariance)
mu1 = [5;0];
mu2 = [5;5];
p = 0.9;
sigma1 = 10*[1,p;p,1];
sigma2 = 10*[1,p;p,1];
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
wheuristic = m2 - m1;
% global mean
m = 1/2*(m1 + m2);
% threshold
c = -w'*m;
cheuristic = (m1 - m2)'*m;
% generate 
s1 = linspace(min([X1(1,:),X2(1,:)]),max([X1(1,:),X2(1,:)]),100);
s2 = -1/w(2)*(c+w(1)*s1);
s2heuristic = -1/wheuristic(2)*(cheuristic+wheuristic(1)*s1);

%
lw = 5;
figure(1), clf
plot(X1(1,:),X1(2,:),'mo','MarkerFaceColor','m','LineWidth',lw), hold on
plot(X2(1,:),X2(2,:),'go','MarkerFaceColor','g','LineWidth',lw)
% a = wlsq/(wlsq'*wlsq)*wlsq'*X1;
% b = wlsq/(wlsq'*wlsq)*wlsq'*X2;
% scatter(a(1,:),a(2,:),'ro','filled'), hold on
% scatter(b(1,:),b(2,:),'bo','filled')
% plot(s1,s2,'--k','LineWidth',2)
% plot(s1,s2heuristic,'-k','LineWidth',2)
minAx1 = min([X1(1,:),X2(1,:)]);
minAx2 = min([X1(2,:),X2(2,:)]);
maxAx1 = max([X1(1,:),X2(1,:)]);
maxAx2 = max([X1(2,:),X2(2,:)]);
axis([minAx1, maxAx1, minAx2, maxAx2])
grid on, axis square
xlabel('feature 1'), ylabel('feature 2')
hold off

%%

figure(11), clf
plot(X1(1,:),X1(2,:),'mo','MarkerFaceColor','m','LineWidth',lw), hold on
plot(X2(1,:),X2(2,:),'go','MarkerFaceColor','g','LineWidth',lw)
a = wheuristic/(wheuristic'*wheuristic)*wheuristic'*X1;
b = wheuristic/(wheuristic'*wheuristic)*wheuristic'*X2;
scatter(a(1,:),a(2,:),'ro','filled'), hold on
scatter(b(1,:),b(2,:),'bo','filled')
% plot(s1,s2,'--k','LineWidth',2)
plot(s1,s2heuristic,'-k','LineWidth',2)
minAx1 = min([X1(1,:),X2(1,:)]);
minAx2 = min([X1(2,:),X2(2,:)]);
maxAx1 = max([X1(1,:),X2(1,:)]);
maxAx2 = max([X1(2,:),X2(2,:)]);
axis([minAx1, maxAx1, minAx2, maxAx2])
grid on, axis square
xlabel('feature 1'), ylabel('feature 2')
hold off

%%

figure(2), clf
plot(X1val(1,:),X1val(2,:),'md','MarkerFaceColor','m','LineWidth',lw), hold on
plot(X2val(1,:),X2val(2,:),'gd','MarkerFaceColor','g','LineWidth',lw)
plot(s1,s2,'--k','LineWidth',2)
plot(s1,s2heuristic,'-k','LineWidth',2)
axis([minAx1, maxAx1, minAx2, maxAx2])
xlabel('feature 1'), ylabel('feature 2')
grid on, axis square
hold off
legend('class 1','class 2','LDA','heuristic')
legend('location','SouthEast')

%%
figure(3)
subplot(1,2,1)
bar( 1:K, sign(w'*X1val+c),'FaceColor','m','EdgeColor','w' ), hold on
bar( K+1:2*K, sign(w'*X2val+c),'FaceColor','g','EdgeColor','w' ), hold on
title('LDA')
axis square, hold off
axis([0 2*K -1 1])
xticks([1 50 100]), yticks([-1 1])
ylabel('class')

subplot(1,2,2)
bar( 1:K, sign(wheuristic'*X1val+cheuristic),'FaceColor','m','EdgeColor','w' ), hold on
bar( K+1:2*K, sign(wheuristic'*X2val+cheuristic),'FaceColor','g','EdgeColor','w' ), hold on
title('heuristic')
axis square, hold off
axis([0 2*K -1 1])
xticks([1 50 100]), yticks([-1 1])
ylabel('class')
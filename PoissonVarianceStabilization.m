%% configuration
fsize = 32;
set(0, 'DefaultAxesFontSize', fsize);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');

clc
I = [5 15 25 35 45];
cmap = gray(length(I)) * 0.75;
cmap = ones(3,1)*linspace(0.2,0.8,length(I)); cmap = cmap';
% n = 0:3*max(I);
n = 0:max(I)+20;
%%
figure(1), clf
subplot(1,3,1)
cost = 0;
for k = 1:length(I)
    p1 = I(k).^n ./ gamma(n+1) .* exp(-I(k));

    hold on
    mu = I(k);
    s = sqrt(I(k));

    p2 = 1/sqrt(2*pi*s^2) * exp(-(n-mu).^2 / (2*s^2));
    plot(n,p2,'-','color',cmap(k,:),'lineWidth',2)
    axis square
    grid on

    cost = cost + norm(p1 - p2, 2);
end

for k = 1:length(I)
    p1 = I(k).^n ./ gamma(n+1) .* exp(-I(k));
    stem(n,p1,'Color',cmap(k,:),'MarkerFaceColor',cmap(k,:))

end

xlabel('$n$')
ylabel('$P(n | I)$')
h = legend('$5$', '$15$','$25$','$35$','$45$');
h.Location = 'NorthOutside';
h.Orientation ='horizontal';
disp(cost)

%
numSamples = 10000;
rng(0)
X = poissrnd(100, [numSamples, 1]);
Y = poissrnd(400, [numSamples, 1]);

% figure(2), clf
subplot(1,3,2)
histogram(X,15,'FaceColor',cmap(2,:),'FaceAlpha',0.25)
axis square
hold on
histogram(Y,15,'FaceColor',cmap(5,:),'FaceAlpha',0.25)
grid on
axis([0 500 0 2500])
h = legend('$X$','$Y$')
h.Location = 'NorthOutside';
h.Orientation ='horizontal';
xlabel('RV value')
xticks([10 20 30 40 50]*10)
ylabel('abs. frequency')

subplot(1,3,3)
histogram(2*sqrt(X),15,'FaceColor',cmap(2,:),'FaceAlpha',0.25)
axis square
hold on
histogram(2*sqrt(Y),15,'FaceColor',cmap(5,:),'FaceAlpha',0.25)
grid on
axis([0 50 0 2500])

xlabel('$2\sqrt{ RV value }$')
xticks([10 20 30 40 50])
ylabel('abs. frequency')
f1 = 2200*exp(-(n-2*sqrt(100)).^2/(2));
plot(n, f1,'k--','lineWidth',2,'color',[cmap(1,:),0.75])
f1 = 2200*exp(-(n-2*sqrt(400)).^2/(2));
plot(n, f1,'k-','lineWidth',2,'color',[cmap(3,:),0.75])

h = legend('$2\sqrt X$','$2\sqrt Y$');
h.Location = 'NorthOutside';
h.Orientation ='horizontal';

%%

var(2*sqrt(X(:)))
var(2*sqrt(Y(:)))
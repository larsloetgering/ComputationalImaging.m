%% configuration
fsize = 32;
set(0, 'DefaultAxesFontSize', fsize);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');

%
clc
I = [5 15 25 35 45];
cmap = gray(length(I)) * 0.75;
cmap = ones(3,1)*linspace(0.2,0.8,length(I)); cmap = cmap';
% n = 0:3*max(I);
n = 0:max(I)+20;

figure(1), clf
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
    stem(p1,'Color',cmap(k,:),'MarkerFaceColor',cmap(k,:))

end

xlabel('$n$')
ylabel('$P(n | I)$')
legend('$I=5$', '$I=15$','$I=25$','$I=35$','$I=45$')

disp(cost)
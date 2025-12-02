clear
%%
fsize = 32;
set(0, 'DefaultAxesFontSize', fsize);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');
%%
% physical parameters
L = 3;y1 = 1;y2 = 1;y3 = 1;
n1 = 1;n2 = 1.33;n3 = 1;
costFun = @(x1,x2) (n1*sqrt( x1.^2 + y1^2 ) + ...
                    n2*sqrt( x2.^2 + y2^2 ) + ...
                    n3*sqrt((L-x1-x2).^2 + y3^2));

x1 = linspace(0.5,1.5);
x2 = linspace(0.5,1.5);
im = costFun(x1,x2');

% loop
r = [1;1];
numIter = 5;
stepSize = 1/4;

%%
contourLevels = linspace(min(im(:)), max(im(:)), 50);
figure(1), clf
contourf(x1,x2,im, contourLevels)
axis image; cb = colorbar;
cmap = turbo(1000); 
cmap(1,:) = [0 0 0];
colormap(cmap);
cb.Label.Interpreter = 'latex';
cb.Label.String = 'cost $$\mathcal{L}(x_1,x_2)$$';
cb.Ticks = [4.68, 5.2]; % Specify the tick values
cb.TickLabels = {'Low', 'High'};
xlabel('$x_1$');ylabel('$x_2$')
hold on; plot(1,1,'wo','lineWidth',2,'MarkerSize',10)
cmap = flipud(gray(numIter+1))*0.8;
cmap = [cmap, 0.5*ones(size(cmap,1),1)]; 
cmap(1,:) = [0 0 0.8 1];
cmap(end,:) = [0 0.8 0 1];

[minCost,idx] = min(im(:));
[ii,jj] = ind2sub(size(im),idx);
figure(1)
hold on
plot(x2(jj), x1(ii) ,'xw','lineWidth',3,'MarkerSize',13)
hold off
% figure(2), clf
%%
for loop = 1:numIter
    %%
    figure(2)
    hold on
    plot([0,r(1),r(1)+r(2),L], [0,y1,y1+y2,y1+y2+y3],...
        'o-','color',cmap(loop,:),'lineWidth',2)
    hold off
    %%
    r_old = r;
    gradx = n1*r(1) / sqrt(r(1)^2+y1^2) - ...
       n3*(L-r(1)-r(2))/sqrt((L-r(1)-r(2))^2+y3^2);
    grady = n2*r(2) / sqrt(r(2)^2+y2^2) - ...
       n3*(L-r(1)-r(2))/sqrt((L-r(1)-r(2))^2+y3^2);
    r = r - stepSize * [gradx;grady];
    %%
    figure(1)
    hold on
    plot(r(1),r(2),'wo','lineWidth',2,'MarkerSize',7,'MarkerFaceColor','w')
    plot([r_old(1), r(1)],[r_old(2), r(2)],'w--','lineWidth',2,'MarkerSize',10)
    hold off
    %%
end
%%
figure(2)
hold on
plot([0,r(1),r(1)+r(2),L], [0,y1,y1+y2,y1+y2+y3],'o-','color',cmap(loop+1,:),'lineWidth',2)
hold off
hold on
yticks([0 1 2 3]);
plot([0 L],[1 1],'k-','lineWidth',2)
plot([0 L],[2 2],'k-','lineWidth',2)

fsize2 = fsize;
xlabel('$x$','FontSize',fsize2)
ylabel('$y$','FontSize',fsize2)
text(2.7,0.85,'$n_1$','FontSize',fsize2)
text(2.7,1.15,'$n_2$','FontSize',fsize2)
text(2.7,2.15,'$n_3$','FontSize',fsize2)
axis square, grid on
h = gca;
h.GridLineWidth = 2;
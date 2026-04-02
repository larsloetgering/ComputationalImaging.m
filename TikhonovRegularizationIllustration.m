%% configuration
clear
fsize = 32;
lw = 2;
set(0, 'DefaultAxesFontSize', fsize);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');

%%
rng(0)
x = -5:5;
xq = -5:0.1:5;
G = [randi(5,[1 5]),0:5];
G(1:2) = 0; G(end-1:end) = 0;
G = interp1(x,G,xq,'nearest');

G = G/max(G(:));

figure(1), clf
plot(xq,G,'k','linewidth',lw)
axis square
grid on
xticks(x)

h = exp(-xq.^2/(2*(2/2.355)^2));

hold on
plot(xq,h,'--','color',[1 1 1]*0.25,'linewidth',lw)

%% compute convolved signal (ignore noise here)

I = real( ifftshift(ifft(fft(h).*fft(G))) );
I = I / max(I(:));

figure(1)
plot(xq,I,'-xk','linewidth',lw)
legend('ground truth','conv. kernel','data')
h = gca;
h.GridLineWidth = 1.5;

%%

H = fft(h);
figure(2), clf
hold on
epsilon = [1e-7 1e-3];
marker = {'-','--'}
cmap = gray(size(epsilon,2))*0.25;

for loop = 1:length(epsilon)
    Gestimate = ifftshift(ifft( conj(H).*fft(I)./(abs(H).^2+epsilon(loop)) ));
    Gestimate = Gestimate / max(Gestimate(:));
    plot( xq,Gestimate, marker{loop}, 'color', cmap(loop,:),'linewidth',1.5 )
end
axis([-inf inf -0.1 1])
axis square
grid on
box on
xticks(x)
h = gca;
h.GridLineWidth = 1.5;
legend('$\epsilon = 10^{-7}$','$\epsilon = 10^{-3}$')






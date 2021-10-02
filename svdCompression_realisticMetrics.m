set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultFigureColor', 'w');
clear
im = single(imread('../testImages/Origami.jpg'));
im = im(:,:,2); % select green channel (out of rgb image)

[U,S,V] = svd(im);
N = min(size(im));
K = 180;

imCom = U(:,1:K)*S(1:K,1:K)*V(:,1:K)';
C = 1-( numel(U(:,1:K)) + numel(S(1:K,1:K)) + numel(V(:,1:K)))/numel(im);
D = 1/(1-C);

figure(1)
imagesc(imCom)
axis image off 
colormap gray
title({['K = ',num2str(K)],...
    ['C = ', num2str(round(1000*(C))/10),'% (',num2str(round(D)),'x)']})


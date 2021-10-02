set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultFigureColor', 'w');
clear
im = single(imread('../testImages/Origami.jpg'));
im = im(:,:,2); % select green channel (out of rgb image)

[U,S,V] = svd(im);
N = min(size(im));
K = 18;

imCom = U(:,1:K)*S(1:K,1:K)*V(:,1:K)';
    
figure(1)
imagesc(imCom)
axis image off 
colormap gray
title(['C = ', num2str(round(1000*(1-K/N))/10),'% (',num2str(round(N/K)),'x)'])


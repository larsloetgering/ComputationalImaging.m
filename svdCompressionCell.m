set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultFigureColor', 'w');
clear
im = single(imread('../testImages/Origami.jpg'));
im = im(:,:,2); % select green channel (out of rgb image)

tic
% choose c = (2,3,4,5,6<slow!!)8,10,12,15,16,20,24,30
c = 30;
im2 = reshape_cXc(im, c);
[U,S,V] = svd(im2);
M = size(im2,1);
N = size(im2,2);
K = 10; % choose K <= c
imCom = U(:,1:K)*S(1:K,1:K)*V(:,1:K)';
C = 1-K*(M+K+N)/(M*N);
D = 1/(1-C);
% compressed image
imCom = reshape_cXc_inv(imCom,size(im,1),size(im,2),c);

figure(1)
imagesc(imCom)
axis image off
colormap gray
title({['c = ', num2str(c),', K = ',num2str(K)],...
    ['C = ', num2str(round(1000*(C))/10),'% (',num2str(round(D)),'x)']})
toc

function B = reshape_cXc(A,c)
B = zeros(c^2,size(A,1)*size(A,2)/c^2);
counter = 0;
for k = 1:c:size(A,1)
    for l = 1:c:size(A,2)
        counter = counter+1;
        B(:,counter) = reshape(A(k:k+c-1,l:l+c-1),[c^2,1]);
    end
end
end

function A = reshape_cXc_inv(B,M,N,c)
A = zeros(M,N);
counter = 0;
for k = 1:c:M
    for l = 1:c:N
        counter = counter+1;
        A(k:k+c-1,l:l+c-1) = reshape(B(:,counter),[c,c]);
    end
end
end
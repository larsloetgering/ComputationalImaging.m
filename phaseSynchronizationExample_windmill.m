set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultFigureColor', 'w');
clear, close all
addpath(genpath('utils'))

rgbImage = imread('../testImages/Windmill.png');

im_red = single(rgbImage(:,:,1)); % select red channel (out of rgb image)
im_red = im_red / max(im_red(:));
im_red = im_red .* exp(-1i * pi * im_red);

im_green = single(rgbImage(:,:,2)); % select green channel (out of rgb image)
im_green = im_green / max(im_green(:));
im_green = im_green .* exp(-1i * pi * im_green);

exp_iPhi = phaseSynchronization(im_red, im_green);

figure(1)
subplot(2,2,1)
imshow(rgbImage)
title('original rgb image')
axis off

subplot(2,2,3)
hsvplot(im_red)
title('red channel')
axis off

subplot(2,2,4)
hsvplot(im_green)
title('green channel')
axis off

subplot(2,2,2)
hsvcolorbarplot(exp_iPhi*im_red,'colorbar','test')
title({'red channel','synchronized'})
axis off

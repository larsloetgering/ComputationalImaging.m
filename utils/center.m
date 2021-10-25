function [ys, xs]  = center( I )
% center image
% [r,ys,xs]  = center( I )
% input: 
% I: gray image
% output:
% ys: center row
% xs: center column
[M, N] = size(I);
y = linspace(-M/2, M/2, M);
x = linspace(-N/2, N/2, N);
[X,Y] = meshgrid(x, y);
% normalize
I = I / sum2(I);
xs = sum2(I .* X);
ys = sum2(I .* Y);
return
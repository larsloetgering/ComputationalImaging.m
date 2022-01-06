N = size(Z,3);
C = zeros(N,N);
for n=1:N
    for m=1:N
        C(m,n)=sum(sum(Z(:,:,n).*Z(:,:,m)));
    end
end

N = 2^8;
x = linspace(-1,1,N);
dx = x(2) - x(1);
[X,Y] = meshgrid(x);
[THETA, R] = cart2pol(X,Y);

[Z1, J] = zernike(R,THETA,2,4); % 0,2 vs 0,0
[Z2, J] = zernike(R,THETA,2,2);
sum(sum(Z1.*Z1)) 
sum(sum(Z2.*Z2)) 
sum(sum(Z1.*Z2)) / sum(sum(Z2.*Z2)) % 0.5 percent error

% how to improve:
% Gram-Schmidt Orthogonalization
% Calculate on a finer grid and bin:

N = 2^10;
x = linspace(-1,1,N);
dx = x(2) - x(1);
[X,Y] = meshgrid(x);
[THETA, R] = cart2pol(X,Y);

[Z1, J] = zernike(R,THETA,2,4); % nicely sampled
[Z2, J] = zernike(R,THETA,2,2);
Z1 = squeeze(sum(sum(reshape(Z1,[4,256,4,256]),1),3)); % bin 4x4
Z2 = squeeze(sum(sum(reshape(Z2,[4,256,4,256]),1),3));
Z1 = Z1 ./ norm(Z1(:));
Z2 = Z2 ./ norm(Z2(:));

        figure%(1)
        imagesc(x,x,Z2(:,:,end))  % scaled image display
        xticks([-1 0 1]), yticks([-1 0 1])
        axis image off
        colormap(cmap)        
        title(['n = ',num2str(n), ', m = ',num2str(m), ', J = ',num2str(J)])

sum(sum(Z1.*Z1)) 
sum(sum(Z2.*Z2)) 
sum(sum(Z1.*Z2)) % 0.5 percent error

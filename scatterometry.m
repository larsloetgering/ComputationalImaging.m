% scatterometry

addpath(genpath('utils'))
% clc, clear
set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultFigureColor', 'w');
rng(1)
% photons per diffraction pattern
photonBudget = 1e3;
% number of diffraction patterns
K = 50;
% generate meshgrid
N = 64;
x = linspace(-N/2,N/2,N);
[X,Y] = meshgrid(x);
width = 5;

figure(1)
subplot(1,2,1)
xi = 1+rand(1,1);
hsvplot(exp(1i * pi * exp(- (X.^2 + Y.^2)/(2*width^2))/xi ))

subplot(1,2,2)
xi = 2+2*rand(1,1);
eta = 2+2*rand(1,1);
theta = randi(180);
hsvplot(exp(1i * pi * ...
        exp(- ((cos(theta)*X+sin(theta)*Y).^2/xi + ...
               (-sin(theta)*X+cos(theta)*Y).^2*xi)/(2*width^2))))
           
%% generate K data points, each for training and validation

data = zeros(N,N,K);
valdata = zeros(N,N,K);
% class 1
for k = 1:K/2
    % training data
    xi = 1+rand(1,1);
    data(:,:,k) = fft2c( exp(1i * pi * exp(- (X.^2 + Y.^2)/(2*width^2))/xi ) );
    % validation data
    xi = 1+rand(1,1);
    valdata(:,:,k) = fft2c( exp(1i * pi * exp(- (X.^2 + Y.^2)/(2*width^2))/xi ) );
end
% class 2
for k = K/2+1:K
    % training data
    xi = 2+2*rand(1,1);
    theta = randi(180);
    data(:,:,k) = fft2c(...
                  exp(1i * pi * ...
                  exp(-((cos(theta)*X+sin(theta)*Y).^2/xi + ...
                       (-sin(theta)*X+cos(theta)*Y).^2*xi)/(2*width^2))));
    % validation data
    xi = 2+2*rand(1,1);
    theta = randi(180);
    valdata(:,:,k) = fft2c(...
                  exp(1i * pi * ...
                  exp(-((cos(theta)*X+sin(theta)*Y).^2/xi + ...
                       (-sin(theta)*X+cos(theta)*Y).^2*xi)/(2*width^2))));
end

% convert to intensity
data = abs(data).^2;
valdata = abs(valdata).^2;
% add noise
data = data / sum(data(:))*photonBudget*K;
valdata = valdata / sum(valdata(:))*photonBudget*K;
data(:) = poissrnd(data(:));
valdata(:) = poissrnd(valdata(:));

%% inspect data

figure(3)
subplot(1,2,1)
imagesc(log10([data(:,:,1),data(:,:,2);...
               data(:,:,3),data(:,:,4)]+1e-3))
axis image off, colormap hot

subplot(1,2,2)
imagesc(log10([data(:,:,K/2+1),data(:,:,K/2+2);
               data(:,:,K/2+3),data(:,:,K/2+4)]+1e-3))
axis image off, colormap hot

%% extract features

x1 = zeros(3, K/2);
x2 = zeros(3, K/2);
valx1 = zeros(3, K/2);
valx2 = zeros(3, K/2);
r = sqrt(X.^2 + Y.^2);
phi=atan2(Y,X);

% moment kernels
kern1 = exp(1i*2*phi);
kern2 = 1;
kern3 = r;
% calculate moments (class 1)
% training
x1(1,:) = 1/N^2 * squeeze( real(sum(sum( kern1.*data(:,:,1:K/2) ))));
x1(2,:) = 1/N^2 * squeeze(sum(sum( kern2.*data(:,:,1:K/2) )));
x1(3,:) = sqrt(1/N^2 * squeeze(sum(sum( kern3.*data(:,:,1:K/2) ))));
% validation
valx1(1,:) = 1/N^2 * squeeze(real(sum(sum( kern1.*valdata(:,:,1:K/2) ))));
valx1(2,:) = 1/N^2 * squeeze(sum(sum( kern2.*valdata(:,:,1:K/2) )));
valx1(3,:) = sqrt(1/N^2 * squeeze(sum(sum( kern3.*valdata(:,:,1:K/2) ))));

% calculate moments (class 2)
% training
x2(1,:) = 1/N^2 * squeeze(real(sum(sum( kern1.*data(:,:,K/2+1:end) ))));
x2(2,:) = 1/N^2 * squeeze(sum(sum( kern2.*data(:,:,K/2+1:end) )));
x2(3,:) = sqrt(1/N^2 * squeeze(sum(sum( kern3.*data(:,:,K/2+1:end) ))));
% validation
valx2(1,:) = 1/N^2 * squeeze(real(sum(sum( kern1.*valdata(:,:,K/2+1:end) ))));
valx2(2,:) = 1/N^2 * squeeze(sum(sum( kern2.*valdata(:,:,K/2+1:end) )));
valx2(3,:) = sqrt(1/N^2 * squeeze(sum(sum( kern3.*valdata(:,:,K/2+1:end) ))));

figure(4)
clf
subplot(1,2,1)
plot3(x1(1,:),x1(2,:),x1(3,:),'mo','MarkerFaceColor','m'), hold on
plot3(x2(1,:),x2(2,:),x2(3,:),'bo','MarkerFaceColor','b')

m1 = mean(x1,2);
m2 = mean(x2,2);
% between class separation matrix
Sb = (m2 - m1)*(m2 - m1)';
% within class covariance
Sw = (x1 - m1)*(x1 - m1)' + ...
     (x2 - m2)*(x2 - m2)';
% separation surface normal
w = Sw\(m2 - m1);
% w = (m2 - m1);
% global mean
m = 1/2*(m1 + m2);
% threshold
c = -w'*m;

% plot separation plane (with random samples)
P = 20;
[s1,s2] = meshgrid(linspace( min([x1(1,:),x2(1,:)]),max([x1(1,:),x2(1,:)]), P ),...
                   linspace( min([x1(2,:),x2(2,:)]),max([x1(2,:),x2(2,:)]), P ));
s1 = s1(:); s2 = s2(:);
s3 = -1/w(3)*(c+w(1)*s1+w(2)*s2);

figure(4)
subplot(1,2,1)
plot3(s1,s2,s3,'yo')
% legend('round','elliptic','boundary')
hold off
xlabel('feature 1')
ylabel('feature 2')
zlabel('feature 3')
grid on
axis square
% legend('round','elliptic','boundary')
view(3,5)

subplot(1,2,2)
plot3(x1(1,:),x1(2,:),x1(3,:),'mo','MarkerFaceColor','m'), hold on
plot3(x2(1,:),x2(2,:),x2(3,:),'bo','MarkerFaceColor','b')
plot3(s1,s2,s3,'yo')
hold off
grid on
axis square
xlabel('feature 1')
ylabel('feature 2')
zlabel('feature 3')
% legend('round','elliptic','boundary')
view(95,40)

%%  show validation results
figure(5)
plot3(valx1(1,:),valx1(2,:),valx1(3,:),'rx','LineWidth',2),hold on
plot3(valx2(1,:),valx2(2,:),valx2(3,:),'gx','LineWidth',2)
plot3(s1,s2,s3,'yo'), hold off
legend('round','elliptic','boundary')
grid on
axis square
xlabel('feature 1')
ylabel('feature 2')
zlabel('feature 3')

figure(6)
subplot(1,2,1)
bar( sign(w'*valx1+c),'FaceColor','k','EdgeColor','w' )
axis([0 K/2+1 -1 1])
title('round validation')
axis square
subplot(1,2,2)
bar( sign(w'*valx2+c),'FaceColor','k','EdgeColor','w' )
axis([0 K/2+1 -1 1])
axis square
title('elliptic validation')


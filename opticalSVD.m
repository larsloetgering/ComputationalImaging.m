% configuration
addpath(genpath('utils'))
set(0, 'DefaultAxesFontSize', 24);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
%% sampling grid and physical parameters
% physical quantities
wavelength = 632.8e-9;  % wavelength
% detector coodinates
Nd = 2^7;               % number of samples in detector plane
dxd = 5e-6;             % detector pixel size
Ld = Nd * dxd;          % detector size
xd = linspace(-Nd/2,Nd/2,Nd)*dxd;% 1D coordinates (detector)
[Xd, Yd] = meshgrid(xd);% 2D coordinates (detector)
z = 10e-2;              % source-receiver distance
% source coordinates
dxs = wavelength*z/Ld/4;  % source sampling
Ns = Nd;        % number of samples in source plane
Ls = Ns * dxs;  % source field of view
xs = linspace(-Ns/2,Ns/2,Ns)*dxs;% 1D coordinates (source)
[Xs, Ys] = meshgrid(xs);% 2D coordinates (source)

%% get source - object eigenmodes
% generate pinhole
diameter = Ls/3;
source = circ(Xs, Ys, diameter);
source = normconv2(source, [1,1;1,1]);
figure(1); imagesc(source); colormap gray 
axis image; title('source') % show result
% find non-zero source indices
[row, col] = find( source > 0 );  
% number of source points 
numSourcePoints = size(row,1);           
% (mutually uncorrelated) spherical wavelets
sphericalWavelets = zeros(Ns, Ns, numSourcePoints,'single'); 
% determines how often source wavelets are shown
figureUpdateFrequency = 50;         
for loop = 1:numSourcePoints
   % displacement
   R = sqrt( (Xd - Xs( row(loop), col(loop) )).^2 + ...
             (Yd - Ys( row(loop), col(loop) )).^2 + z.^2 );
   % Greens function (Rayleigh-Sommerfeld) for each point source
   sphericalWavelets(:,:,loop) = source(row(loop),col(loop))*...
                  exp(1i * 2*pi / wavelength .* R) .* z ./ R.^2;
    
   if mod(loop,figureUpdateFrequency)==0
       figure(2)
       hsvplot(sphericalWavelets(:,:,loop)); 
       axis off; 
       title({'spherical wavelets in entrance pupil',...
               num2str(loop)}); 
       drawnow
   end
end 

%% receiver aperture
receiverAperture = normconv2(circ(Xd, Yd, Ld/2), [1,1;1,1]);
sphericalWavelets = bsxfun(@times, sphericalWavelets, receiverAperture);

%% get orthogonal modes
maxNumModes = 500; % maximum number of modes to be calculated
wrank = min(numSourcePoints, maxNumModes); % # modes for tsvd 
[U,S,V] = tsvd(reshape(sphericalWavelets,...
                    [Nd^2, numSourcePoints]), wrank); 
s = diag(S);
purity = sqrt(sum(s.^2))/sum(s)^2;

%% orthogonal source and receiver modes 
w = zeros(Ns, Ns, wrank, 'single');
u = zeros(Nd, Nd, wrank, 'single');
idx = find( source > 0 );
for k = 1:wrank
    u(:,:,k) = reshape( U(:, k), Nd * [1,1] );
    temp = 0*source;
    temp(idx) = V(:,k);
    w(:,:,k) = temp;

end

%% show source and receiver modes
n = round(diameter/dxs/2) + 10;
figure(3)
subplot(1,3,1)
hsvmodeplot(w(end/2-n:end/2+n,end/2-n:end/2+n,1:min(wrank, 25)))
title('source modes'), axis off

subplot(1,3,2)
hsvmodeplot(u(:,:,1:min(wrank, 25)))
title('receiver modes'), axis off

subplot(1,3,3)
semilogy(diag(S(1:100,1:100))/trace(S),'ko','MarkerFaceColor','k')
axis square, grid on
title('singular values')
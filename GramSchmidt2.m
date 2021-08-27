clear
clc
set(0, 'DefaultAxesFontSize', 30);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
addpath(genpath('/home/user/Dropbox/Codes/fracPty/utils/export_fig2018'))
% matlab code: Gram-Schmidt

rng(666 , 'twister')
n = 10;
% B = rand(n);
B = hilb(n);
precision = 'double'; % or 'single'
Qc = zeros(size(B), precision);

% classical Gram-Schmidt (CGS)
for j = 1:size(B,2)
    q = B(:,j);
    % note: the following loop will be executed only for j>1
    for k = 1:j-1
        q = q - (Qc(:,k)'*q)*Qc(:,k);
    end
        
    Qc(:,j) = q / sqrt(q'*q);
end
% errorCGS = norm(eye(n)-Qc*Qc',2);
errorCGS = abs(abs(det(Qc))-1);
disp(['error CGS: ', num2str(errorCGS)])

% modified Gram-Schmidt (MGS)
Qm = B;
for j = 1:size(B,2)
    Qm(:,j) = Qm(:,j)/sqrt(Qm(:,j)'*Qm(:,j));
    for k = (j+1):n
        Qm(:,k) = Qm(:,k) - (Qm(:,j)'*Qm(:,k))*Qm(:,j);
    end
end
% errorMGS = norm(eye(n)-Qm*Qm',2);
errorMGS = abs(abs(det(Qm))-1);
disp(['error MGS: ', num2str(errorMGS)])

figure(1)
imagesc(Qc)
axis image
yticklabels('')
% title('Q_{CGS,d}')
colorbar 
colormap hot
export_fig(['Qc',precision,'.png'])

figure(2)
% subplot(2,3,2)
imagesc(Qm)
axis image
yticklabels('')
% title('Q_{MGS,d}')
colorbar
colormap hot
export_fig(['Qm',precision,'.png'])

% colorbar('location','NorthOutside')

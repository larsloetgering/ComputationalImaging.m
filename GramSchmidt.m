clear
clc
set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultFigureColor', 'w');
% matlab code: Gram-Schmidt

n = 3;
% rng(666 , 'twister'), B = rand(n); % alternative B
B = hilb(n);
precision = 'double'; % also try 'single'
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
errorMGS = abs(abs(det(Qm))-1);
disp(['error MGS: ', num2str(errorMGS)])

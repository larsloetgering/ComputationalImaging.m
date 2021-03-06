clc
% index "s": small
% index "l": large
tic
N = 4;
w = exp(1i * 2*pi / N);
Isub = eye(N/2);
Fsub = fft(Isub);
dsub = w.^(0:N/2-1)';
Dsub = diag( dsub );
% matrix-based implementation (slow)
% compute ground truth
Flarge = fft( eye(N) );
% permutation matrix
Plarge = zeros(N);
idx = sub2ind([N N], 1:N, [1:2:(N-1), 2:2:N]);
Plarge( idx ) = 1;
% compute fft step
Glarge = [Isub, Dsub; Isub, -Dsub]*...
         [Fsub, zeros(N/2); zeros(N/2), Fsub]*Plarge;
toc
% check Frobenius norm 
disp(norm(Flarge - Glarge, 'fro'))
% more efficient implementation 
% > FFT recursion directly applied to random vector v
rng(0)
v = rand(N,1);
% notice that index 1 in Matlab corresponds 
% to 0 when using zero-based idexing
tic
veven = v(1:2:N-1);  % even
vodd = v(2:2:N);     % odd
veven_ft = fft(veven);
vodd_ft = fft(vodd);
v_ft = [veven_ft + dsub.*vodd_ft; veven_ft - dsub.*vodd_ft];
toc
% check L2 norm 
disp(norm(v_ft - fft(v), 2))
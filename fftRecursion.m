clc
% index "s": small
% index "l": large
tic
N = 8;
w = exp(1i * 2*pi / N);
Isub = eye(N/2);
Fsub = fft(Isub);
dsub = w.^(0:N/2-1)';
Dsub = diag( dsub );
% matrix-based implementation (slow)
% compute ground truth
Fl = fft( eye(N) );
% compute fft step
Pl = zeros(N);
idx = sub2ind([N N], 1:N, [1:2:(N-1), 2:2:N]);
Pl( idx ) = 1;
Gl = [Isub, Dsub; Isub, -Dsub]*...
     [Fsub, zeros(N/2); zeros(N/2), Fsub]*Pl;
toc
% check Frobenius norm 
disp(norm(Fl - Gl, 'fro'))
% more efficient implementation 
% > FFT recursion directly applied to random vector v
rng(0)
v = rand(N,1);
% notice that index 1 in Matlab corresponds 
% to 0 when using zero-based idexing
tic
veven = v(1:2:N-1); 
vodd = v(2:2:N);
veven_ft = fft(veven);
vodd_ft = fft(vodd);
v_ft = [veven_ft + dsub.*vodd_ft; veven_ft - dsub.*vodd_ft];
toc
% check L2 norm 
disp(norm(v_ft - fft(v), 2))
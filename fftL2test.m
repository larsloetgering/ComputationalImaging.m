clc % clear command window
% generate signal
x = (1:3)';
% perform (normalized) fft
y = fftL2(x);
% perform inverse (normalized) fft
z = ifftL2(y);
% display L2 squared
disp(x'*x)
disp(y'*y)
disp(z'*z)

function v = fftL2(u)
% assumes input is 1D vector
v = fft(u)/sqrt( numel(u) );
end

function u = ifftL2(v)
% assumes input is 1D vector
u = ifft(v)*sqrt( numel(v) );
end
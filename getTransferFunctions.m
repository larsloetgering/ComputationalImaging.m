function [PTF, ATF] = getTransferFunctions(source, pupil)
% ATF: absorption transfer function
% PTF: phase transfer function
% assumes object exp(mu-j*phi) ~ 1 + mu - j*phi

term1 = ifft2( fft2( source .* conj(pupil) ) .* ...
         conj( fft2( conj(pupil) ) ) );
term1 = rot90( term1, 2 );
term2 = ifft2( fft2( source .* pupil ) .* ...
         conj( fft2( pupil ) ) );

term1 = fftshift(term1);
term2 = fftshift(term2);

normalization = sum( source .* abs(pupil).^2, [1,2] );
ATF = (term1 + term2)/normalization;
PTF = 1j*(term1 - term2)/normalization;
end
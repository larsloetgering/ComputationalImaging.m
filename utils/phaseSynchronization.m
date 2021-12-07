function exp_iPhi = phaseSynchronization(A, B)
% A, B: complex-valued images of same dimension
% exp_iPhi: minimizes the cost 
% L = norm( exp_iPhi * A(:) - B(:), 2 )
a = A(:); % flatten image A to vector
b = B(:); % flatten image A to vector
% compute phase synch. term
exp_iPhi = sqrt( (a'*b) / (b'*a) );
% compare L2 norms of candidate solutions 
L2_1 = norm(  exp_iPhi * a - b, 2 );
L2_2 = norm( -exp_iPhi * a - b, 2 );
if L2_2 < L2_1
    exp_iPhi = -exp_iPhi;
end

end
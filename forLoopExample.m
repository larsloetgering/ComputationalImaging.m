% maximum power
N = 10;
% base
q = 0.5;
% sum
s = 0;
for k = 0:N
    s = s + q^k;
end

disp(['for-loop result: ', num2str(s)])
disp(['analytic result: ', num2str( (1-q^(N+1))/(1-q) )])
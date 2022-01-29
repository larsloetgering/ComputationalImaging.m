%% Collatz
tic
set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
%%
% collatzFun(670617279)
%%
inputRange = 16;
computationTime = zeros(1,inputRange);

for k = 1:inputRange

    count = 0;
    n = k;
    while n ~= 1
        
        count = count + 1;
        n = collatzFun(n);
    end
    computationTime(k) = count;
end

cmap = 0.9*jet(inputRange);
[~,idx] = sort(computationTime,'ascend');


figure(1), clf
% scatter(1:inputRange, computationTime,[], cmap,'filled')
scatter(idx, computationTime(idx),[], cmap,'filled')
axis square
xlabel('input n')
ylabel('Collatz iteration counts')
grid on

powersOfTwo = 2.^(0:floor(log2(inputRange)));
hold on
plot(powersOfTwo, computationTime(powersOfTwo),'k-')
hold off
toc
% h = gca;
% h.XScale = 'log';

function n = collatzFun(n)
% function to execute a single collatz iteration
if mod(n,2) ==  0
    n = n/2;
else
    n = 3*n + 1;
end
end

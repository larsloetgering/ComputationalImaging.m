function [S, index_list] = NFINDR(D, r)
% input: data matrix D (two-dimensional)
% output: spectra S
% r: rank (or number of spectral components)
% assumptions: 
% D >=0, rows = pixel, channels = columns
% see https://doi.org/10.1117/12.366289

num_px = size(D,1);
% normalize by L1 norm if needed
if abs( sum(D(:)) - num_px ) > 1e-7
    D = D./(sum(D, 2) + eps);
end
% define volume function (constant factors ignored)
getVol = @(X)(det(X*X'));
% note: getVol is uses XX' instead of X'X, because 
% below E is tranposed as compared to original paper

% initialize
index_list = randperm( num_px, r );
E = [D(index_list,:), ones(r,1,'like',D)];
maxVolume = getVol(E);

% optimization loop
for index_loop = 1:r
    test_index_list = index_list;
    for px_loop = 1:num_px
        test_index_list(index_loop) = px_loop;
        E = [D(test_index_list,:), ones(r,1)];
        testVolume = getVol(E);
        if testVolume > maxVolume
            maxVolume = testVolume;
            index_list = test_index_list;
        end
    end
end
% output
S = D(index_list,:);
end
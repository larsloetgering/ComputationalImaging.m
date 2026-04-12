function r = var_stab(x)
r = x - mean(x, [1,2]);
r = r ./ std(r, [], [1,2]);
end
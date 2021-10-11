function r = circ(X, Y, D)
% r = circ(X, Y, D)
% D: diameter
% X, Y: meshgrid (in x/y-direction)

r = single((X.^2+Y.^2) <= D^2/4);

end
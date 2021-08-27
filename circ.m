function z = circ(x, y, D)
% function z = circ(x, y, D)

r = sqrt(x.^2+y.^2);
z = single(r <= D/2);

end
A1 = [1 2 7; 1 4 5; 3 2 5];
A2 = [1 1 -2; 1 -2 1; -2 1 1];
% un/comment (ctrl + r/t) each of the following examples
% example 1:
% A = A1;
% example 2:
% A = A2;
% example 3:
A = A2 + eye(size(A2));
v = rand(3,1,'single');

numIter = 100;

% power method
for k = 1:numIter
    v = A*v;
    v = v / sqrt(v'*v);
end
% eig solution
[X,D] = eig(A);

% show results
disp(v)
disp(diag(D))
disp(X)
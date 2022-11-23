set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultFigureColor', 'w');
A1 = [1 2 7;... 
      1 4 5;... 
      3 2 5];
A2 = [1 1 -2;... 
      1 -2 1;... 
      -2 1 1];
A3 = A2 + eye(size(A2));
A4 = [1 2 3; 0 1 2; 0 0 1];
% un/comment (ctrl + r/t) each of the following examples
% example 1:
% A = A1;
% example 2:
% A = A2;
% example 3:
A = A3;
% example 4:
% A = A4; 
% initialize start vector
v = rand(3,1,'single');

numIter = 10;

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
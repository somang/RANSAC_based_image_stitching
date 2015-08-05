%% 1
% z = alpha*x + beta*y + gamma in matrix form.
% [Zi] = [Xi Yi 1]*[alpha beta gamma]'

%% 2
% N can be 3 at its minimum, for (xi,yi,zi) where i = 1,...,N 
% Since there are Three unknowns and Three equations.

%% 3
x = [1:500]';
a = 5;
b = 5;
c = 5;

y = [1:500]';
y = y + randn(size(y));

z = a*x + b*y + c;
z = z + randn(size(z));

%% 4
A = [x y ones(size(z))]; % Ax = b
estimates = A\z % x = A\b

%% 5
absolute_error = abs([a - estimates(1); b - estimates(2); c - estimates(3)])
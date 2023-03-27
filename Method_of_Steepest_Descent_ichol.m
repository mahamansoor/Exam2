function [ x, niters ] = Method_of_Steepest_Descent_ichol( A, b, x0 )
% Define the inputs: A, b, x0, and M
% A is the system matrix
% b is the right-hand side vector
% x0 is the initial guess for the solution
% M is the incomplete Cholesky factorization of A
% L is the lower-triangular matrix obtained from M
% The parameter maxit is the maximum number of iterations

function [x, niters] = Method_of_Steepest_Descent_ichol(A, b, x0)

% Compute the incomplete Cholesky factorization of A and obtain L
L = sparse(ichol(sparse(A)));

% Define initial residual r0 and set k to 0
r = b - A*x0;
k = 2000;

% Start loop
while norm(r) > 1e-7
% Calculate p_k and q_k
p = L' \ (L \ r);
q = A*p;

% Calculate alpha_k and update x_k
alpha = (p'*r) / (p'*q);
x = x0 + alpha*p;
end
% Calculate new residual r and update k and x0
r = r - alpha*q;
k = k + 1;
x0 = x;
end

% Output the number of iterations taken
niters = k;

end

function [x, niters] = PCG(A, b, x0, tol)
% PCG: Preconditioned Conjugate Gradient Method with Incomplete
% Cholesky Factorization
%
% Input:
% A: n x n symmetric positive definite matrix
% b: n x 1 right-hand side vector
% x0: n x 1 initial guess
% tol: convergence tolerance
%
% Output:
% x: n x 1 solution vector
% niters: number of iterations
% Incomplete Cholesky factorization
L = sparse(ichol(sparse(A), struct('type', 'ict', 'droptol', 1e-2)));
% Initialize variables
x = x0;
r = b - Ax;
z = L'(L\r);
p = z;
niters = 0;
while norm(r) > tol
niters = niters + 1;
q = Ap;
alpha = (p'z)/(p'q);
x = x + alphap;
r_new = r - alphaq;
z_new = L'(L\r_new);
beta = (z_new'*q)/(z'q);
p = z_new + betap;
r = r_new;
z = z_new;
end
end

function [x, niters] = PCG(A, b, L)

% Initialize variables
x = zeros(size(b));
r = b - A*x;
z = L'\(L\r);
p = z;
q = A*p;
niters = 0;

% Loop until convergence or maximum iterations
while norm(r) > 0 && niters < length(b)
    niters = niters + 1;
    alpha = (p'*z)/(p'*q);
    x = x + alpha*p;
    r_new = r - alpha*q;
    z_new = L'\(L\r_new);
    beta = (z_new'*r_new)/(z'*r);
    p = z_new + beta*p;
    q = A*p;
    r = r_new;
    z = z_new;
end

end
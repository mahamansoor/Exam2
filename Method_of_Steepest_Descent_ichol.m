function [ x, niters ] = Method_of_Steepest_Descent_ichol( A, b, x0 )
% We cab define the inputs: A, b, x0, and M
% where A is the system matrix
% b is the right-hand side vector
% x0 is the initial variable for the solution
% M is the incomplete Cholesky factorization of A
% L is the lower-triangular matrix obtained from M
% The parameter maxit is the maximum number of iterations
% Compute the incomplete Cholesky factorization of A and obtain L
% the below function performs the incomplete cholesky factorization of A. ichol references the lower triangule of A 
% and produces lower traingular factors. 
L = sparse(ichol(sparse(A)));
% Define initial residual r0 and set k to 0. here we are caluclating our residual value by substratcing matrix A*x0 by b value 
r = b - A*x0;
k = 2000;
% number of iterations are set to 2000 for our steepest gradient descent 
% This step initiates the Start loop process. 
while norm(r) > 1e-7 % tolerance is set to 1e-7
% Calculate p_k and q_k 
p = L' \ (L \ r);
% we can multiply A * p and set it equal to q 
q = A*p;
% we can then calculate alpha_k and update x_k
alpha = (p'*r) / (p'*q);
%In this step we can set up x value equal to x0 + alpha *p where alpha value is used to choose the step with the deepest descent, 
x = x0 + alpha*p;
end % here we have ended the loop and after this step we will solve the residuals for our convergence 
% Calculate new residual r and update k and x0
r = r - alpha*q;
k = k + 1; % on this step, we set niters = niters + 1 or k = k+ 1
x0 = x;
end
% Output the number of iterations taken
niters = k; % this step will give us our number of iterations
end % at this end we have ended the matlab code for steepest descent algorithm 

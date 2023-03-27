function [x, niters] = PCG(A, b, x0, tol)
% PCG: Preconditioned Conjugate Gradient Method with Incomplete
% Cholesky Factorization
% in this algorithm we are initializing the following variables 
% Input:
% A: n x n symmetric positive definite matrix
% b: n x 1 right-hand side vector
% x0: n x 1 initial variable 
% tol: convergence tolerance
% to get output, we are using the following equation 
% x: n x 1 solution vector
% niters: number of iterations
% with the help of Incomplete Cholesky factorization
L = sparse(ichol(sparse(A), struct('type', 'ict', 'droptol', 1e-2)));
% Initialize variables as follows: - 
x = x0; % we can set x as x0 
r = b - Ax; % for residual we can substratc b from Ax
z = L'(L\r); % then z is our preconditioned residual for conjugate gradient method 
p = z; % we will initialize our pre-conditioned residual for p 
niters = 500; % we will number of iterations as 500 
while norm(r) >  1e-6 %  % this step is used to set norm of vector r less than tolerance value. our tolerance value is 1e-6
niters = niters + 1; % this step will add 1 every time to the iteration number 
q = Ap; % here we are multiplying A by p 
alpha = (p'z)/(p'q); % we have tried to determine the alpha value by dividing the variables for preconditioned residual by conditioned variable 
x = x + alphap; % we will set our x value to x + aplha p values 
r_new = r - alphaq; % we will then our new residual value to r - alpha q
z_new = L'(L\r_new); % we will then initialzed a new z value for preconditioning 
beta = (z_new'*q)/(z'q); % we will then calulate our beta variable similar to like we did for alpha value 
p = z_new + betap; % we will then initialize our p value where we can add z new variable and beta p 
r = r_new; % we will then set our r value to residual new value 
z = z_new; % same for z value 
end
end
% the matlab code below can be used to output niters function 
function [x, niters] = PCG(A, b, L)
% Initialize variables
x = zeros(size(b)); % we initialize x, residual r, z as preconditioning variable, p for preconditioning and q as multiplier for A by p 
r = b - A*x;
z = L'\(L\r);
p = z;
q = A*p;
niters = 500; % here we set our iteration value to 500
% We will then write the matlab code to check for convergence and maximum iterations
% Loop until convergence or maximum iterations
while norm(r) > 0 && niters < length(b) % here we will make the while loop 
    niters = niters + 1;
    alpha = (p'*z)/(p'*q);
    x = x + alpha*p;
    r_new = r - alpha*q; % we will initialize r_new value and z_new values 
    z_new = L'\(L\r_new); 
    beta = (z_new'*r_new)/(z'*r); % we will then set beta values 
    p = z_new + beta*p; % we will then initialize p value 
    q = A*p;
    r = r_new; % here we will output our r and z values to check for the relationship between convergence and maximum number of iterations
    z = z_new;
end
end

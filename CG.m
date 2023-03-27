function [ x, niters ] = CG( A, b, x0 )
% we can implement the Conjugate Gradient Method for solving Ax = b using
% A, b and x0
% we can define some variables as 
% A is the n x n matrix
% b is our n x 1 vector
% x0: n x 1 vector, hence this will be our initial value for x Ax = b
% niters: number of iterations required to converge
% k = niters
    n = length(b); % we will set n value to the dimension of length b 
    x0 = zeros(n, 1); % we will set x0 value 
    r = b - A * x0; % we will then make an equation to calculate our residual value 
    p = r; % we will initialize p = r 
    rsold = r' * r;  % we will then calculate our residual 
    niters = 0; % niter is the iteration counter 
% we will then form a for while loop for length b 
    for i = 1:length(b)
%while rsold > eps^2 && niters < n
        Ap = A * p; % we will make an equation to calculate Ap
        alpha = rsold / (p' * Ap); %computes the step size
        x = x + alpha * p; % updating the solution / iteration
        r = r - alpha * Ap; % caluclation of a residual
        rsnew = r' * r; % this is our new residual value 
        if sqrt(rsnew) < 1e-6 % here we set our tolerance value 
            break
        end
        p = r + (rsnew / rsold) * p;
        rsold = rsnew;   %updates the residual norm and increments the
        % iteration counter niters 
        niters = niters + 1; % finally we will print the iteration value to check for convergence in our solution 
    end
end
% The conjugate gradient function takes two inputs: the matrix A and 
% the vector b. In addition, it returns two important outputs including 
%vector x and the number of iterations which we named as niters in this 
%function.

function [ x, niters ] = CG( A, b, x0 )
% we can implement the Conjugate Gradient Method for solving Ax = b using
% A, b and x0
% we can define some variables as 
% A is the n x n matrix
% b is our n x 1 vector
% x0: n x 1 vector, hence this will be our initial value for x Ax = b
% niters: number of iterations required to converge
% k = niters
    n = length(b); 
    x0 = zeros(n, 1);
    r = b - A * x0;
    p = r;
    rsold = r' * r; 
    niters = 0;

    
    for i = 1:length(b)
%while rsold > eps^2 && niters < n
        Ap = A * p; %
        alpha = rsold / (p' * Ap); %computes the step size
        x = x + alpha * p; % updating the solution / iteration
        r = r - alpha * Ap; % caluclation of a residual
        rsnew = r' * r;
        if sqrt(rsnew) < 1e-10
            break
        end
        p = r + (rsnew / rsold) * p;
        rsold = rsnew;   %updates the residual norm and increments the
        % iteration counter niters 
        niters = niters + 1; 
    end
end
% The conjugate gradient function takes two inputs: the matrix A and 
% the vector b. In addition, it returns two important outputs including 
%vector x and the number of iterations which we named as niters in this 
%function
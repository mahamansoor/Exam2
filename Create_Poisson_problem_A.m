function [ A ] = Create_Poisson_problem_A( N )
%This matlab code will create a sparse matrix A of N^2 by N^2. This matrix
%A corresponds to the Possion matrix A(N). This code will set the diagonal
%entries of A to 4, the entries above and below the diagonal are set equal
%to -1 (except for the last column) and the entries to the left and right
%of diagonal are set to -1 (except for the last row). Hence, with this
%matlab code, the resulting matrix A is symmetric and positive definite. 
% N: size of the N x N mesh (meaning there are N^2 interior points)
% h: distance between grid points
h = 1 / (N+1); % we will set the distance between the grid points 
% The above variable will creates a (N^2) x (N^2) matrix of zeros
A = zeros(N^2);
% Populate the matrix with the stencil values and a matrix of zeros with
% dimensions (N^2) x N^2 is created. 
for i=1:N^2 % in this step we will apply the for loop 
    % we populate this matrix by stensil values for the poisson equation
    % with the help of nested loops. 
A(i,i) = 4; % here we will set the matrix A values 
% Check if the current point is on the left or right boundary
if mod(i-1,N) == 0 % we will then if mod (i-1, N) is equal to 0
% The point is on the left boundary
% Only include the stencil values to the right, up, and down of the
% current point (don't include the value to the left)
if i > N % the if statement is used to check if i > N
A(i,i-N) = -1; % this equation will initialize the A(i, i-N) equal to -1 
end
if i <= N*(N-1) % another if statement is applied to set the matrix diagonal entries 
A(i,i+N) = -1;
end
A(i,i+1) = -1;
elseif mod(i,N) == 0
% The point is on the right boundary
% Only include the stencil values to the left, up, and down of the
% current point (don't include the value to the right)
if i > N
A(i,i-N) = -1;
end
if i <= N*(N-1)
A(i,i+N) = -1;
end
A(i,i-1) = -1;
else
% The point is not on a boundary
% Include all stencil values
A(i,i-N) = -1;
A(i,i-1) = -1;
A(i,i+1) = -1;
A(i,i+N) = -1;
end
end
% Scale the matrix by h^2
A = (1 / h^2) * A; %h^2 values are used to scale the matrix in order to create the Poisson matrix A 
end
end % This code creates a matrix with size N with boundary points, size of grid, stencil values to creste poisson matrix A 
% can also be done by the second code 
function A = Create_Poisson_problem_A(N)
% we can create the matrix A for the Poisson problem on a unit square.
% we can use following variables for our input 
%   N: The size of the N x N mesh (meaning there are N^2 interior points).
% we can then compute the grid spacing.
h = 1/(N+1);
% we will then construct the diagonal entries.
diag_entries = -4*ones(N*N,1);
% we then construct the off-diagonal entries for the x direction.
x_entries = ones(N*N,1);
for i = 2:N
    x_entries((i-1)*N+1:i*N) = 1-h^2;
end
% Construct the off-diagonal entries for the y direction.
y_entries = ones(N*N,1);
for i = 1:N-1
    y_entries(i*N+1:(i+1)*N) = 1-h^2;
end
% Assemble the matrix A.
A = spdiags([diag_entries, y_entries, x_entries, y_entries, x_entries],...
            [-N,-1,0,1,N], N*N, N*N);
end % finally we will print the poisson matrix A(N)

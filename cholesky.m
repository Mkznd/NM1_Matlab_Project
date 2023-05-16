function [X, L, d] = cholesky(A, b)
% https://www.mathworks.com/help/simulink/mdl_gd/hi/matlab-functions.html
% FUNCTION NAME:
%   cholesky
%
% DESCRIPTION:
%   Cholesky factorization to find the solution to system of equations of
%   the form Ax = b and calculate the determinant of A
%
% INPUT:
%   A - (matrix Nx5) 
%   b - (matrix Nx1) RHS of the system 
%
% OUTPUT:
%   X - (matrix) Solution of the given system
%   d - (double) Calculated determinant of the given matrix
%
% ASSUMPTIONS AND LIMITATIONS:
%   A - SPD + pentadiagonal matrix
%   b - vector Nx1
%   Both A and b have N rows

[n, m] = size( A );
L = zeros( n, m );
% Perform decomposition
% For a formula and an example i referred to:
% https://atozmath.com/example/CONM/GaussEli.aspx?q=CD2&q1=E1
for i=1:n
   % populate diagonal entries
   L(i, i) = sqrt(A(i, i) - L(i, :)*L(i, :)');
   for k=(i + 1):n
      % populate lower entries
      L(k, i) = (A(k, i) - L(i,:)*L(k ,:)')/L(i, i);
   end
end
% Calculating the determinant
detL = 1.0;
detLT = 1.0;
for k = 1:n
    detL = detL * L(k,k);
end
for k = 1:n
    detLT = detLT * L(k,k)';
end
d = detL*detLT;
% Solving the system of equations
y = zeros(n, 1);
    y(1) = b(1) / L(1, 1);
    for k = 2:n
        y(k) = (b(k) - L(k, 1:(k-1)) * y(1:(k-1))) / L(k, k);
    end
    
    X = zeros(n, 1);
    X(n) = y(n) / L(n, n);
    for k = (n-1):-1:1
        X(k) = (y(k) - L((k+1):n, k)' * X((k+1):n)) / L(k, k);
    end
end
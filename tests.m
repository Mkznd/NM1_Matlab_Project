lower_rand = -20; % lower bound of matrix elements
upper_rand = 20; % upper bound of matrix elements
matrix_size = 5; % matrix size

% Generate SPD matrix
ii = rand(matrix_size);
A = ii*ii.';

%Generate_solution
x_known = randi([lower_rand, upper_rand], matrix_size, 1);

%Generate RHS
b = A*x_known;

% Find what we need to run the tests
det_known = det(A);
L_known = chol(A, "lower");
I = eye(matrix_size);

% Use the function
[x_calculated, L_calculated, det_calculated] = cholesky(A, b);

% Run the tests
condition_number = norm(inv(A)) * norm(A)
solution_error = norm(x_calculated - x_known)/ norm(x_known)
forward_stability_error = norm(x_calculated - x_known)/ (norm(x_known) * condition_number)
inverse_error = norm(I - x_calculated.' * A)/(norm(A)*norm(x_calculated))
determinant_error = abs((det_calculated - det_known)/(det_known))
factorization_error = norm(A - L_calculated*L_calculated.') / norm(A)
%%%%%%%%%%%
% tests three models: simplex, dentzig-wolfe and linprog on different
% number of cases

% outputs objective function value and runtime index for each model and store in list
%%%%%%%%%%%


% construct net of case_count that will be evaluated 
increment1 = 2;
increment2 = 10;
case_counts = [1:increment1: 99, 100:increment2:300];

%%%%% for testing %%%%%%
%case_test = [1, 50, 100];
%case_counts = case_test;

% parameter setup for the dantzig wolfe model
M = 10000;
termi_limit = 1e-6;


% log the object values at optimality
obj_simplex_lst = zeros(length(case_counts),1);
obj_dantzig_lst = zeros(length(case_counts),1);
obj_linprog_lst = zeros(length(case_counts),1);

% log the solutions of x s at optimality
x_simplex_lst = zeros(length(case_counts),1);
x_dantzig_lst = zeros(length(case_counts),1);
x_linprog_lst = zeros(length(case_counts),1);

% log the runtime to solve the LP
runtime_simplex_lst = zeros(length(case_counts),1);
runtime_dantzig_lst = zeros(length(case_counts),1);
runtime_linprog_lst = zeros(length(case_counts),1);

for n = case_counts
    % generate the matrix A, b, c
    [A_i,b_i,c_i] = coefficient_generator(n);
    
    % idx for logging obj_val, x_opt and runtime
    i = find(n == case_counts, 1);

    % simple method
    tic;
    [obj_simplex, x_simplex] = simplex_solver(M, A_i,b_i,c_i); 
    obj_simplex_lst(i, :) = obj_simplex;
    %x_simplex_lst(i, :) = x_simplex;
    runtime_simplex_lst(i) = toc;
    
    % dantzig wolfe
    tic;
    [obj_dantzig, x_dantzig] = dantzig_solver(M, A_i, b_i, c_i, n, termi_limit);
    obj_dantzig_lst(i, :) = obj_dantzig;
    %x_dantzig_lst(i, :) = x_dantzig;
    runtime_dantzig_lst(i) = toc;
    
    % linprog
    tic;
    [x_linprog, obj_linprog] = linprog(c_i,A_i,b_i,[],[],zeros(1,3+6*n),[]);
    obj_linprog_lst(i, :) = obj_linprog;
    %x_linprog_lst(i, :) = x_linprog;
    runtime_linprog_lst(i) = toc;

end
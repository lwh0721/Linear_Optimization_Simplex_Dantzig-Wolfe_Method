function [obj_opt, x_opt] = simplex_solver(M, A, b, c)
    % simplex method model
    %
    % input: 
    %       A, b ,c for std LP
    %       min c*x
    %       s.t. A*x <= b
    %            x >= 0
    %
    %
    % output: 
    %       obj_opt: objective value associated with optimal solution,
    %       minimum cost/ maximum profit in this case 
    %       x_opt: optimal solution x

   
    %%%%%%%
    % standardize the problem by adding slack variables
    
    %%%%%%%%%%%%%%%%%  A,b  %%%%%%%%%%%%%%%%%%

    % extract the size of the A matrix
    [m_A, n_A] = size(A);

    % add an identity matrix to A for the slack variables
    I_A = eye(m_A);
    A_std = [A, I_A];

    % convert the RHD of the constraints to positive and update the A and
    % b matrix
    position_negative = (b > 0) - (b <= 0);
    % update the A and b matrix where the b values are negative
    A_std = A_std.*position_negative;
    b_std = b.*position_negative; %%%%%%%%%%%%%%%%%%%% b %%%%%%%%%%%%%%%%
    
    % introduce the big M variables
    I_bigM = eye(m_A); % number of big M variables is equal to the number of constraints 
    % append indentity matrix for big M variable 
    A_std = [A_std, I_bigM]; 
    
    %%%%%%%%%%%%%%%%%  c  %%%%%%%%%%%%%%%%%%

    % add 1*m zeros for the slack variables in c
    c_slack = zeros(1,m_A);
    c_std = [c, c_slack];

    % add coeffecients of the slack variables to the objective coefficients
    c_bigM = M*ones(1,m_A);  % number of big M variables is equal to the number of constraints 
    % append coefficient matrix for slack variable 
    c_std = [c_std, c_bigM]; 
    
    %Update the dimensions of A of the std LP
    [m_A, n_A] = size(A_std);
    

    %%%%%
    % Step 0: initialization
    
    %%%%%%% x %%%%%%%%%%

    % solve A*x = b  for the initial basic solutions
    % extract the portion of A associated with slack variables
    A_xb = A_std(:,(n_A-m_A+1):end);
    % solve for x_B
    x_B = linsolve(A_xb,b_std);

    % for initial solutions, x_N are zeros 
    % number of elements in x_N
    xn_count = n_A-length(x_B);
    %set x_N to xn_count*1 zeros
    x_N = zeros(xn_count,1);
    
    %%%%%%% N_bar, B_bar, N, B, c_N, c_B %%%%%%%%

    % N_bar, B_bar
    N_bar = 1:length(x_N);
    x_B_range = 1:length(x_B);
    B_bar = length(x_N) + x_B_range;

    % B and N matrix
    B = A_std(:,B_bar);
    N = A_std(:,N_bar);
    
    % c_B and c_N
    c_B = c_std(B_bar);
    c_N = c_std(N_bar);
    
    % calculate initial objective function value
    obj_opt = c_N*x_N + c_B*x_B;
    
    % track number of iterations
    iteration = 0;

    while iteration >= 0

        %%%%%%
        % Step1: check if reached optimality 
        % calculate r_N
       
        cb_b_inverse = linsolve(B.',c_B.');
        r_N = c_N - cb_b_inverse.'*N;

        % check if there's negative r_N values
        if all(r_N >= 0) 
        % if all r_N are non_negative
        % reached optimality, break out of the loop
            break
        end

        % else there's nagative r_n values
        % select the first non_positive element in the non-basis as
        % entering variable
        negative_rn_idx = N_bar(r_N < 0);
        entering_index = negative_rn_idx(1);

        %%%%%%%%
        % Step 2: calculate d
        N_q = N(:,N_bar == entering_index);

        % calcualte -b inverse n_q
        dq_basis = linsolve(B, -N_q);
        % set very small elements in d to zero to correct inaccuracy
        dq_basis(abs(dq_basis) < 0.00001) = 0; 

        % calculate e_q
        I_eq = eye(length(x_N)); 
        % update the entering variable position to 1
        e_q = I_eq(:,N_bar == entering_index);  


        % check dq for negative vables 

        if all(dq_basis >= 0)
        % if all the values are non_negative
        % then the problem is unbounded
            obj_opt = 1234567890;
            
            % size of the A matrix 
            [m, n] = size(A);
            % set the optimal solution to 0 
            x_opt = zeros(n,1);
            
            % construct the extreme direction
            for j = 1:n
                % set the ouput solution x to d_q

                if any(B_bar == j)
                    % the bais portion of dq
                    x_opt(j) = dq_basis(B_bar == j);

                else
                    % the non_basis portion correspondind to e_q
                    x_opt(j) = e_q(N_bar == j);

                end
            end

            return
        end
        


        %%%%%%%%
        % Step 3: calculate d

        % index of basis variable with negative d1 value
        negative_d_index = B_bar(dq_basis < 0);

        % minimum ratio test 
        x_j = x_B(dq_basis < 0);
        d_j = dq_basis(dq_basis < 0);
        alpha_lst = -x_j ./ d_j;

        % minimum value is alpha 
        alpha = min(alpha_lst);
        
        %%%%%%%%%
        % Step4: update x
        x_B = x_B + alpha * dq_basis;
        x_N = x_N + alpha * e_q;
        
        %%%%%%%%
        % Step 5: update x, N_bar, b_bar, N, B, c_B, c_N

        % index of the exiting variable
        exiting_index = negative_d_index(alpha_lst == alpha);
        exiting_index = min(exiting_index);
        
        % update x_B and x_N
        x_B(B_bar == exiting_index) = x_N(N_bar == entering_index);
        x_N(N_bar == entering_index) = 0;
        
        % update B and N
        B_exiting = B(:,B_bar == exiting_index);
        B(:,B_bar == exiting_index) = N_q;
        N(:,N_bar == entering_index) = B_exiting;

        % update c_B and c_N
        c_exiting = c_B(B_bar == exiting_index);
        c_B(B_bar == exiting_index) = c_N(N_bar == entering_index);
        c_N(N_bar == entering_index) = c_exiting;
        
        % update B_bar and N_bar
        B_bar(B_bar == exiting_index) = entering_index;
        N_bar(N_bar == entering_index) = exiting_index;


           
        % update the objective function value
        obj_opt = c_N*x_N + c_B*x_B;

        % log the iteration count
        iteration = iteration + 1;
        
    end

    
    % form the optimal solution x 
    x_opt = zeros(n_A,1);
    
    % update the x_N portion and the x_B portion individually
    for j = 1:n_A
        if any(B_bar == j)
            % the basis portion of x
            x_opt(j) = x_B(B_bar == j);

        else
            % the non_basis portion of x
            x_opt(j) = x_N(N_bar == j);

        end

    end

end
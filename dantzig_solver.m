function [obj_opt, x_opt] = dantzig_solver(M, A, b, c, case_count, termi_limit)
    % dantzig_wolfe model
    %
    % input: 
    %       M: the value of the big M
    %       A, b ,c for std LP
    %       min c*x
    %       s.t. A*x <= b
    %            x >= 0
    %
    %       case_count = number of cases to be considered, consistent with
    %       number of cases when generating the matrix
    %       termi_limit: apply threshold on metrics to terminate the
    %       algorithm 
    %
    %
    % output: 
    %       obj_opt: objective value associated with optimal solution,
    %       minimum cost/ maximum profit in this case 
    %       x_opt: optimal solution x
    

    
    % initialize the objects (optimal solution and object value for return
    % from the function 
    obj_opt = [];
    x_opt = [];
    
    %%%%%%%%%%%%%%
    % Extract the constraints for the master problem 

    % extract the size of the A matrix
    [m_A, n_A] = size(A);

    % initialize objects for extracting the A and b matrix associated with
    % the master problem
    A_master = zeros(3*case_count, n_A);
    b_master = zeros(3*case_count, 1);
    
    % loop thru each case to extract the A and b from the generated A and b
    % matrix 
    for i = 1:case_count

        % index for extracting linking constraints from A and b
        idx_original_l = 4*i-2;
        idx_original_u = 4*i;
        idx_master_l = 3*i-2;
        idx_master_u = 3*i;
        
        % A and b matrix for the master problem 
        A_master((idx_master_l):(idx_master_u),:) = A((idx_original_l):(idx_original_u),:);
        b_master((idx_master_l):(idx_master_u)) = b((idx_original_l):(idx_original_u));
    end
    
    %%%%%%%%%%%%%%%%%%
    % extract constraints for the subproblems 
    
    % initialize objects for extracting the A and b matrix associated with
    % the subproblems
    A_sub = cell(1,case_count+1);
    b_sub = cell(1,case_count+1);
    c_sub = cell(1,case_count+1);
    
    %extract constraints for the first subproblem
    A_sub{1} = A(1,1:3);
    b_sub{1} = b(1);
    c_sub{1} = c(1:3);

    % loop thru each case to extract the A and b and c from the generated A and b
    % matrix
    for i = 1:case_count
        % index for extracting constraints for the subproblems
        idx1 = 4*i+1;
        idx2 = (6*i)-2;
        idx3 = (6*i)+3;

        %extract constraints for the rest of the subproblems
        A_sub{i+1} = A(idx1, idx2:idx3);
        b_sub{i+1} = b(idx1);
        c_sub{i+1} = c(idx2:idx3);
    end
    
    %%%%%%%%%
    % Step 0: initilization
    
    % initialize objects for extreme points of each
    % subproblem and track which subproblem the extreme points are from 
    v = cell(1,case_count+1); % extreme points 
    sub_index = zeros(1,case_count+1); % index of subproblem that the extreme points are generated from 
    
    % track where each extreme points are from and initialize all these
    % extreme points to 0
    for i = 1:(case_count+1)
        sub_index(i) = i;
        v{i} = [zeros(length(c_sub{i}),1)];
      
    end
    
    %%%%%%%%%%
    % initialize the constraint for the master problem


    %%%%%%%%%%%%%%%% A, b %%%%%%%%%%%%%
    
    % extract the size of the A_master matrix for indexing
    [m_A_master, n_A_master] = size(A_master);

    % introduce the slack variables to convert the master problem from
    % inequality to equality
    I_slack = eye(m_A_master); % number of slack variables is equal to the number of constraints in the master problem
    % append indentity matrix for slack variable 
    A_master = [A_master, I_slack]; 

    % The RHS matrix of the constraints for the master problem are all non
    % positive 
    % multiply matrix A and b to convert the RHS to non-negative
    A_master = -1 * A_master;
    b_master = -1 * b_master;  %%%%%%%%%%%  b %%%%%%%%%%%%
    
    % introduce the big M variables
    I_bigM = eye(m_A_master); % number of big M variables is equal to the number of constraints in the master problem
    % append indentity matrix for big M variable 
    A_master = [A_master, I_bigM]; 
    
    %%%%%%%%%%%  c %%%%%%%%%%%%

    % add coeffecients of the slack variables to the objective coefficients
    c_slack = zeros(1,m_A_master); % number of slack variables is equal to the number of constraints in the master problem
    % append coefficient matrix for slack variable 
    c = [c, c_slack];
    
    % add coeffecients of the big M variables to the objective coefficients 
    c_bigM = M*ones(1,m_A_master);  % number of M variables is equal to the number of constraints in the master problem
    % append coefficient matrix for big M variable 
    c = [c, c_bigM];

    %%%%%%%%%%%  L %%%%%%%%%%%%
    
    % extract the L matrix from the A_master
    % initialize objects for L matrix 
    L = cell(1,case_count+1);

    %extract constraints for the first subproblem
    L{1} = A_master(:,1:3);

    % loop thru each case to extract the L from A_master
    for i = 1:case_count
        L{i+1} = A_master(:,((6*i)-2):((6*i)+3));
    end
    

    %%%%%%%%%%%%% 
    % initial basis matrix 
    
    % the index of the basis set 
    B_idx_l = n_A_master+m_A_master+1;
    B_idx_u = n_A_master+m_A_master+m_A_master;
    B_idx = [sub_index, B_idx_l:B_idx_u];
    
    % initialize the basis matrix with 0
    B = zeros(m_A_master+case_count+1, length(B_idx));

    % initialize the fb with 0
    f_B = zeros(length(B_idx),1);

    % loop thru to update and B matrix and fb matrix 
    for idx_fb = 1:length(B_idx)
        % the portion associated with lambda and mu variables 
        if B_idx(idx_fb) <= (case_count+1)
            % coefficients for the lamba values
            l_v = L{idx_fb} * v{idx_fb};
            % coefficients for the s variables
            s = zeros(case_count+1,1);
            s(idx_fb) = 1;
            % update the B matrix 
            B(:,idx_fb) = [l_v; s];
            % calculate the fb value
            f_B(idx_fb) = c_sub{idx_fb}*v{idx_fb};

        % the portion associated with big M variables 
        else
            % coefficient for the big M variables 
            B_left = A_master(:,B_idx(idx_fb));
            B_right = zeros(case_count+1,1);
            % update the B matrix 
            B(:,idx_fb) = [B_left; B_right];
            % calculate the fb value
            f_B(idx_fb) = M;
        end
    end
    
    %%%%%%%%%%% set the initial value for x 
    
    % the left portion are zeros
    % set the x associates with big M variables to b_master
    x_B = [ones(case_count+1,1);b_master];
    
    % calculate the initial object funtion value
    obj_opt = f_B'*x_B;

    % track number of iterations
    iteration = 0;

    while iteration >= 0
        
        %%%%%%%% 
        % Step1: calculate pi

        pi = linsolve(B', f_B);
        
        % split pi into pi1 and pi2 
        pi1 = pi(1:m_A_master);
        pi2 = pi(m_A_master+1:end);
        
        %%%%%
        %Step2: check if reached optimality 

        % initialize the r* with 0
        r_star = zeros(1,case_count+1);
        
        % track the extreme direction 
        extreme_found = 0; % whether the extreme direction is found 
        entering_idx = 0; % index of the entering variable

        %log the fb value for updating 
        fb_exit = 0;
        
        % log the solution to the subproblems
        x_SP = cell(1,case_count+1);
        
        % initialize a 
        a = zeros(m_A_master+case_count+1,1);
        
        for j = 1:(case_count+1)
            %In certain cases, obj is producing nub_masterers such as 2e-14
            %instead of 0. To allow the algorithm to terminate, we will set
            %these to 0 based on the tolerance input.

            % c matrix for the subproblem 
            c_SP = c_sub{j} - pi1' * L{j};
            % for c values smaller than the threshold, set to 0 to allow
            % the optimization to terminate
            c_SP(abs(c_SP) < termi_limit) = 0;

            %[r, dir] = linprog(c_SP, [], [], A_sub{j}, b_sub{j}, [], [], options);

            % solve the subproblem 
            [obj_sub, x_sub] = simplex_solver(M, A_sub{j}, b_sub{j}, c_SP); 

            % in case where the subproblem is unbounded
            % break and return the x_sub
            if obj_sub == 1234567890
                % update the basis portion of a 
                a(1:m_A_master) = L{j} * x_sub;
                % update status variable 
                extreme_found = 1;
                % log the index of the entering variable
                entering_idx = j;
                % update the extreme directions
                x_SP{j} = x_sub;
                % log the fb_value 
                fb_exit = c_sub{j} * x_sub;

                % update r*
                r_star(j) = obj_sub;
                break

            % else update r* and and the subproblem solution
            else
                r_star(j) = obj_sub - pi2(j);
                x_SP{j} = x_sub;
            end
            
        end
        
        %%%%%
        %Step3: calculate r min and a

        if extreme_found == 0
        % if did not find extreme direction
        % use the smallest r value from the subproblems 
            r_min = min(r_star);
            if (r_min >= 0 || abs(r_min) <= termi_limit)
                % in the case where minimum r value is non neghative 
                % or its smaller than the threshold value
                % break out of the loop 
                break
            end
            % index associated with the minimum r value is the entering
            % variable
            r_min_idx = find(r_star == r_min, 1);
            entering_idx = r_min_idx;
            % calculate a
            a(1:m_A_master) = L{r_min_idx} * x_SP{r_min_idx};
            % update e vector in the a matrix 
            a(m_A_master+r_min_idx) = 1;
            % log the fb_value 
            fb_exit = c_sub{r_min_idx} * x_SP{r_min_idx};
        end
        
        %%%%%
        %Step4: calculate d

        % solve Bd = -a for direction d
        d = linsolve(B,-a);
        % apply the threshold and anything lower than the threshold is set
        % to 0 in the matrix d
        d(abs(d) < termi_limit) = 0; 

        % check direction d for negative values
        % if all_non negative, the problem is unbounded 
        if all(d >= 0)
            disp("LP is unbounded");
            break
        end
        
        %%%%%
        %Step5: calculate alpha
        %extract the index of d_L that is negative 
        d_negative_index = find(d < 0);
        d_l = d(d_negative_index);
        x_l = x_B(d_negative_index);

        % calculate x_l/ d_l for all negative d's
        alphas = -x_l./d_l;
        % find the minimum value and set as alpha
        alpha = min(alphas);
        
        % find the index for updating the B matrix 
        exiting_index = find(alphas == alpha, 1); % the index of the exiting varible
        B_update_index = d_negative_index(exiting_index);
        
        %%%%%
        %Step6: update x

        x_B = x_B + alpha * d;
        
        %%%%%
        %Step7: update basis and x, extreme points, fb, obj values

        % update x
        x_B(B_update_index) = alpha;

        % update the basis matrix
        B(:,B_update_index) = a;
        % update the basis index
        B_idx(B_update_index) = entering_idx;
        
        % update fb
        f_B(B_update_index) = fb_exit;

        % update the extreme direction
        v{B_update_index} = x_SP{entering_idx};
        % update the objective function value
        obj_opt = f_B' * x_B;


        % log the iteration count
        iteration = iteration + 1;
    end
    
    % form the optimal solution x 
    for idx_x = 1:(case_count+1)
        % the index of the x associated with the n cases 
        primary_idx = find(B_idx == idx_x);

        % first case append 3 zeros
        if idx_x == 1
            x_values = zeros(3,1);
        % other cases append 6 zeros
        else
            x_values = zeros(6,1);
        end

        % for the idx_x index , calculate x_B * v
        for l = primary_idx
            contribution = x_B(l) * v{l};
            % append to the zeros
            x_values = x_values + contribution;
        end

        % output x
        x_opt = [x_opt; x_values];

    end
end
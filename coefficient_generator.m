function [A, b, c] = coefficient_generator(n)
    % input: 
    %       n: number of cases to be considered 
    % output: 
    %       A, b ,c for std LP
    %       min c*x
    %       s.t. A*x = b
    %            x >= 0

    %%%%%%%%%%%%%%%%%%%%%    A   %%%%%%%%%%%%%%%%%%%%%%%

    %%%%%% x1, x2, x3 %%%%%%
    % the constraint on total planting area is common across all the
    % cases,construct the A for the first constraint 
    
    %%%%%%% for the first constraint total planting area 
    % x1 + x2 + x3 <= 500
    A_common_x = [1, 1, 1]; 
    

    %%%%%%% for the constraint 2-5 
    %       - 2.5x1 - y1 + w1 <= -200
    %       -3x2 - y2 + w2 <= -240    
    %       -20x3 + w3 + w4 <= 0 
    %       w3 <= 6000
    
    % for each case the coefficient for x1, x2, x3 takes into consideration
    % the yields 
    A_x = A_common_x;

    % loop through the n cases
    for i = 1:n 
        % generate random number in the range of 0.5 -1.5 for yield 
        yield = 0.5 + rand(1); 
        % coefficient of x1, x2, x3 for constraint 2-4
        yields = -diag(yield .* [2.5, 3, 20]); 
        % coefficient of x1, x2, x3 for constraint 5
        fill = zeros(1,3);
        A_x = [A_x; yields; fill];
    end

    %%%%%% y1, y2 ,w1, w2, w3, w4 %%%%%%

    % only constraint 2-5 applies 

    % now look at the rest of the constraints, these 4 constraints will be
    % applied for each case 
    %       - 2.5x1 - y1 + w1 <= -200
    %       -3x2 - y2 + w2 <= -240    
    %       -20x3 + w3 + w4 <= 0 
    %       w3 <= 6000

    % notice that the x1 and x2 variable are common for all cases but the 
    % y1, y2 ,w1, w2, w3, w4 variables are different for each case
    % (subscripts for each case)

    % construct the A matrix for the y1, y2 ,w1, w2, w3, w4 variables
    % for an individual case
    
    %     y1 y2 w1 w2  w3 w4
    A_yw = [-1, 0, 1, 0, 0, 0; %Wheat Production Constraint
          0, -1, 0, 1, 0, 0; %Corn Production Constraint
          0, 0, 0, 0, 1, 1; %Beet Production Constraint
          0, 0, 0, 0, 1, 0]; %High Price Beet Sales Constraint

    % concat all n cases for variable y1, y2 ,w1, w2, w3, w4 in constaints
    % 2-5
    A_yw_n = kron(eye(n),A_yw);

    % build a 1*6n matrix to fill up the space above matrix A_yw_n for
    % constraint 1
    A_fill = zeros(1, 6*n);
    
    % concat the fill matrix and the n cases yw matrix to for the entire yw
    % portion of A matrix
    A_yw_all = [A_fill; A_yw_n];

    %%%%%%%% concat the x portion with the yw portion to form A matrix
    A = [A_x, A_yw_all];
    
    %%%%%%%%%%%%%%%%%%%%%    B   %%%%%%%%%%%%%%%%%%%%%%%
    % constraint 1 is common for all cases 
    b_common = [500]; 

    % constraint 2-5 repeated n times for n cases 
    b_n= repmat([-200; -240; 0; 6000],n,1); % concat along axis = 0

    %%%%%%%% concat and form the b matrix
    b = [b_common; b_n];
    
    %%%%%%%%%%%%%%%%%%%%%    C   %%%%%%%%%%%%%%%%%%%%%%%
    % cost coefficient for x1, x2, x3 are the same for n cases 
    c_common = [150, 230, 260]; 

    % cost coefficient for y1, y2, w1, w2, w3, w4 repated n times and take
    % average of all n cases 
    c_n = repmat([238, 210, -170, -150, -36, -10],1,n)/n; % concat along axis = 1

    %%%%%%%% concat and form the c matrix
    c = [c_common, c_n];

end
    

    

    




    
    

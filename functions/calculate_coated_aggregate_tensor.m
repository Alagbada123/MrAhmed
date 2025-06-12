function Ainf = calculate_coated_aggregate_tensor(Chom_cp, stiffness_coated, fractions_coated, I, J, K)
% CALCULATE_COATED_AGGREGATE_TENSOR Calculates strain concentration tensor for 
% coated aggregate with multiple coating layers using Herve-Zaoui solution
%
% Inputs:
%   Chom_cp         - Homogenized stiffness of cement paste (matrix)
%   stiffness_coated - Single structure from stiffness.coated array
%   fractions_coated - Single structure from fractions.coated array
%   I, J, K         - Fundamental tensors
%
% Output:
%   Ainf            - Structure containing strain concentration tensors for core and coatings

    % Get number of coatings
    n_coatings = length(stiffness_coated.coating);
    
    % Calculate cumulative volumes for radius calculation
    V_cumulative = zeros(n_coatings + 2, 1);  % +2 for core and matrix
    V_cumulative(1) = fractions_coated.core;
    
    % Add volumes of coatings
    for i = 1:n_coatings
        V_cumulative(i+1) = V_cumulative(i) + fractions_coated.coating(i);
    end
    % Add cement paste volume
    V_cumulative(end) = V_cumulative(end) + fractions.cement_paste;
    
    % Calculate radii
    R = zeros(n_coatings + 2, 1);
    for i = 1:length(V_cumulative)
        R(i) = (1/2)*6^(1/3)*V_cumulative(i)^(1/3)/pi^(1/3);
    end
    
    % Prepare input matrix for Herve-Zaoui solution
    input_matrix = zeros(n_coatings + 2, 3);
    
    % Core properties
    input_matrix(1,:) = [stiffness_coated.core.k, stiffness_coated.core.mu, R(1)];
    
    % Coating properties
    for i = 1:n_coatings
        input_matrix(i+1,:) = [
            stiffness_coated.coating(i).k, ...
            stiffness_coated.coating(i).mu, ...
            R(i+1)
        ];
    end
    
    % Matrix properties (from Chom_cp)
    k_matrix = Chom_cp(1,1) - 4/3*(Chom_cp(6,6)/2);
    mu_matrix = Chom_cp(6,6)/2;
    input_matrix(end,:) = [k_matrix, mu_matrix, R(end)];
    
    % Apply Herve-Zaoui solution
    sol_int = fun_HZ_int(input_matrix);
    
    % Extract strain concentration tensors for all phases
    % Core
    Ainf.core = sol_int(1,5)*J + sol_int(1,7)*K;
    
    % Coatings
    Ainf.coating = [];
    for i = 1:n_coatings
        Ainf.coating(i).A = sol_int(i+1,5)*J + sol_int(i+1,7)*K;
    end
end
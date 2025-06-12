function stiffness = calculate_stiffness_matrices(materialProperties, J, K)
% CALCULATE_STIFFNESS_MATRICES Computes stiffness matrices for all phases
%
% Inputs:
%   materialProperties - Structure containing material properties
%   J, K              - Volumetric and deviatoric projection tensors
%
% Output:
%   stiffness         - Structure containing stiffness matrices for all phases

    % Calculate clinker stiffness
    stiffness.clinker = 3*materialProperties.clinker.k*J + 2*materialProperties.clinker.mu*K;
    
    % Calculate hydrate stiffness
    stiffness.hydrate = 3*materialProperties.hydrate.k*J + 2*materialProperties.hydrate.mu*K;
    
    % Calculate uncoated aggregate stiffnesses
    uncoated_names = fieldnames(materialProperties.uncoated_aggregates);
    for i = 1:length(uncoated_names)
        agg_name = uncoated_names{i};
        agg = materialProperties.uncoated_aggregates.(agg_name);
        stiffness.uncoated_aggregates.(agg_name) = 3*agg.k*J + 2*agg.mu*K;
    end
    
    % Calculate coated aggregate stiffnesses if present
    if isfield(materialProperties, 'coated_aggregates')
        coated_names = fieldnames(materialProperties.coated_aggregates);
        for i = 1:length(coated_names)
            agg_name = coated_names{i};
            agg = materialProperties.coated_aggregates.(agg_name);
            stiffness.coated_aggregates.(agg_name) = 3*agg.k*J + 2*agg.mu*K;
        end
    end
    
    % % Pore stiffness (zero)
    % stiffness.pore = zeros(6,6);

    k_H2O = 2.3;
    mu_H2O = 0;
    stiffness.pore = 3*k_H2O*J+2*mu_H2O*K;
    
    % Verify all stiffness matrices are 6x6
    stiffness = verify_stiffness_dimensions(stiffness);
end

function stiffness = verify_stiffness_dimensions(stiffness)
    % Verify dimensions of all stiffness matrices
    fields = fieldnames(stiffness);
    for i = 1:length(fields)
        if isstruct(stiffness.(fields{i}))
            % Recursive check for nested structures
            stiffness.(fields{i}) = verify_stiffness_dimensions(stiffness.(fields{i}));
        else
            assert(all(size(stiffness.(fields{i})) == [6,6]), ...
                ['Stiffness matrix for ', fields{i}, ' must be 6x6']);
        end
    end
end
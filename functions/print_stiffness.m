function print_stiffness(stiffness)
% PRINT_STIFFNESS Prints stiffness matrices for all phases in matrix form
%
% Input:
%   stiffness - Structure containing stiffness matrices for all phases

    fprintf('\n=== Stiffness Matrices Summary ===\n\n');
    
    % Print clinker stiffness
    fprintf('Clinker Stiffness:\n');
    print_matrix(stiffness.clinker);
    fprintf('\n');
    
    % Print hydrate stiffness
    fprintf('Hydrate Stiffness:\n');
    print_matrix(stiffness.hydrate);
    fprintf('\n');
    
    % Print uncoated aggregate stiffnesses
    fprintf('Uncoated Aggregate Stiffnesses:\n');
    uncoated_names = fieldnames(stiffness.uncoated_aggregates);
    for i = 1:length(uncoated_names)
        fprintf('%s:\n', uncoated_names{i});
        print_matrix(stiffness.uncoated_aggregates.(uncoated_names{i}));
        fprintf('\n');
    end
    
    % Print coated aggregate stiffnesses if present
    if isfield(stiffness, 'coated_aggregates')
        fprintf('Coated Aggregate Stiffnesses:\n');
        coated_names = fieldnames(stiffness.coated_aggregates);
        for i = 1:length(coated_names)
            fprintf('%s:\n', coated_names{i});
            print_matrix(stiffness.coated_aggregates.(coated_names{i}));
            fprintf('\n');
        end
    end
    
    % Print pore stiffness
    fprintf('Pore Stiffness:\n');
    print_matrix(stiffness.pore);
    fprintf('\n');
end

function print_matrix(C)
    % Helper function to print 6x6 stiffness matrix in matrix form
    if all(C(:) == 0)
        fprintf('  [Zero matrix]\n');
    else
        fprintf('  [\n');
        for i = 1:6
            fprintf('   ');
            for j = 1:6
                if abs(C(i,j)) < 1e-10
                    fprintf('    0.00');
                else
                    fprintf('%7.2f', C(i,j));
                end
            end
            fprintf('\n');
        end
        fprintf('  ]\n');
    end
end
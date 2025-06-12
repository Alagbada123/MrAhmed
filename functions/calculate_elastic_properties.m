function elastic = calculate_elastic_properties(C, phase_name, print_results)
% CALCULATE_ELASTIC_PROPERTIES Calculates elastic properties from stiffness tensor
%
% Inputs:
%   C            - Homogenized stiffness tensor
%   phase_name   - Name of the phase (optional, for printing)
%   print_results - Boolean flag to control printing (optional)
%
% Output:
%   elastic     - Structure containing elastic properties (k, mu, nu, E)

    if nargin < 2
        phase_name = '';
        print_results = false;
    elseif nargin < 3
        print_results = false;
    end

    % Calculate elastic properties
    elastic.k = C(1,1) - 4/3*(C(6,6)/2);    % Bulk modulus
    elastic.mu = C(6,6)/2;                   % Shear modulus
    elastic.nu = (3*elastic.k - 2*elastic.mu)/(6*elastic.k + 2*elastic.mu);  % Poisson's ratio
    elastic.E = 9*elastic.k*elastic.mu/(3*elastic.k + elastic.mu);           % Young's modulus
    
    % Print results if requested
    if print_results
        if ~isempty(phase_name)
            fprintf('\n%s Properties:\n', phase_name);
        else
            fprintf('\nElastic Properties:\n');
        end
        fprintf('Bulk Modulus (K): %.2f GPa\n', elastic.k);
        fprintf('Shear Modulus (μ): %.2f GPa\n', elastic.mu);
        fprintf('Young''s Modulus (E): %.2f GPa\n', elastic.E);
        fprintf('Poisson''s Ratio (ν): %.3f\n\n', elastic.nu);
    end
end
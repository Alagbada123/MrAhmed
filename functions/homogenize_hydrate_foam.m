function [Chom_hf, EEinfty_hf_iso, elastic_hf] = homogenize_hydrate_foam(fractions, stiffness, params, I, J, K, foam_type)
% HOMOGENIZE_HYDRATE_FOAM Computes homogenized properties of hydrate foam using
% self-consistent scheme
%
% Inputs:
%   fractions  - Structure containing volume fractions
%   stiffness  - Structure containing stiffness matrices
%   params     - Analysis parameters including tolerance
%   I, J, K    - Fundamental tensors
%   foam_type  - Type of hydrate foam ('normal' or 'ITZ')
%
% Outputs:
%   Chom_hf        - Homogenized stiffness tensor
%   EEinfty_hf_iso - Strain concentration tensor
%   elastic_hf     - Structure containing elastic properties (K, G, E, nu)

    fprintf('\n=== Hydrate Foam Homogenization (%s Foam) ===\n', foam_type);
    
    % Select hydrate foam type
    switch foam_type
        case 'normal'
            foam = fractions.hydrate_foam;
        case 'ITZ'
            foam = fractions.ITZ_hydrate_foam;
        otherwise
            error('Invalid foam_type. Use ''normal'' or ''ITZ''.');
    end
    
    % Initialize with hydrate stiffness weighted by volume fraction
    C0 = foam.hydrates * stiffness.hydrate;
    deviation = 1;
    
    % Iterative self-consistent scheme
    iteration = 0;
    while deviation > params.tolerance && iteration < 100
        iteration = iteration + 1;
        
        % Calculate Hill tensor for spherical pores
        P_pore = fun_P_sphere_iso(C0);
        
        % Calculate strain concentration tensors
        Ainf_pore = inv(I + P_pore*(stiffness.pore - C0));
        Ainf_hyd = fun_Ainf_needle_iso(stiffness.hydrate, C0);
        
        % Calculate overall strain concentration tensor
        EEinfty_hf = inv(foam.porosity*Ainf_pore + ...
                        foam.hydrates*Ainf_hyd);
        
        % Update strain concentration tensor for hydrates
        A_hyd = Ainf_hyd * EEinfty_hf;
        
        % Calculate homogenized stiffness
        Chom_hf_new = foam.hydrates * stiffness.hydrate * A_hyd;
        
        % Check convergence
        deviation = abs(norm(Chom_hf_new - C0)/norm(Chom_hf_new));
        
        % Print iteration progress
        fprintf('Iteration %d: Convergence = %.2e\n', iteration, deviation);
        
        % Update for next iteration
        C0 = Chom_hf_new;
    end
    
    % Store final results
    Chom_hf = C0;
    EEinfty_hf_iso = EEinfty_hf;
    
    % Calculate elastic properties
    elastic_hf.k = Chom_hf(1,1) - 4/3*(Chom_hf(6,6)/2);
    elastic_hf.mu = Chom_hf(6,6)/2;
    elastic_hf.nu = (3*elastic_hf.k - 2*elastic_hf.mu)/(6*elastic_hf.k + 2*elastic_hf.mu);
    elastic_hf.E = 9*elastic_hf.k*elastic_hf.mu/(3*elastic_hf.k + elastic_hf.mu);
    
    % Print results
    fprintf('\n%s Hydrate Foam Properties:\n', foam_type);
    fprintf('Bulk Modulus (K): %.2f GPa\n', elastic_hf.k);
    fprintf('Shear Modulus (μ): %.2f GPa\n', elastic_hf.mu);
    fprintf('Young''s Modulus (E): %.2f GPa\n', elastic_hf.E);
    fprintf('Poisson''s Ratio (ν): %.3f\n\n', elastic_hf.nu);
end

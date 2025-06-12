function [Chom_cp, EEinfty_cp, elastic_cp] = homogenize_cement_paste(fractions, Chom_hf, stiffness, I, J, K, paste_type)
% HOMOGENIZE_CEMENT_PASTE Computes homogenized properties of cement paste using
% Mori-Tanaka scheme
%
% Inputs:
%   fractions  - Structure containing volume fractions
%   Chom_hf    - Homogenized stiffness of hydrate foam
%   stiffness  - Structure containing stiffness matrices
%   I, J, K    - Fundamental tensors
%   paste_type - Type of cement paste ('normal' or 'ITZ')
%
% Outputs:
%   Chom_cp    - Homogenized stiffness tensor of cement paste
%   EEinfty_cp - Strain concentration tensor
%   elastic_cp - Structure containing elastic properties (K, G, E, nu)

    fprintf('\n=== Cement Paste Homogenization (%s Cement Paste) ===\n', paste_type);
    
    % Select paste type
    switch paste_type
        case 'normal'
            paste_composition = fractions.cement_paste_composition;
            hydrate_foam_total = fractions.hydrate_foam.total_in_paste;
        case 'ITZ'
            paste_composition = fractions.ITZ_composition;
            hydrate_foam_total = fractions.ITZ_hydrate_foam.total_hydrate_foam;
        otherwise
            error('Invalid paste_type. Use ''normal'' or ''ITZ''.');
    end

    % Use homogenized hydrate foam properties as matrix
    C0 = Chom_hf;
    
    % Calculate Hill tensor for spherical clinker inclusions
    P_sph = fun_P_sphere_iso(C0);
    
    % Calculate strain concentration tensors
    Ainf_clin = inv(I + P_sph*(stiffness.clinker - C0));
    Ainf_hf = I;  % Identity tensor for matrix phase
    
    % Calculate overall strain concentration tensor
    EEinfty_cp = inv(paste_composition.clinker*Ainf_clin + ...
                     hydrate_foam_total*Ainf_hf);
    
    % Calculate phase-specific strain concentration tensors
    A_clin = Ainf_clin * EEinfty_cp;
    A_hf = Ainf_hf * EEinfty_cp;
    
    % Calculate homogenized stiffness
    Chom_cp = paste_composition.clinker*stiffness.clinker*A_clin + ...
              hydrate_foam_total*Chom_hf*A_hf;
    
    % Calculate elastic properties
    elastic_cp.k = Chom_cp(1,1) - 4/3*(Chom_cp(6,6)/2);
    elastic_cp.mu = Chom_cp(6,6)/2;
    elastic_cp.nu = (3*elastic_cp.k - 2*elastic_cp.mu)/(6*elastic_cp.k + 2*elastic_cp.mu);
    elastic_cp.E = 9*elastic_cp.k*elastic_cp.mu/(3*elastic_cp.k + elastic_cp.mu);
    
    % Print results
    fprintf('\n%s Cement Paste Properties:\n', paste_type);
    fprintf('Bulk Modulus (K): %.2f GPa\n', elastic_cp.k);
    fprintf('Shear Modulus (μ): %.2f GPa\n', elastic_cp.mu);
    fprintf('Young''s Modulus (E): %.2f GPa\n', elastic_cp.E);
    fprintf('Poisson''s Ratio (ν): %.3f\n\n', elastic_cp.nu);
end

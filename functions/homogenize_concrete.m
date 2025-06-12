function [Chom_conc, EEinfty_conc, elastic_conc] = homogenize_concrete(fractions, Chom_cp, stiffness, I, J, K)
% HOMOGENIZE_CONCRETE Computes homogenized properties of concrete using
% Mori-Tanaka scheme for uncoated aggregates and Herve-Zaoui solution for coated ones

    fprintf('\n=== Concrete Homogenization (Mori-Tanaka) ===\n');
    
    % Use homogenized cement paste properties as matrix
    C0 = Chom_cp;
    
    % Calculate Hill tensor for spherical inclusions
    P_sph = fun_P_sphere_iso(C0);
    
    % Initialize strain concentration tensors sum
    sum_f_A = fractions.cement_paste*I;  % Matrix phase
    
    % Handle uncoated aggregates
    uncoated_names = fieldnames(fractions.uncoated_aggregates);
    for i = 1:length(uncoated_names)
        agg_name = uncoated_names{i};
        % Calculate strain concentration tensor
        Ainf_uncoated = inv(I + P_sph*(stiffness.uncoated_aggregates.(agg_name) - C0));
        % Add contribution
        sum_f_A = sum_f_A + fractions.uncoated_aggregates.(agg_name)*Ainf_uncoated;
    end
    
    % Store Herve-Zaoui solutions for coated aggregates
    hz_solutions = struct();
    
    % Handle coated aggregates if present
    if isfield(fractions, 'coated_aggregates')
        coated_names = fieldnames(fractions.coated_aggregates);
        for i = 1:length(coated_names)
            agg_name = coated_names{i};
            
            % Calculate radii from volume fractions
            R1 = (1/2)*6^(1/3)*fractions.coated_aggregates.(agg_name).core^(1/3)/pi^(1/3);
            R2 = (1/2)*6^(1/3)*(fractions.coated_aggregates.(agg_name).core + ...
                  fractions.ITZ)^(1/3)/pi^(1/3);
            R3 = (1/2)*6^(1/3)*(fractions.coated_aggregates.(agg_name).core + ...
                  fractions.ITZ + fractions.cement_paste)^(1/3)/pi^(1/3);
            
            % Extract elastic properties for Herve-Zaoui solution
            k_core = stiffness.coated_aggregates.(agg_name)(1,1) - 4/3*(stiffness.coated_aggregates.(agg_name)(6,6)/2);
            mu_core = stiffness.coated_aggregates.(agg_name)(6,6)/2;
            
            % Calculate bulk cement paste properties from Chom_cp
            k_cp = Chom_cp(1,1) - 4/3*(Chom_cp(6,6)/2);
            mu_cp = Chom_cp(6,6)/2;
            E_cp = 9*k_cp*mu_cp/(3*k_cp + mu_cp);
            nu_cp = (3*k_cp - 2*mu_cp)/(2*(3*k_cp + mu_cp));
            
            % Calculate ITZ properties
            ITZ_reduction = 0.78;  % 20% reduction in E
            E_ITZ = E_cp * ITZ_reduction;
            nu_ITZ = nu_cp;  % Same Poisson's ratio as bulk paste
            
            % Calculate k_ITZ and mu_ITZ using fun_kmu_from_Enu
            [k_ITZ, mu_ITZ] = fun_kmu_from_Enu(E_ITZ, nu_ITZ);
            
            % Apply Herve-Zaoui solution
            input_matrix = [k_core, mu_core, R1;      % Core aggregate
                           k_ITZ, mu_ITZ, R2;         % ITZ layer
                           k_cp, mu_cp, R3];          % Bulk cement paste
                        
            % Store solution for later use
            hz_solutions.(agg_name) = fun_HZ_int(input_matrix);
            
            % Calculate strain concentration tensors
            Ainf_core = hz_solutions.(agg_name)(1,5)*J + hz_solutions.(agg_name)(1,7)*K;
            Ainf_coating = hz_solutions.(agg_name)(2,5)*J + hz_solutions.(agg_name)(2,7)*K;
            
            % Add contributions
            sum_f_A = sum_f_A + fractions.coated_aggregates.(agg_name).core*Ainf_core + ...
                     fractions.coated_aggregates.(agg_name).coating*Ainf_coating;
        end
    end
    
    % Add condition number check before matrix inversion
    cond_num = cond(sum_f_A);
    if cond_num > 1e15
        warning('Poor conditioning in strain concentration tensor calculation. Condition number: %.2e', cond_num);
    end
    
    % Calculate overall strain concentration tensor
    EEinfty_conc = inv(sum_f_A);
    
    % Check for NaN values
    if any(isnan(EEinfty_conc(:)))
        warning('NaN values detected in strain concentration tensor');
        disp('sum_f_A matrix:');
        disp(sum_f_A);
    end
    
    
    % Calculate homogenized stiffness
    Chom_conc = fractions.cement_paste*Chom_cp*EEinfty_conc;
    
    % Add uncoated aggregate contributions
    for i = 1:length(uncoated_names)
        agg_name = uncoated_names{i};
        Ainf_uncoated = inv(I + P_sph*(stiffness.uncoated_aggregates.(agg_name) - C0));
        A_uncoated = Ainf_uncoated*EEinfty_conc;
        Chom_conc = Chom_conc + ...
            fractions.uncoated_aggregates.(agg_name)*stiffness.uncoated_aggregates.(agg_name)*A_uncoated;
    end
    
    % Add coated aggregate contributions
    if isfield(fractions, 'coated_aggregates')
        for i = 1:length(coated_names)
            agg_name = coated_names{i};
            
            % Use stored Herve-Zaoui solution
            A_core = (hz_solutions.(agg_name)(1,5)*J + hz_solutions.(agg_name)(1,7)*K)*EEinfty_conc;
            A_coating = (hz_solutions.(agg_name)(2,5)*J + hz_solutions.(agg_name)(2,7)*K)*EEinfty_conc;
            
            % Add contributions
            Chom_conc = Chom_conc + ...
                fractions.coated_aggregates.(agg_name).core*stiffness.coated_aggregates.(agg_name)*A_core + ...
                fractions.coated_aggregates.(agg_name).coating*Chom_cp*A_coating;
        end
    end
    
    % Calculate elastic properties
    elastic_conc.k = Chom_conc(1,1) - 4/3*(Chom_conc(6,6)/2);
    elastic_conc.mu = Chom_conc(6,6)/2;
    elastic_conc.nu = (3*elastic_conc.k - 2*elastic_conc.mu)/(6*elastic_conc.k + 2*elastic_conc.mu);
    elastic_conc.E = 9*elastic_conc.k*elastic_conc.mu/(3*elastic_conc.k + elastic_conc.mu);
    
    % Print results
    fprintf('\nConcrete Properties:\n');
    fprintf('Bulk Modulus (K): %.2f GPa\n', elastic_conc.k);
    fprintf('Shear Modulus (μ): %.2f GPa\n', elastic_conc.mu);
    fprintf('Young''s Modulus (E): %.2f GPa\n', elastic_conc.E);
    fprintf('Poisson''s Ratio (ν): %.3f\n\n', elastic_conc.nu);
end
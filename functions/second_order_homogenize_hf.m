function [Chom_hf_strength] = second_order_homogenize_hf(fractions, stiffness, materialProperties, Chom_hf, analysisParams, I, J, K)
    % SECOND_ORDER_HOMOGENIZE_HF 
    % Computes second-order homogenization of hydrate foam for strength prediction
    % using the already computed Chom_hf as reference stiffness
    
    %% 1. Initialize parameters and extract properties
    % Extract volume fractions
    f_hyd = fractions.hydrate_foam.hydrates;
    f_por = fractions.hydrate_foam.porosity;
    
    % Extract stiffness matrices
    C_hyd = stiffness.hydrate;
    C_por = stiffness.pore;
    
    % Extract elastic parameters
    k_hyd = materialProperties.hydrate.k;
    mu_hyd = materialProperties.hydrate.mu;
    
    % Analysis parameters
    tolerance_2 = analysisParams.tolerance_2;
    fapp = 0.001;  % Perturbation parameter
    f_hyd_diff = fapp/5;

    % Use Stroud points from parameters
    stroud_azi = analysisParams.stroud_azi;
    stroud_zeni = analysisParams.stroud_zeni;
    
    % Iteration parameters
    deviation = 1;
    countermark = 0;
    iteration = 8;
    
    % % Print initial information
    % fprintf('\nStarting Second-Order Homogenization for Hydrate Foam\n');
    % fprintf('Volume Fractions: Hydrates = %.3f, Pores = %.3f\n', f_hyd, f_por);
    % fprintf('Perturbation Parameters: fapp = %.3f, f_hyd_diff = %.3f\n', fapp, f_hyd_diff);
    % 
    %% 2. Second order deviatoric strength
    % Use Chom_hf as reference stiffness
    C0 = Chom_hf;
    
    % Perturbed hydrate properties for deviatoric part
    mu_hyd_diff = (1 + fapp) * mu_hyd;
    C_hyd_diff = 3*k_hyd*J + 2*mu_hyd_diff*K;
    
    % % Print deviatoric calculation header
    % fprintf('\n--- Deviatoric Strength Calculation ---\n');
    % fprintf('Iteration   Error       C1111        C3333\n');
    % fprintf('----------------------------------------\n');
    
    % Iterative solution for deviatoric part
    while deviation > tolerance_2 && countermark < iteration
        % fprintf('%4d      %10.2e  %10.2e  %10.2e\n', ...
        %     countermark, deviation, C0(1,1), C0(3,3));
        countermark = countermark + 1;
        
        % Spherical pores
        P_p = P_isotrans_sph(C0);
        Ainf_p = inv(I + P_p*(C_por - C0));
        
        % Acicular hydrates, orientated isotropically in all space directions
        sumAinf_hyd = zeros(6,6);
        for i = 1:15  % Using the 15-point Stroud integration
            % Transformation matrices using pre-defined Stroud points
            azi = stroud_azi(i);   % Using the pre-defined stroud_azi
            zeni = stroud_zeni(i); % Using the pre-defined stroud_zeni
            Q4 = fun_Q4_bp(azi, zeni); 
            Q4t = transpose(Q4);
            
            % Transform stiffness to needle orientation
            C0_aniso = Q4*C0*Q4t;
            P_hyd_i_e3 = fun_P_ellipsoid_transiso(C0_aniso, 1, 1e-20);
            % Transform back to global coordinates
            P_hyd_i = Q4t * P_hyd_i_e3 * Q4;
            
            % Calculate and sum concentration tensors
            Ainf_hyd_i = inv(I + P_hyd_i*(C_hyd - C0));
            sumAinf_hyd = sumAinf_hyd + Ainf_hyd_i;
        end
        sumAinf_hyd = 1/15 * sumAinf_hyd;
        
        % Acicular hydrate with larger mu, oriented in e3
        P_hyd_e3 = fun_P_ellipsoid_transiso(C0, 1, 1e-20);
        Ainf_hyd_diff = inv(I + P_hyd_e3*(C_hyd_diff - C0));
        Ainf_hyd_minus = inv(I + P_hyd_e3*(C_hyd - C0));
        
        % Strain concentration tensors
        EEinfty_hf = inv(f_por*Ainf_p + f_hyd*sumAinf_hyd + ...
                        f_hyd_diff*Ainf_hyd_diff - f_hyd_diff*Ainf_hyd_minus);
        
        A_p = Ainf_p*EEinfty_hf;
        A_hyd = sumAinf_hyd*EEinfty_hf;
        A_hyd_diff = Ainf_hyd_diff*EEinfty_hf;
        A_hyd_minus = Ainf_hyd_minus*EEinfty_hf;
        
        % Homogenized stiffness - SELF CONSISTENT
        Chom_hf_dev_transiso = f_hyd*C_hyd*A_hyd + ...
                              f_hyd_diff*C_hyd_diff*A_hyd_diff - ...
                              f_hyd_diff*C_hyd*A_hyd_minus;
        
        % Update and check convergence
        C0_old = C0;
        C0 = Chom_hf_dev_transiso;
        deviation = abs(norm(C0 - C0_old)/norm(C0));
    end
    
    % % Print deviatoric results
    % fprintf('Deviatoric calculation completed in %d iterations\n', countermark);
    % if deviation <= tolerance_2
    %     fprintf('Converged to tolerance of %.2e\n', tolerance_2);
    % else
    %     fprintf('Warning: Maximum iterations reached before convergence\n');
    % end
    
    %% 3. Second order volumetric strength
    % Reset reference stiffness to Chom_hf
    C0 = Chom_hf;
    
    % Perturbed hydrate properties for volumetric part
    k_hyd_diff = (1 + fapp) * k_hyd;
    C_hyd_diff = 3*k_hyd_diff*J + 2*mu_hyd*K;
    
    % Reset iteration parameters
    deviation = 1;
    countermark = 0;
    
    % % Print volumetric calculation header
    % fprintf('\n--- Volumetric Strength Calculation ---\n');
    % fprintf('Iteration   Error       C1111        C3333\n');
    % fprintf('----------------------------------------\n');
    
    % Iterative solution for volumetric part
    while deviation > tolerance_2 && countermark < iteration
        % fprintf('%4d      %10.2e  %10.2e  %10.2e\n', ...
        %     countermark, deviation, C0(1,1), C0(3,3));
        countermark = countermark + 1;
        
        % Spherical pores
        P_p = P_isotrans_sph(C0);
        Ainf_p = inv(I + P_p*(C_por - C0));
        
        % Acicular hydrates, orientated isotropically in all space directions
        sumAinf_hyd = zeros(6,6);
        for i = 1:15  % Using the 15-point Stroud integration
            % Transformation matrices using pre-defined Stroud points
            azi = stroud_azi(i);   % Using the pre-defined stroud_azi
            zeni = stroud_zeni(i); % Using the pre-defined stroud_zeni
            Q4 = fun_Q4_bp(azi, zeni); 
            Q4t = transpose(Q4);
            
            % Transform stiffness
            C0_aniso = Q4*C0*Q4t;
            P_hyd_i_e3 = fun_P_ellipsoid_transiso(C0_aniso, 1, 1e-20);
            P_hyd_i = Q4t * P_hyd_i_e3 * Q4;
            
            % Calculate and sum concentration tensors
            Ainf_hyd_i = inv(I + P_hyd_i*(C_hyd - C0));
            sumAinf_hyd = sumAinf_hyd + Ainf_hyd_i;
        end
        sumAinf_hyd = 1/15 * sumAinf_hyd;
        
        % Acicular hydrate with larger k, oriented in e3
        P_hyd_e3 = fun_P_ellipsoid_transiso(C0, 1, 1e-20);
        Ainf_hyd_diff = inv(I + P_hyd_e3*(C_hyd_diff - C0));
        Ainf_hyd_minus = inv(I + P_hyd_e3*(C_hyd - C0));
        
        % Strain concentration tensors
        EEinfty_hf = inv(f_por*Ainf_p + f_hyd*sumAinf_hyd + ...
                        f_hyd_diff*Ainf_hyd_diff - f_hyd_diff*Ainf_hyd_minus);
        
        A_p = Ainf_p*EEinfty_hf;
        A_hyd = sumAinf_hyd*EEinfty_hf;
        A_hyd_diff = Ainf_hyd_diff*EEinfty_hf;
        A_hyd_minus = Ainf_hyd_minus*EEinfty_hf;
        
        % Homogenized stiffness - SELF CONSISTENT
        Chom_hf_vol_transiso = f_hyd*C_hyd*A_hyd + ...
                              f_hyd_diff*C_hyd_diff*A_hyd_diff - ...
                              f_hyd_diff*C_hyd*A_hyd_minus;
        
        % Update and check convergence
        C0_old = C0;
        C0 = Chom_hf_vol_transiso;
        deviation = abs(norm(C0 - C0_old)/norm(C0));
    end
    
    % % Print volumetric results
    % fprintf('Volumetric calculation completed in %d iterations\n', countermark);
    % if deviation <= tolerance_2
    %     fprintf('Converged to tolerance of %.2e\n', tolerance_2);
    % else
    %     fprintf('Warning: Maximum iterations reached before convergence\n');
    % end
    % 
    %% 4. Package results and print summary
    Chom_hf_strength = struct(...
        'volumetric', Chom_hf_vol_transiso, ...
        'deviatoric', Chom_hf_dev_transiso);
    
    % % Print final results summary
    % fprintf('\nFinal Results Summary:\n');
    % fprintf('Deviatoric C1111: %.4e, C3333: %.4e\n', ...
    %     Chom_hf_strength.deviatoric(1,1), Chom_hf_strength.deviatoric(3,3));
    % fprintf('Volumetric C1111: %.4e, C3333: %.4e\n', ...
    %     Chom_hf_strength.volumetric(1,1), Chom_hf_strength.volumetric(3,3));
    % fprintf('----------------------------------------\n\n');
    
   end
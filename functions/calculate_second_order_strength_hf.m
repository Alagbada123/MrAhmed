function [Chom_hf_dev, Chom_hf_vol] = calculate_second_order_strength_hf(fractions, stiffness, materialProperties, Chom_hf, analysisParams, I, J, K)
% CALCULATE_SECOND_ORDER_STRENGTH_HF Computes second-order strength tensors for hydrate foam
%
% Inputs:
%   fractions          - Structure containing volume fractions
%   stiffness          - Structure containing stiffness matrices
%   materialProperties - Structure containing material properties
%   Chom_hf           - Homogenized stiffness of hydrate foam
%   analysisParams    - Analysis parameters including tolerances and Stroud points
%   I, J, K           - Fundamental tensors
%
% Outputs:
%   Chom_hf_dev - Deviatoric strength tensor
%   Chom_hf_vol - Volumetric strength tensor

    fprintf('\n=== Second-Order Strength Calculation (Hydrate Foam) ===\n');
    
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
    
    %% 2. Calculate perturbed properties
    % Deviatoric perturbation
    mu_hyd_diff = (1 + fapp) * mu_hyd;
    C_hyd_diff_dev = 3*k_hyd*J + 2*mu_hyd_diff*K;
    
    % Volumetric perturbation
    k_hyd_diff = (1 + fapp) * k_hyd;
    C_hyd_diff_vol = 3*k_hyd_diff*J + 2*mu_hyd*K;

    %% 3. Deviatoric part calculation
    fprintf('\nCalculating deviatoric part:\n');
    C0 = Chom_hf;
    deviation = 1;
    iteration = 0;
    
    while deviation > tolerance_2 && iteration < 8
        iteration = iteration + 1;
        
        % Calculate Hill tensor for spherical pores
        P_pore = P_isotrans_sph(C0);
        Ainf_pore = inv(I + P_pore*(C_por - C0));
        
        % Calculate strain concentration tensors for hydrates
        sumAinf_hyd = zeros(6,6);
        for i = 1:15
            % Transformation matrices
            azi = stroud_azi(i); 
            zeni = stroud_zeni(i);
            Q4 = fun_Q4_bp(azi, zeni); 
            Q4t = transpose(Q4);
            
            % Transform stiffness to needle orientation
            C0_aniso = Q4*C0*Q4t;
            P_hyd_i_e3 = fun_P_ellipsoid_transiso(C0_aniso, 1, 1e-20);
            P_hyd_i = Q4t * P_hyd_i_e3 * Q4;
            
            % Calculate and sum strain concentration tensors
            Ainf_hyd_i = inv(I + P_hyd_i*(C_hyd - C0));
            sumAinf_hyd = sumAinf_hyd + Ainf_hyd_i;
        end
        sumAinf_hyd = sumAinf_hyd/15;
        
        % Calculate strain concentration tensor for perturbed hydrate
        P_hyd_e3 = fun_P_ellipsoid_transiso(C0, 1, 1e-20);
        Ainf_hyd_diff = inv(I + P_hyd_e3*(C_hyd_diff_dev - C0));
        Ainf_hyd_minus = inv(I + P_hyd_e3*(C_hyd - C0));
        
        % Calculate overall strain concentration tensor
        EEinfty_hf = inv(f_por*Ainf_pore + f_hyd*sumAinf_hyd + ...
                        f_hyd_diff*Ainf_hyd_diff - f_hyd_diff*Ainf_hyd_minus);
        
        % Calculate phase-specific strain concentration tensors
        A_hyd = sumAinf_hyd*EEinfty_hf;
        A_hyd_diff = Ainf_hyd_diff*EEinfty_hf;
        A_hyd_minus = Ainf_hyd_minus*EEinfty_hf;
        
        % Calculate homogenized stiffness
        Chom_hf_dev_transiso = f_hyd*C_hyd*A_hyd + ...
                              f_hyd_diff*C_hyd_diff_dev*A_hyd_diff - ...
                              f_hyd_diff*C_hyd*A_hyd_minus;
        
        % Update
        C0_old = C0;
        C0 = Chom_hf_dev_transiso;
        deviation = abs(norm(C0-C0_old)/norm(C0));
        fprintf('Iteration %d (Deviatoric): Convergence = %.2e\n', iteration, deviation);
    end
    Chom_hf_dev = Chom_hf_dev_transiso;
    
    %% 4. Volumetric part calculation
    fprintf('\nCalculating volumetric part:\n');
    C0 = Chom_hf;
    deviation = 1;
    iteration = 0;
    
    while deviation > tolerance_2 && iteration < 8
        iteration = iteration + 1;
        
        % Calculate Hill tensor for spherical pores
        P_pore = P_isotrans_sph(C0);
        Ainf_pore = inv(I + P_pore*(C_por - C0));
        
        % Calculate strain concentration tensors for hydrates
        sumAinf_hyd = zeros(6,6);
        for i = 1:15
            % Transformation matrices
            azi = stroud_azi(i); 
            zeni = stroud_zeni(i);
            Q4 = fun_Q4_bp(azi, zeni); 
            Q4t = transpose(Q4);
            
            % Transform stiffness to needle orientation
            C0_aniso = Q4*C0*Q4t;
            P_hyd_i_e3 = fun_P_ellipsoid_transiso(C0_aniso, 1, 1e-20);
            P_hyd_i = Q4t * P_hyd_i_e3 * Q4;
            
            % Calculate and sum strain concentration tensors
            Ainf_hyd_i = inv(I + P_hyd_i*(C_hyd - C0));
            sumAinf_hyd = sumAinf_hyd + Ainf_hyd_i;
        end
        sumAinf_hyd = sumAinf_hyd/15;
        
        % Calculate strain concentration tensor for perturbed hydrate
        P_hyd_e3 = fun_P_ellipsoid_transiso(C0, 1, 1e-20);
        Ainf_hyd_diff = inv(I + P_hyd_e3*(C_hyd_diff_vol - C0));
        Ainf_hyd_minus = inv(I + P_hyd_e3*(C_hyd - C0));
        
        % Calculate overall strain concentration tensor
        EEinfty_hf = inv(f_por*Ainf_pore + f_hyd*sumAinf_hyd + ...
                        f_hyd_diff*Ainf_hyd_diff - f_hyd_diff*Ainf_hyd_minus);
        
        % Calculate phase-specific strain concentration tensors
        A_hyd = sumAinf_hyd*EEinfty_hf;
        A_hyd_diff = Ainf_hyd_diff*EEinfty_hf;
        A_hyd_minus = Ainf_hyd_minus*EEinfty_hf;
        
        % Calculate homogenized stiffness
        Chom_hf_vol_transiso = f_hyd*C_hyd*A_hyd + ...
                              f_hyd_diff*C_hyd_diff_vol*A_hyd_diff - ...
                              f_hyd_diff*C_hyd*A_hyd_minus;
        
        % Update
        C0_old = C0;
        C0 = Chom_hf_vol_transiso;
        deviation = abs(norm(C0-C0_old)/norm(C0));
        fprintf('Iteration %d (Volumetric): Convergence = %.2e\n', iteration, deviation);
    end
    Chom_hf_vol = Chom_hf_vol_transiso;
end
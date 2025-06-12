function [Chom_cp_strength, diffQ_vol_cp, diffQ_dev_cp] = second_order_homogenize_cp(fractions, stiffness, materialProperties, Chom_cp, Chom_hf_strength, analysisParams, I, J, K)
    % SECOND_ORDER_HOMOGENIZE_CP 
    % Computes second-order homogenization of cement paste for strength prediction
    % using Chom_hf_strength from second_order_homogenize_hf
    
    %% 1. Initialize parameters and extract properties
    % Extract volume fractions
    f_cem = fractions.cement_paste_composition.clinker;
    f_hf = fractions.hydrate_foam.total_in_paste;
    
    % Extract stiffness matrices
    C_cem = stiffness.clinker;
    
    % Extract parameters for difference quotients
    fapp = 0.001;  % Perturbation parameter
    k_hyd = materialProperties.hydrate.k;
    mu_hyd = materialProperties.hydrate.mu;
    
    %% 2. Second order deviatoric strength
    % Use deviatoric part from hydrate foam strength
    C0 = Chom_hf_strength.deviatoric;
    
    % spherical clinkers, SCMs, and inert fillers
    P_sph = P_isotrans_sph(C0);
    Ainf_cem = inv(I + P_sph*(C_cem - C0));
    Ainf_hf = I;
    
    % strain concentration tensors
    EEinfty_cp_dev_transiso = inv(f_cem*Ainf_cem + f_hf*Ainf_hf);
    A_cem_dev_transiso = Ainf_cem*EEinfty_cp_dev_transiso;
    A_hf_dev_transiso = Ainf_hf*EEinfty_cp_dev_transiso;
    
    % Homogenized stiffness (MORI-TANAKA)
    Chom_cp_dev_transiso = f_cem*C_cem*A_cem_dev_transiso + f_hf*Chom_hf_strength.deviatoric*A_hf_dev_transiso;
    
    %% 3. Second order volumetric strength
    % Use volumetric part from hydrate foam strength
    C0 = Chom_hf_strength.volumetric;
    
    % spherical clinkers, SCMs, and inert fillers
    P_sph = P_isotrans_sph(C0);
    Ainf_cem = inv(I + P_sph*(C_cem - C0));
    Ainf_hf = I;
    
    % strain concentration tensors
    EEinfty_cp_vol_transiso = inv(f_cem*Ainf_cem + f_hf*Ainf_hf);
    A_cem_vol_transiso = Ainf_cem*EEinfty_cp_vol_transiso;
    A_hf_vol_transiso = Ainf_hf*EEinfty_cp_vol_transiso;
    
    % Homogenized stiffness (MORI-TANAKA)
    Chom_cp_vol_transiso = f_cem*C_cem*A_cem_vol_transiso + f_hf*Chom_hf_strength.volumetric*A_hf_vol_transiso;
    
    %% 4. Package results
    Chom_cp_strength = struct(...
        'volumetric', Chom_cp_vol_transiso, ...
        'deviatoric', Chom_cp_dev_transiso);
    
    %% 5. Calculate difference quotients
    diffQ_vol_cp = 1/(fapp*k_hyd) * (Chom_cp_vol_transiso - Chom_cp);
    diffQ_dev_cp = 1/(fapp*mu_hyd) * (Chom_cp_dev_transiso - Chom_cp);
end
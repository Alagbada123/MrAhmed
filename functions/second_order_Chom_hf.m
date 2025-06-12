function Chom_hf_deriv = second_order_Chom_hf(fractions, stiffness, materialProperties, Chom_hf, analysisParams, I, J, K, mode)
% SECOND_ORDER_Chom_hf
% Performs second-order expansions for the hydrate foam stiffness.
% 
% INPUTS:
%   fractions         - Struct with hydrate foam volume fractions
%   stiffness         - Struct with hydrate and pore stiffness matrices
%   materialProperties- Struct with hydrate properties (k, mu)
%   Chom_hf           - First-order homogenized stiffness (6x6 matrix)
%   analysisParams    - Struct with analysis parameters (e.g., tolerance_2)
%   I, J, K           - Identity and projection tensors (6x6 matrices)
%   mode              - 'normal' or 'ITZ', specifies which fractions to use
% 
% OUTPUT:
%   Chom_hf_deriv     - Struct with deviatoric, volumetric corrections, and Ainf_hyd_e3

% Extract fractions based on mode
if strcmp(mode, 'normal')
    f_hyd = fractions.hydrate_foam.hydrates;
    f_por = fractions.hydrate_foam.porosity;
elseif strcmp(mode, 'ITZ')
    f_hyd = fractions.ITZ_hydrate_foam.hydrates;
    f_por = fractions.ITZ_hydrate_foam.porosity;
else
    error('Invalid mode specified. Use "normal" or "ITZ".');
end

% Extract stiffness and properties
C_hyd = stiffness.hydrate;
C_por = stiffness.pore;
k_hyd = materialProperties.hydrate.k;
mu_hyd = materialProperties.hydrate.mu;

% Extract parameters
stroud_azi = analysisParams.stroud_azi;
stroud_zeni = analysisParams.stroud_zeni;
tolerance_2 = analysisParams.tolerance_2;

% Perturbation factor
fapp = 0.001;
f_hyd_diff = fapp / 5;

%% Deviatoric Expansion
C0_dev = Chom_hf;
mu_hyd_diff = (1 + fapp) * mu_hyd;
C_hyd_diff_d = 3 * k_hyd * J + 2 * mu_hyd_diff * K;
deviation = 1;
iterCount = 0;
maxIter = 8;

while (deviation > tolerance_2) && (iterCount < maxIter)
    iterCount = iterCount + 1;

    % Spherical pores
    P_p = P_isotrans_sph(C0_dev);
    Ainf_p = inv(I + P_p * (C_por - C0_dev));

    % Sum over 15 directions for acicular hydrates
    sumAinf_hyd = zeros(6, 6);
    nDir = 15;
    for i = 1:nDir
        azi = stroud_azi(i);
        zeni = stroud_zeni(i);
        Q4 = fun_Q4_bp(azi, zeni);
        Q4t = transpose(Q4);

        % Transform to local coordinate
        C0_aniso = Q4 * C0_dev * Q4t;
        P_hyd_i_e3 = fun_P_ellipsoid_transiso(C0_aniso, 1, 1e-20);
        % Transform back to global coordinate
        P_hyd_i = Q4t * P_hyd_i_e3 * Q4;

        % Strain concentration tensors
        Ainf_hyd_i = inv(I + P_hyd_i * (C_hyd - C0_dev));
        sumAinf_hyd = sumAinf_hyd + Ainf_hyd_i;
    end
    sumAinf_hyd = sumAinf_hyd / nDir;

    % Single needle with bigger shear, oriented in e3
    P_hyd_e3_d = fun_P_ellipsoid_transiso(C0_dev, 1, 1e-20);
    Ainf_hyd_diff = inv(I + P_hyd_e3_d * (C_hyd_diff_d - C0_dev));
    Ainf_hyd_minus_d = inv(I + P_hyd_e3_d * (C_hyd - C0_dev));

    % Total strain concentration
    EEinfty_hf_dev = inv(f_por * Ainf_p + f_hyd * sumAinf_hyd + ...
                         f_hyd_diff * Ainf_hyd_diff - f_hyd_diff * Ainf_hyd_minus_d);

    % Homogenized stiffness
    Chom_hf_dev_transiso = f_hyd * C_hyd * (sumAinf_hyd * EEinfty_hf_dev) + ...
                           f_hyd_diff * C_hyd_diff_d * (Ainf_hyd_diff * EEinfty_hf_dev) - ...
                           f_hyd_diff * C_hyd * (Ainf_hyd_minus_d * EEinfty_hf_dev);

    % Update and convergence
    deviation = abs(norm(Chom_hf_dev_transiso - C0_dev) / max(norm(C0_dev), eps));
    C0_dev = Chom_hf_dev_transiso;
end

%% Volumetric Expansion
C0_vol = Chom_hf;
k_hyd_diff = (1 + fapp) * k_hyd;
C_hyd_diff_v = 3 * k_hyd_diff * J + 2 * mu_hyd * K;
deviation = 1;
iterCount = 0;

while (deviation > tolerance_2) && (iterCount < maxIter)
    iterCount = iterCount + 1;

    % Spherical pores
    P_p = P_isotrans_sph(C0_vol);
    Ainf_p = inv(I + P_p * (C_por - C0_vol));

    % Sum over 15 directions for acicular hydrates
    sumAinf_hyd = zeros(6, 6);
    for i = 1:nDir
        azi = stroud_azi(i);
        zeni = stroud_zeni(i);
        Q4 = fun_Q4_bp(azi, zeni);
        Q4t = transpose(Q4);

        % Transform to local coordinate
        C0_aniso = Q4 * C0_vol * Q4t;
        P_hyd_i_e3 = fun_P_ellipsoid_transiso(C0_aniso, 1, 1e-20);
        % Transform back to global coordinate
        P_hyd_i = Q4t * P_hyd_i_e3 * Q4;

        % Strain concentration tensors
        Ainf_hyd_i = inv(I + P_hyd_i * (C_hyd - C0_vol));
        sumAinf_hyd = sumAinf_hyd + Ainf_hyd_i;
    end
    sumAinf_hyd = sumAinf_hyd / nDir;

    % Single needle with bigger bulk, oriented in e3
    P_hyd_e3_v = fun_P_ellipsoid_transiso(C0_vol, 1, 1e-20);
    Ainf_hyd_diff = inv(I + P_hyd_e3_v * (C_hyd_diff_v - C0_vol));
    Ainf_hyd_minus = inv(I + P_hyd_e3_v * (C_hyd - C0_vol));

    % Total strain concentration
    EEinfty_hf_vol = inv(f_por * Ainf_p + f_hyd * sumAinf_hyd + ...
                         f_hyd_diff * Ainf_hyd_diff - f_hyd_diff * Ainf_hyd_minus);

    % Homogenized stiffness
    Chom_hf_vol_transiso = f_hyd * C_hyd * (sumAinf_hyd * EEinfty_hf_vol) + ...
                           f_hyd_diff * C_hyd_diff_v * (Ainf_hyd_diff * EEinfty_hf_vol) - ...
                           f_hyd_diff * C_hyd * (Ainf_hyd_minus * EEinfty_hf_vol);

    % Update and convergence
    deviation = abs(norm(Chom_hf_vol_transiso - C0_vol) / max(norm(C0_vol), eps));
    C0_vol = Chom_hf_vol_transiso;
end

%% Hill Tensor for e3-Oriented Hydrate
P_hyd_e3 = fun_P_ellipsoid_transiso(Chom_hf, 1, 1e-20);
Ainf_hyd_e3 = inv(I + P_hyd_e3 * (C_hyd - Chom_hf));

%% Package Results
Chom_hf_deriv = struct('volumetric', Chom_hf_vol_transiso, 'deviatoric', Chom_hf_dev_transiso, 'Ainf_hyd_e3', Ainf_hyd_e3);
end

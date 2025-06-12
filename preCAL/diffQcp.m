%% PRECALCULATION of the difference Quotients for second order stress concentration to hydrates
clear all; clc; close all

% --- NEW: Add path for helper functions ---
% Ensure functions like fun_kmu_from_Enu, fun_CfromEnu, fun_Cfromkmu,
% fun_P_sphere_iso, fun_Ainf_needle_iso, P_isotrans_sph, fun_Q4_bp,
% fun_P_ellipsoid_transiso are in this folder or already on the MATLAB path.
addpath('../functions'); % Looks for functions one level up from preCAL;
% --- Assuming stroud_points script defines stroud_azi, stroud_zeni ---
stroud_points; % Keep this if needed for anisotropic calculations

%% Input Parameters (Updated Section)
% --- Grid to calculate on ---
wc_list = [0.50, 0.55];       % NEW water-cement ratios
xi_list = 0.05:0.05:1;         % NEW fixed hydration degree list
Fpor_list = 1.0:0.01:1.7 ;     % NEW porosity factor list

% --- Densities ---
rhoH2O = 1000;      % [kg/m^3] Water density
rhohyd = 2073;      % [kg/m^3] Hydrate density (Example value)
rhocem = 3150;      % [kg/m^3] Cement density (Example value)

% --- Phase Stiffness Properties ---
% Cement Clinker (assuming GPa for E)
Eclin = 139.9;      % [GPa] Young's Modulus
nuclin = 0.3;       % Poisson's ratio
[kclin, muclin] = fun_kmu_from_Enu(Eclin, nuclin); % Bulk and Shear Moduli
Cclin = fun_CfromEnu(Eclin, nuclin); % Stiffness Tensor (Voigt notation)

% Hydrates (assuming GPa for E)
Ehyd = 29.15786664; % [GPa] Young's Modulus
nuhyd = 0.24;       % Poisson's ratio
[khyd, muhyd] = fun_kmu_from_Enu(Ehyd, nuhyd); % Bulk and Shear Moduli
Chyd = fun_CfromEnu(Ehyd, nuhyd); % Stiffness Tensor (Voigt notation)

% Pores (filled with water, assuming GPa for k)
k_H2O = 2.3;        % [GPa] Bulk Modulus of Water
mu_H2O = 0;         % Shear Modulus of Water (fluid)
Cp = fun_Cfromkmu(k_H2O, mu_H2O); % Stiffness Tensor (Voigt notation)
% Note: Original code seemed to assume empty pores (Cp=0 implicitly in SC part)
%       If pores are empty, use: Cp = zeros(6,6);

% --- Homogenization Precision ---
tol_iso = 1e-10;     % Tolerance for isotropic homogenization
tol_aniso = 1e-8;    % Tolerance for anisotropic homogenization
maxcounter = 5;      % Max iterations for anisotropic homogenization

% --- Difference Quotient Calculation ---
fapp = 0.001;        % Small perturbation factor
khyddiff = (1 + fapp) * khyd; % Perturbed hydrate bulk modulus
muhyddiff = (1 + fapp) * muhyd; % Perturbed hydrate shear modulus
fhyddiff = fapp / 5;   % Small volume fraction for perturbation method

% --- Output ---
filename = 'precalc_cpITZOD4_updated.mat'; % Consider a new filename

%% Initialization (Recommended Structure - With Fpor)
% Define structure template for results
calc_struct=struct('E',NaN,'nu',NaN,'diffQvol',NaN,'diffQdev',NaN);
% Define main precalc structure template
precalc_struct=struct('wc',NaN,'xi',NaN,'Fpor',NaN, ...   % Inputs (Added Fpor)
                      'calc_cp',calc_struct, ...         % Results for cement paste
                      'calc_hf',struct('E',NaN,'nu',NaN), ... % Simplified results for hydrate foam
                      'precision',NaN(1,3), ...         % Tolerances, maxiter
                      'numerics',NaN(1,2), ...          % fapp, fhyddiff
                      'vol_cp',struct('cem',NaN,'hf',NaN),... % Vol frac WITHIN cement paste
                      'vol_hf',struct('hyd',NaN,'por',NaN)); % Vol frac WITHIN hydrate foam


% --- Load existing data or initialize output structure ---
if exist(filename, 'file') == 2
    load(filename)
    disp(['Warning: ',filename,' exists and is loaded. Appending new results if any.'])
    % Basic check if size matches new lists - might need more robust checking
    current_size = size(outputITZ_cell);
    if current_size(1)~=length(wc_list) || current_size(2)~=length(xi_list) || current_size(3)~=length(Fpor_list)
       disp('Warning: Loaded data size does not match current parameter lists.')
       % Decide how to handle mismatch: error, re-initialize, or attempt merge (complex)
       % For simplicity, we'll re-initialize here if sizes mismatch drastically
       if prod(current_size) == 0 || ... % If loaded is empty
          current_size(1) > length(wc_list) || ... % Or significantly different
          current_size(2) > length(xi_list) || ...
          current_size(3) > length(Fpor_list)
           disp('Re-initializing output structure due to size mismatch.');
           outputITZ_cell=repmat({precalc_struct}, length(wc_list),length(xi_list),length(Fpor_list));
       else
           disp('Attempting to use loaded data. Ensure consistency.');
           % Optional: Add code here to resize/pad 'outputITZ_cell' if needed
       end
    end
else
    disp([filename,' does not exist. Initializing new output structure.'])
    outputITZ_cell=repmat({precalc_struct}, length(wc_list),length(xi_list),length(Fpor_list));
end

% --- Define standard identity tensors (if not defined elsewhere) ---
% It's good practice to define these explicitly if not loaded
% Check if they exist from 'input_phase_props.m' or define them
if ~exist('I','var') || ~exist('J','var') || ~exist('K','var')
    I = eye(6);
    J = zeros(6,6); J(1:3,1:3) = 1/3;
    K = zeros(6,6); K(1:3,1:3) = diag([2/3, 2/3, 2/3]) - 1/3; K(4:6,4:6) = eye(3);
    disp('Defined standard identity tensors I, J, K.');
end


%% Main Calculation Loops
for wcit = 1:length(wc_list)
    wc = wc_list(wcit);
    % xi_list is now fixed, defined in the input section
    for xiit = 1:length(xi_list)
        xi = xi_list(xiit);
        % Check if xi exceeds physical limit for the current wc (optional but good)
        xi_ult_physical = wc / 0.42; % Approx. physical limit
        if xi > xi_ult_physical + 1e-6 % Add tolerance for floating point comparison
             disp(['Skipping calculation for wc=',num2str(wc),' xi=',num2str(xi),...
                   ' (xi exceeds physical limit ', num2str(xi_ult_physical),')']);
             continue; % Skip to the next iteration
        end

        for Fit = 1:length(Fpor_list)
            Fpor = Fpor_list(Fit);

            % --- Check if results already exist ---
            % Use try-catch for robustness if structure fields might not exist yet
            calc_exist = false; % Assume not calculated
            try
                if ~isnan(outputITZ_cell{wcit,xiit,Fit}.calc_cp.diffQvol)
                    calc_exist = true;
                end
            catch
                % Field doesn't exist or error accessing - means not calculated
                calc_exist = false;
            end

            if ~calc_exist % --- Calculate if results don't exist ---
                disp(['Calculating for wc=',num2str(wc),' xi=',num2str(xi),' Fpor=',num2str(Fpor)]);

                % --- Powers model -> cement paste volumes ---
                % Using the density variables defined earlier
                wceff = wc; % Effective w/c ratio

                % --- Calculate initial volume fractions based on Powers model ---
                denominator = (1 + rhocem / rhoH2O * wceff);
                fcem_PA = ((1 - xi) / denominator); % Cement fraction in initial paste
                fhyd_PA = (1.42 * rhocem * xi / (rhohyd * denominator)); % Hydrate fraction in initial paste
                fH2O_PA = (rhocem * (wceff - 0.42 * xi) / (rhoH2O * denominator)); % Water fraction in initial paste
                fH2O_PA = max(0, fH2O_PA); % Ensure water fraction isn't negative if xi is high
                fair_PA = 1 - fcem_PA - fhyd_PA - fH2O_PA; % Air fraction (remaining)
                fair_PA = max(0, fair_PA); % Ensure air fraction isn't negative

                % --- Calculate Cement Paste (cp) and Hydrate Foam (hf) fractions ---
                % This uses the exact formulas from the original script.
                % Their physical derivation might be specific to the model's assumptions.
                f_porosity_initial = fH2O_PA + fair_PA;
                f_solids_initial = 1 - f_porosity_initial; % = fcem_PA + fhyd_PA

                if abs(f_solids_initial) < 1e-10 % Avoid division by zero if no solids
                    warning('Initial solid fraction is near zero. Setting cp/hf fractions to NaN.');
                    fcem_cp = NaN;
                    fhf_cp = NaN;
                else
                    % Original formula for fcem_cp
                    fcem_cp = (1 - f_porosity_initial * Fpor) / f_solids_initial * fcem_PA;
                    % Original formula for fhf_cp (ensures fhf_cp = 1 - fcem_cp IF interpreted relative to a specific volume)
                    % However, let's ensure they sum correctly based on the implied volumes.
                    % If fcem_cp = Vcem / V_new and fhf_cp = Vhf / V_new, then fhf_cp should be calculated independently or ensure normalization.
                    % Let's stick to the original script's sequence:
                    fhf_cp = 1 - fcem_cp; % This assumes fcem_cp + fhf_cp = 1, defining the relative fractions.
                end

                % Original formulas for hydrate foam fractions
                if abs(fhf_cp) < 1e-10 % Avoid division by zero
                     warning('fhf_cp is near zero. Setting hf fractions to NaN.');
                     fpor_hf = NaN;
                     fhyd_hf = NaN;
                else
                     % Original formula for fpor_hf
                     fpor_hf = (fH2O_PA + fair_PA) * Fpor / fhf_cp; % Scaled porosity relative to fhf_cp volume?
                     % Original formula for fhyd_hf
                     fhyd_hf = 1 - fpor_hf; % Assumes fpor_hf + fhyd_hf = 1 within the hf phase.
                end

                % --- Checks for physical validity ---
                % Check that there is no negative fraction!!
                if ~isnan(fcem_cp) && fcem_cp < -1e-9 % Allow for small numerical inaccuracies
                    warning('Porosity factor Fpor=%.2f results in negative fcem_cp=%.2e. Check inputs/model. Setting fractions to NaN.', Fpor, fcem_cp);
                     fcem_cp = NaN; fhf_cp=NaN; fpor_hf=NaN; fhyd_hf=NaN; % Invalidate all subsequent steps
                     % Consider using 'continue' to skip to next loop iteration instead
                     % continue;
                end
                 % Optional: Add checks for other fractions being negative or > 1+eps if needed
                 if ~isnan(fpor_hf) && fpor_hf < -1e-9
                     warning('Negative fpor_hf=%.2e calculated. Setting hf fractions to NaN.', fpor_hf);
                     fpor_hf=NaN; fhyd_hf=NaN;
                 end
                 if ~isnan(fhyd_hf) && fhyd_hf < -1e-9
                     warning('Negative fhyd_hf=%.2e calculated. Setting hf fractions to NaN.', fhyd_hf);
                     fpor_hf=NaN; fhyd_hf=NaN;
                 end


                %% ELASTICITY HYDRATE FOAM (Self-Consistent)
                if fhyd_hf < 1e-9 && fpor_hf < 1e-9 % Handle case if HF is effectively empty
                    Chom_hf = zeros(6,6); % Or some other default
                    khom_hf = 0; muhom_hf = 0; vhom_hf = 0; Ehom_hf = 0;
                    disp('Hydrate foam phase empty, setting properties to zero.');
                elseif fpor_hf > 1-1e-9 % Handle case of pure pores
                    Chom_hf = Cp;
                    khom_hf = k_H2O; muhom_hf = mu_H2O;
                    vhom_hf = (3*khom_hf-2*muhom_hf)/(6*khom_hf+2*muhom_hf);
                    Ehom_hf = 9*khom_hf*muhom_hf/(3*khom_hf+muhom_hf);
                    if isnan(vhom_hf) || isinf(vhom_hf), vhom_hf=0.5; end % Handle fluid case for Poisson's ratio
                    if isnan(Ehom_hf) || isinf(Ehom_hf), Ehom_hf=0; end % Handle fluid case for Young's modulus
                    disp('Hydrate foam phase is pure (scaled) pores.');
                elseif fhyd_hf > 1-1e-9 % Handle case of pure hydrates
                    Chom_hf = Chyd;
                    khom_hf = khyd; muhom_hf = muhyd;
                    vhom_hf = nuhyd; Ehom_hf = Ehyd;
                    disp('Hydrate foam phase is pure hydrates.');
                else % Proceed with homogenization
                    C0 = fhyd_hf * Chyd; % Initial guess (Voigt bound)
                    abweichung = 1;
                    iter_sc = 0;
                    max_iter_sc = 100; % Add iteration limit for SC

                    while abweichung > tol_iso && iter_sc < max_iter_sc
                        iter_sc = iter_sc + 1;
                        % Spherical pores (using Cp = water/fluid stiffness now)
                        P_p = fun_P_sphere_iso(C0); % Eshelby tensor for sphere in matrix C0
                        % Need to handle inversion carefully if (Cp-C0) is singular
                        % Using try-catch or checking condition number might be needed
                        try
                           Ainf_p = inv(I + P_p * (Cp - C0));
                        catch ME
                           warning('Matrix inversion failed for Ainf_p in SC iso: %s. Using pseudo-inverse.', ME.identifier);
                           Ainf_p = pinv(I + P_p * (Cp - C0)); % Use pseudo-inverse as fallback
                        end


                        % Acicular hydrates (isotropic orientation average)
                        Ainf_hyd = fun_Ainf_needle_iso(Chyd, C0); % Average for needles

                        % Strain concentration tensors (Self-Consistent)
                        EEinfty_hf = inv(fpor_hf * Ainf_p + fhyd_hf * Ainf_hyd);
                        A_p = Ainf_p * EEinfty_hf;
                        A_hyd = Ainf_hyd * EEinfty_hf;

                        % Homogenized stiffness - SELF CONSISTENT
                        % Note: Pores now contribute if Cp is not zero
                        Chom_hf = fpor_hf * Cp * A_p + fhyd_hf * Chyd * A_hyd;

                        % Update and check convergence
                        C0_old = C0;
                        C0 = Chom_hf;
                        % Check norms are non-zero before dividing
                        norm_C0 = norm(C0);
                        if norm_C0 < 1e-12 % Avoid division by zero if stiffness is near zero
                            abweichung = norm(C0 - C0_old);
                        else
                            abweichung = norm(C0 - C0_old) / norm_C0;
                        end
                        % Ensure Chom_hf remains symmetric and positive definite (optional checks)

                    end % End while loop SC iso
                    if iter_sc == max_iter_sc
                       warning('Max iterations reached for SC isotropic homogenization.');
                    end


                    % Extract effective isotropic properties from Chom_hf
                    muhom_hf = C0(6,6)/2; % Shear modulus G = C44 for isotropic
                    khom_hf = C0(1,1) - 4/3 * muhom_hf; % Bulk modulus K = C11 - 4/3 G
                    % Check for physically unrealistic negative moduli
                    if muhom_hf < 0 || khom_hf < 0
                       warning('Negative moduli calculated for HF: K=%.2f, G=%.2f. Setting to small positive.', khom_hf, muhom_hf);
                       muhom_hf = max(muhom_hf, 1e-6);
                       khom_hf = max(khom_hf, 1e-6);
                    end
                    % Calculate E and nu
                    denominator_E = (3 * khom_hf + muhom_hf);
                    if abs(denominator_E) < 1e-10
                         Ehom_hf = 0; % Avoid division by zero
                    else
                         Ehom_hf = 9 * khom_hf * muhom_hf / denominator_E;
                    end
                    denominator_nu = (2 * (3 * khom_hf + muhom_hf));
                     if abs(denominator_nu) < 1e-10
                         vhom_hf = 0.5; % Incompressible limit often represented as 0.5
                     else
                         vhom_hf = (3 * khom_hf - 2 * muhom_hf) / denominator_nu;
                     end
                    % Clamp Poisson's ratio to valid range [-1, 0.5] if needed
                    vhom_hf = max(-0.99, min(0.5, vhom_hf));

                end % End if/else for HF homogenization types

                %% ELASTICITY CEMENT PASTE (Mori-Tanaka)
                if fcem_cp < 1e-9 % Handle case of pure hydrate foam
                    Chom_cp = Chom_hf;
                    khom_cp = khom_hf; muhom_cp = muhom_hf;
                    vhom_cp = vhom_hf; Ehom_cp = Ehom_hf;
                    disp('Cement paste is pure hydrate foam.');
                elseif fhf_cp < 1e-9 % Handle case of pure cement
                     Chom_cp = Cclin;
                     khom_cp = kclin; muhom_cp = muclin;
                     vhom_cp = nuclin; Ehom_cp = Eclin;
                     disp('Cement paste is pure clinker.');
                else % Proceed with Mori-Tanaka
                    C0 = Chom_hf; % Matrix is the homogenized hydrate foam

                    % Spherical clinkers (inclusions)
                    P_sph = fun_P_sphere_iso(C0); % Eshelby tensor for sphere in matrix C0
                    Ainf_cem = inv(I + P_sph * (Cclin - C0)); % Dilute strain conc. for clinker
                    Ainf_hf = I; % Strain conc. for matrix phase is identity

                    % Strain concentration tensors (Mori-Tanaka)
                    EEinfty_cp = inv(fcem_cp * Ainf_cem + fhf_cp * Ainf_hf);
                    A_cem = Ainf_cem * EEinfty_cp; % Strain conc. for clinker
                    A_hf = Ainf_hf * EEinfty_cp; % Strain conc. for hydrate foam matrix

                    % Homogenized stiffness (MORI-TANAKA)
                    Chom_cp = fcem_cp * Cclin * A_cem + fhf_cp * Chom_hf * A_hf;

                     % Extract effective isotropic properties from Chom_cp
                     muhom_cp = Chom_cp(6,6)/2;
                     khom_cp = Chom_cp(1,1) - 4/3 * muhom_cp;
                     % Check for physically unrealistic negative moduli
                     if muhom_cp < 0 || khom_cp < 0
                        warning('Negative moduli calculated for CP: K=%.2f, G=%.2f. Setting to small positive.', khom_cp, muhom_cp);
                        muhom_cp = max(muhom_cp, 1e-6);
                        khom_cp = max(khom_cp, 1e-6);
                     end

                     % Calculate E and nu
                     denominator_E = (3 * khom_cp + muhom_cp);
                      if abs(denominator_E) < 1e-10
                          Ehom_cp = 0;
                      else
                          Ehom_cp = 9 * khom_cp * muhom_cp / denominator_E;
                      end
                     denominator_nu = (2 * (3 * khom_cp + muhom_cp));
                     if abs(denominator_nu) < 1e-10
                         vhom_cp = 0.5;
                     else
                         vhom_cp = (3 * khom_cp - 2 * muhom_cp) / denominator_nu;
                     end
                     vhom_cp = max(-0.99, min(0.5, vhom_cp)); % Clamp Poisson's ratio
                end % End if/else for CP homogenization


                %% SECOND ORDER STRENGTH CALCULATIONS (Difference Quotients)
                % These calculations involve anisotropy and are more complex.
                % They require functions like P_isotrans_sph, fun_Q4_bp, fun_P_ellipsoid_transiso
                % and the Stroud points (stroud_azi, stroud_zeni)

                % --- DEVIATORIC Perturbation ---
                disp('Starting Deviatoric Perturbation Calculation...');
                Chom_hf_dev_transiso = nan(6,6); % Initialize as NaN
                Chom_cp_dev_transiso = nan(6,6); % Initialize as NaN

                try % Wrap potentially complex calculation in try-catch
                    % 1. Hydrate Foam (Anisotropic SC with perturbed muhyd)
                    C0 = Chom_hf; % Start from isotropic HF result
                    Chyd_diff = 3 * khyd * J + 2 * muhyddiff * K; % Perturbed hydrate stiffness (shear)
                    abweichung = 1;
                    countermark = 0;
                    while abweichung > tol_aniso && countermark < maxcounter
                        %disp(['DEV: Fehler=',num2str(abweichung,6),' > tol_aniso=',num2str(tol_aniso),'  C1111=',num2str(C0(1,1),10),'  C3333=',num2str(C0(3,3),10)])
                        countermark = countermark + 1;

                        % Pores (Spherical, use P_isotrans_sph for consistency)
                        P_p = P_isotrans_sph(C0);
                        Ainf_p = inv(I + P_p * (Cp - C0));

                        % Hydrates (Needle, average over Stroud points)
                        sumAinf_hyd = zeros(6,6);
                        for i = 1:length(stroud_azi) % Use length of Stroud points array
                            azi = stroud_azi(i); zeni = stroud_zeni(i);
                            Q4 = fun_Q4_bp(azi, zeni); Q4t = transpose(Q4); % Rotation matrix
                            C0_aniso = Q4 * C0 * Q4t; % Rotate matrix stiffness
                            P_hyd_i_e3 = fun_P_ellipsoid_transiso(C0_aniso, 1, 1e-20); % Eshelby for needle in rotated frame
                            P_hyd_i = (Q4t * P_hyd_i_e3 * Q4); % Rotate Eshelby tensor back
                            Ainf_hyd_i = inv(I + P_hyd_i * (Chyd - C0)); % Strain conc. for this orientation
                            sumAinf_hyd = sumAinf_hyd + Ainf_hyd_i; % Accumulate
                        end
                        num_stroud = length(stroud_azi); % Number of Stroud points used
                        sumAinf_hyd = (1 / num_stroud) * sumAinf_hyd; % Average

                        % Perturbed hydrate fraction (single orientation, e3)
                        P_hyd_e3 = fun_P_ellipsoid_transiso(C0, 1, 1e-20); % Eshelby for needle along e3
                        Ainf_hyd_diff = inv(I + P_hyd_e3 * (Chyd_diff - C0)); % For perturbed stiffness
                        Ainf_hyd_minus = inv(I + P_hyd_e3 * (Chyd - C0)); % For original stiffness (to subtract)

                        % Strain concentration tensors (Anisotropic SC)
                        % Note: Using fhyddiff volume fraction for the perturbed part
                        M_inv = (fpor_hf * Ainf_p + ...
                                (fhyd_hf - fhyddiff) * sumAinf_hyd + ... % Regular hydrates (adjusted fraction)
                                fhyddiff * Ainf_hyd_diff - ... % Contribution from added perturbed part
                                0); % Original code subtracted Ainf_hyd_minus here, seems unusual for SC? Check derivation.
                                    % Let's follow original structure for now, but be aware:
                        M_inv_orig = (fpor_hf * Ainf_p + ...
                                fhyd_hf * sumAinf_hyd + ...
                                fhyddiff * Ainf_hyd_diff - ...
                                fhyddiff * Ainf_hyd_minus);

                        EEinfty_hf = inv(M_inv_orig); % Use original formulation's structure matrix
                        A_p = Ainf_p * EEinfty_hf;
                        A_hyd = sumAinf_hyd * EEinfty_hf; % Average hydrate concentration
                        A_hyd_diff = Ainf_hyd_diff * EEinfty_hf;
                        A_hyd_minus = Ainf_hyd_minus * EEinfty_hf;

                        % Homogenized stiffness (Anisotropic SC)
                        Chom_hf_dev_transiso = (fpor_hf * Cp * A_p + ... % Pores contribution
                                               (fhyd_hf-fhyddiff) * Chyd * A_hyd + ... % Unperturbed hydrates
                                               fhyddiff * Chyd_diff * A_hyd_diff - ... % Perturbed hydrates
                                               0 ); % Again, original had a subtraction term here
                        % Using original formulation:
                         Chom_hf_dev_transiso_orig = (fpor_hf*Cp*A_p + ... % Added Cp term if pores have stiffness
                                                    fhyd_hf*Chyd*A_hyd + ...
                                                    fhyddiff*Chyd_diff*A_hyd_diff - ...
                                                    fhyddiff*Chyd*A_hyd_minus);


                        % Update and check convergence
                        C0_old = C0;
                        C0 = Chom_hf_dev_transiso_orig; % Update using the calculated value
                         norm_C0 = norm(C0);
                        if norm_C0 < 1e-12
                            abweichung = norm(C0 - C0_old);
                        else
                            abweichung = norm(C0 - C0_old) / norm_C0;
                        end
                    end % End while Aniso SC dev
                    if countermark == maxcounter
                       warning('Max iterations reached for SC anisotropic homogenization (Deviatoric).');
                    end
                    Chom_hf_dev_transiso = C0; % Store the converged result

                    % 2. Cement Paste (Anisotropic MT with perturbed HF matrix)
                    C0 = Chom_hf_dev_transiso; % Matrix is the anisotropic HF result
                    P_sph = P_isotrans_sph(C0); % Eshelby tensor (use anisotropic version)
                    Ainf_cem = inv(I + P_sph * (Cclin - C0)); % Clinker inclusion
                    Ainf_hf = I; % Matrix
                    EEinfty_cp_dev = inv(fcem_cp * Ainf_cem + fhf_cp * Ainf_hf);
                    A_cem_dev = Ainf_cem * EEinfty_cp_dev;
                    A_hf_dev = Ainf_hf * EEinfty_cp_dev;
                    Chom_cp_dev_transiso = fcem_cp * Cclin * A_cem_dev + fhf_cp * C0 * A_hf_dev; % MT result

                catch ME_dev
                   warning('Error during Deviatoric perturbation calculation: %s. Skipping.', ME_dev.message);
                   % Ensure results remain NaN if calculation failed
                   Chom_hf_dev_transiso = nan(6,6);
                   Chom_cp_dev_transiso = nan(6,6);
                end % End try-catch dev


                % --- VOLUMETRIC Perturbation ---
                disp('Starting Volumetric Perturbation Calculation...');
                Chom_hf_vol_transiso = nan(6,6); % Initialize as NaN
                Chom_cp_vol_transiso = nan(6,6); % Initialize as NaN

                 try % Wrap potentially complex calculation in try-catch
                    % 1. Hydrate Foam (Anisotropic SC with perturbed khyd)
                    C0 = Chom_hf; % Start from isotropic HF result
                    Chyd_diff = 3 * khyddiff * J + 2 * muhyd * K; % Perturbed hydrate stiffness (bulk)
                    abweichung = 1;
                    countermark = 0;
                    while abweichung > tol_aniso && countermark < maxcounter
                        %disp(['VOL: Fehler=',num2str(abweichung,6),' > tol_iso=',num2str(tol_aniso),'  C1111=',num2str(C0(1,1),10),'  C3333=',num2str(C0(3,3),10)])
                        countermark = countermark + 1;
                         % --- Repeat anisotropic SC calculation structure as in Deviatoric ---
                         P_p = P_isotrans_sph(C0);
                         Ainf_p = inv(I + P_p * (Cp - C0));
                         sumAinf_hyd = zeros(6,6);
                         for i = 1:length(stroud_azi)
                            azi = stroud_azi(i); zeni = stroud_zeni(i);
                            Q4 = fun_Q4_bp(azi, zeni); Q4t = transpose(Q4);
                            C0_aniso = Q4 * C0 * Q4t;
                            P_hyd_i_e3 = fun_P_ellipsoid_transiso(C0_aniso, 1, 1e-20);
                            P_hyd_i = (Q4t * P_hyd_i_e3 * Q4);
                            Ainf_hyd_i = inv(I + P_hyd_i * (Chyd - C0));
                            sumAinf_hyd = sumAinf_hyd + Ainf_hyd_i;
                         end
                         num_stroud = length(stroud_azi);
                         sumAinf_hyd = (1 / num_stroud) * sumAinf_hyd;
                         P_hyd_e3 = fun_P_ellipsoid_transiso(C0, 1, 1e-20);
                         Ainf_hyd_diff = inv(I + P_hyd_e3 * (Chyd_diff - C0)); % Use VOLUMETRIC Chyd_diff here
                         Ainf_hyd_minus = inv(I + P_hyd_e3 * (Chyd - C0));

                         M_inv_orig = (fpor_hf * Ainf_p + ...
                                fhyd_hf * sumAinf_hyd + ...
                                fhyddiff * Ainf_hyd_diff - ... % Ainf based on VOL Chyd_diff
                                fhyddiff * Ainf_hyd_minus);
                         EEinfty_hf = inv(M_inv_orig);
                         A_p = Ainf_p * EEinfty_hf;
                         A_hyd = sumAinf_hyd * EEinfty_hf;
                         A_hyd_diff = Ainf_hyd_diff * EEinfty_hf; % Uses Ainf based on VOL Chyd_diff
                         A_hyd_minus = Ainf_hyd_minus * EEinfty_hf;

                         Chom_hf_vol_transiso_orig = (fpor_hf*Cp*A_p + ...
                                                    fhyd_hf*Chyd*A_hyd + ...
                                                    fhyddiff*Chyd_diff*A_hyd_diff - ... % Use VOL Chyd_diff here
                                                    fhyddiff*Chyd*A_hyd_minus);

                         C0_old = C0;
                         C0 = Chom_hf_vol_transiso_orig; % Update
                         norm_C0 = norm(C0);
                         if norm_C0 < 1e-12
                             abweichung = norm(C0 - C0_old);
                         else
                             abweichung = norm(C0 - C0_old) / norm_C0;
                         end
                    end % End while Aniso SC vol
                    if countermark == maxcounter
                       warning('Max iterations reached for SC anisotropic homogenization (Volumetric).');
                    end
                    Chom_hf_vol_transiso = C0; % Store converged result

                    % 2. Cement Paste (Anisotropic MT with perturbed HF matrix)
                    C0 = Chom_hf_vol_transiso; % Matrix is the anisotropic HF result
                    P_sph = P_isotrans_sph(C0);
                    Ainf_cem = inv(I + P_sph * (Cclin - C0));
                    Ainf_hf = I;
                    EEinfty_cp_vol = inv(fcem_cp * Ainf_cem + fhf_cp * Ainf_hf);
                    A_cem_vol = Ainf_cem * EEinfty_cp_vol;
                    A_hf_vol = Ainf_hf * EEinfty_cp_vol;
                    Chom_cp_vol_transiso = fcem_cp * Cclin * A_cem_vol + fhf_cp * C0 * A_hf_vol; % MT result

                 catch ME_vol
                    warning('Error during Volumetric perturbation calculation: %s. Skipping.', ME_vol.message);
                    % Ensure results remain NaN if calculation failed
                    Chom_hf_vol_transiso = nan(6,6);
                    Chom_cp_vol_transiso = nan(6,6);
                 end % End try-catch vol


                %% Difference Quotients Calculation
                % Check if prerequisite calculations were successful
                if ~any(isnan(Chom_cp(:))) && ~any(isnan(Chom_cp_vol_transiso(:)))
                    diffQ_vol_cp = (1 / (fapp * khyd)) * (Chom_cp_vol_transiso - Chom_cp);
                else
                    diffQ_vol_cp = nan(6,6); % Assign NaN if inputs are invalid
                    warning('Cannot calculate Volumetric Diff Quotient due to NaNs in stiffness matrices.');
                end

                if ~any(isnan(Chom_cp(:))) && ~any(isnan(Chom_cp_dev_transiso(:)))
                    diffQ_dev_cp = (1 / (fapp * muhyd)) * (Chom_cp_dev_transiso - Chom_cp);
                else
                    diffQ_dev_cp = nan(6,6); % Assign NaN if inputs are invalid
                     warning('Cannot calculate Deviatoric Diff Quotient due to NaNs in stiffness matrices.');
                end

                %% Save results in cell structure (Recommended Structure - using 3D indexing)
                outputITZ_cell{wcit, xiit, Fit}.wc = wc; % Use 3 indices
                outputITZ_cell{wcit, xiit, Fit}.xi = xi;
                outputITZ_cell{wcit, xiit, Fit}.Fpor = Fpor; % Save Fpor
                % Hydrate Foam results (Simplified)
                outputITZ_cell{wcit, xiit, Fit}.calc_hf.E = Ehom_hf;
                outputITZ_cell{wcit, xiit, Fit}.calc_hf.nu = vhom_hf;
                % Hydrate Foam volume fractions
                outputITZ_cell{wcit, xiit, Fit}.vol_hf.hyd = fhyd_hf; % Saved under vol_hf.hyd
                outputITZ_cell{wcit, xiit, Fit}.vol_hf.por = fpor_hf; % Saved under vol_hf.por
                % Cement Paste results
                outputITZ_cell{wcit, xiit, Fit}.calc_cp.E = Ehom_cp;
                outputITZ_cell{wcit, xiit, Fit}.calc_cp.nu = vhom_cp;
                outputITZ_cell{wcit, xiit, Fit}.calc_cp.diffQvol = diffQ_vol_cp;
                outputITZ_cell{wcit, xiit, Fit}.calc_cp.diffQdev = diffQ_dev_cp;
                % Cement Paste volume fractions
                outputITZ_cell{wcit, xiit, Fit}.vol_cp.cem = fcem_cp; % Saved under vol_cp.cem
                outputITZ_cell{wcit, xiit, Fit}.vol_cp.hf = fhf_cp;  % Saved under vol_cp.hf
                % Precision and numerics
                outputITZ_cell{wcit, xiit, Fit}.precision = [tol_iso, tol_aniso, maxcounter];
                outputITZ_cell{wcit, xiit, Fit}.numerics = [fapp, fhyddiff];
                
                % --- Save intermediate results to file ---
                save(filename, 'outputITZ_cell'); % Save the 3D cell array
                disp(['Results saved to ', filename])

            else % --- calc_exist was true ---
                disp(['Skipping calculation for wc=',num2str(wc),' xi=',num2str(xi),' Fpor=',num2str(Fpor),' (already exists).']);
            end % End if calc_exist==0

        end % End Fpor loop
    end % End xi loop
end % End wc loop

disp('--- Calculation Finished ---');
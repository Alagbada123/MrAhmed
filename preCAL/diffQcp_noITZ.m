%% PRECALCULATION of the difference Quotients for second order stress concentration to hydrates
% --- VERSION WITHOUT Fpor ---
clear all; clc; close all

% --- Add path for helper functions ---
addpath('../functions'); % Ensure necessary .m files are here
% --- Assuming stroud_points script defines stroud_azi, stroud_zeni ---
stroud_points; % Keep this if needed for anisotropic calculations

%% Input Parameters (Updated Section - No Fpor)
% --- Grid to calculate on ---
wc_list = [0.50, 0.55, 0.60];      % NEW wc range from template
xi_list = 0.05:0.05:1.0;     % NEW fixed xi range from template

% --- Densities (Keep previously defined values) ---
rhoH2O = 1000;      % [kg/m^3] Water density
rhohyd = 2073;      % [kg/m^3] Hydrate density
rhocem = 3150;      % [kg/m^3] Cement density

% --- Phase Stiffness Properties (Keep previously defined values) ---
% Cement Clinker
Eclin = 139.9; nuclin = 0.3;
[kclin, muclin] = fun_kmu_from_Enu(Eclin, nuclin);
Cclin = fun_CfromEnu(Eclin, nuclin);
% Hydrates
Ehyd = 29.15786664; nuhyd = 0.24;
[khyd, muhyd] = fun_kmu_from_Enu(Ehyd, nuhyd);
Chyd = fun_CfromEnu(Ehyd, nuhyd);
% Pores (filled with water or not, when filled with water, then k_H2O = 2.3)
k_H2O = 2.3; mu_H2O = 0;
Cp = fun_Cfromkmu(k_H2O, mu_H2O);

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
filename = 'precalc_cpwcallSSD2.mat'; % NEW filename from template

%% initialization
% Define structure template for results
calc_struct=struct('E',NaN,'nu',NaN,'diffQvol',NaN,'diffQdev',NaN);
% Define main precalc structure template
precalc_struct=struct('wc',NaN,'xi',NaN, ...             % Inputs
                      'calc_cp',calc_struct, ...         % Results for cement paste
                      'calc_hf',struct('E',NaN,'nu',NaN), ... % Simplified results for hydrate foam (adjust if diffQ needed for hf)
                      'precision',NaN(1,3), ...         % Tolerances, maxiter
                      'numerics',NaN(1,2), ...          % fapp, fhyddiff
                      'vol_cp',struct('cem',NaN,'hf',NaN),... % Vol frac WITHIN cement paste
                      'vol_hf',struct('hyd',NaN,'por',NaN)); % Vol frac WITHIN hydrate foam

% --- Load existing data or initialize output structure (Now 2D) ---
if exist(filename, 'file') == 2
    load(filename)
    disp(['Warning: ',filename,' exists and is loaded. Appending new results if any.'])
    % Basic check if size matches new lists
    current_size = size(output_cell); % Use new variable name
    if current_size(1)~=length(wc_list) || current_size(2)~=length(xi_list)
       disp('Warning: Loaded data size does not match current parameter lists.');
       % Decide how to handle mismatch (re-initialize for simplicity)
       disp('Re-initializing output structure due to size mismatch.');
       output_cell=repmat({precalc_struct}, length(wc_list), length(xi_list)); % Now 2D
    end
else
    disp([filename,' does not exist. Initializing new output structure.'])
    output_cell=repmat({precalc_struct}, length(wc_list), length(xi_list)); % Now 2D
end

% --- Define standard identity tensors (if not defined elsewhere) ---
if ~exist('I','var') || ~exist('J','var') || ~exist('K','var')
    I = eye(6);
    J = zeros(6,6); J(1:3,1:3) = 1/3;
    K = zeros(6,6); K(1:3,1:3) = diag([2/3, 2/3, 2/3]) - 1/3; K(4:6,4:6) = eye(3);
    disp('Defined standard identity tensors I, J, K.');
end


%% Main Calculation Loops (Now 2 Loops: wc, xi)
for wcit = 1:length(wc_list)
    wc = wc_list(wcit);
    % xi_list is fixed, defined in the input section

    for xiit = 1:length(xi_list)
        xi = xi_list(xiit);

        % Optional: Check if xi exceeds physical limit for the current wc
        xi_ult_physical = wc / 0.42; % Approx. physical limit
        if xi > xi_ult_physical + 1e-6 % Add tolerance
             disp(['Skipping calculation for wc=',num2str(wc),' xi=',num2str(xi),...
                   ' (xi exceeds physical limit ', num2str(xi_ult_physical),')']);
             continue; % Skip to the next iteration
        end

        % --- Check if results already exist (using 2D indexing) ---
        calc_exist = false; % Assume not calculated
        try
            % Access using 2 indices now
            if ~isnan(output_cell{wcit, xiit}.calc_cp.diffQvol)
                calc_exist = true;
            end
        catch
            calc_exist = false;
        end

        if ~calc_exist % --- Calculate if results don't exist ---
            disp(['Calculating for wc=',num2str(wc),' xi=',num2str(xi)]);

            % --- Powers model -> cement paste volumes ---
            wceff = wc;
            denominator = (1 + rhocem / rhoH2O * wceff);
            fcem_PA = ((1 - xi) / denominator);
            fhyd_PA = (1.42 * rhocem * xi / (rhohyd * denominator));
            fH2O_PA = (rhocem * (wceff - 0.42 * xi) / (rhoH2O * denominator));
            fH2O_PA = max(0, fH2O_PA);
            fair_PA = 1 - fcem_PA - fhyd_PA - fH2O_PA;
            fair_PA = max(0, fair_PA);

            % --- Calculate CP and HF volume fractions (using template logic) ---
            fcem_cp = fcem_PA;  % Cement fraction in CP = Absolute cement fraction
            fhyd_cp = fhyd_PA;  % Hydrate fraction term used below

            % Hydrate foam includes hydrates and ALL porosity
            fhf_cp = fhyd_cp + fH2O_PA + fair_PA; % Hydrate foam fraction in CP
            fhf_cp = max(0, fhf_cp); % Ensure non-negative

            % --- Fractions within Hydrate Foam (HF = Hydrates + Pores) ---
            if abs(fhf_cp) < 1e-10 % Avoid division by zero if HF fraction is zero
                fpor_hf = NaN; % Cannot determine fractions within HF
                fhyd_hf = NaN;
                % If fhf_cp is zero, it implies fhyd_cp, fH2O_PA, fair_PA are all zero.
                % This means fcem_cp = 1 (pure clinker case).
                if abs(fcem_cp-1.0) > 1e-9
                   warning('fhf_cp is zero, but fcem_cp is not 1. Check volume fraction logic.');
                end
            else
                fpor_hf = (fH2O_PA + fair_PA) / fhf_cp; % Porosity fraction within HF phase
                fhyd_hf = fhyd_cp / fhf_cp;           % Hydrate fraction within HF phase
            end

            % --- Normalize and check HF fractions ---
             if ~isnan(fpor_hf) && ~isnan(fhyd_hf)
                 fpor_hf = max(0, fpor_hf); fhyd_hf = max(0, fhyd_hf);
                 hf_sum = fpor_hf + fhyd_hf;
                 if abs(hf_sum - 1.0) > 1e-9 && hf_sum > 1e-9 % Check if sum is close to 1
                     warning('HF fractions do not sum to 1 (Sum=%.4f). Normalizing.', hf_sum);
                     fpor_hf = fpor_hf / hf_sum;
                     fhyd_hf = fhyd_hf / hf_sum;
                 elseif hf_sum < 1e-9 % If sum is zero, fractions should be zero
                     fpor_hf=0; fhyd_hf=0;
                 end
             else
                 % If fractions are NaN, likely due to fhf_cp being zero.
                 % Handle subsequent calculations appropriately.
                 warning('Cannot calculate HF fractions as fhf_cp is near zero.');
             end

             % --- Check CP fractions ---
             % Based on this logic: fcem_cp + fhf_cp = fcem_PA + fhyd_PA + fH2O_PA + fair_PA = 1.0
             cp_sum = fcem_cp + fhf_cp;
             if abs(cp_sum - 1.0) > 1e-9
                  warning('CP fractions do not sum to 1 (Sum=%.4f). Check logic.', cp_sum);
                  % Decide on handling: normalize? error? For now, proceed.
             end

            %% ELASTICITY HYDRATE FOAM (Self-Consistent)
            % --- Handle cases where fractions might be NaN ---
            if isnan(fhyd_hf) || isnan(fpor_hf)
                Chom_hf = nan(6,6);
                khom_hf = NaN; muhom_hf = NaN; vhom_hf = NaN; Ehom_hf = NaN;
                disp('Skipping HF homogenization due to invalid volume fractions.');
            elseif fpor_hf > 1-1e-9 % Handle case of pure pores
                Chom_hf = Cp;
                khom_hf = k_H2O; muhom_hf = mu_H2O; vhom_hf=0.5; Ehom_hf=0;
                disp('Hydrate foam phase is pure pores.');
            elseif fhyd_hf > 1-1e-9 % Handle case of pure hydrates
                Chom_hf = Chyd;
                khom_hf = khyd; muhom_hf = muhyd; vhom_hf = nuhyd; Ehom_hf = Ehyd;
                disp('Hydrate foam phase is pure hydrates.');
            else % --- Proceed with homogenization ---
                C0 = fhyd_hf * Chyd; % Initial guess
                abweichung = 1; iter_sc = 0; max_iter_sc = 100;
                while abweichung > tol_iso && iter_sc < max_iter_sc
                    iter_sc = iter_sc + 1;
                    try
                        P_p = fun_P_sphere_iso(C0);
                        Ainf_p = inv(I + P_p * (Cp - C0));
                    catch ME
                        warning('Matrix inversion failed for Ainf_p in SC iso: %s. Using pseudo-inverse.', ME.identifier);
                        Ainf_p = pinv(I + P_p * (Cp - C0));
                    end
                    Ainf_hyd = fun_Ainf_needle_iso(Chyd, C0);
                    EEinfty_hf = inv(fpor_hf * Ainf_p + fhyd_hf * Ainf_hyd);
                    A_p = Ainf_p * EEinfty_hf;
                    A_hyd = Ainf_hyd * EEinfty_hf;
                    Chom_hf = fpor_hf * Cp * A_p + fhyd_hf * Chyd * A_hyd; % SC formula
                    C0_old = C0; C0 = Chom_hf;
                    norm_C0 = norm(C0);
                    if norm_C0 < 1e-12, abweichung = norm(C0 - C0_old);
                    else, abweichung = norm(C0 - C0_old) / norm_C0; end
                end % End while SC iso
                if iter_sc == max_iter_sc, warning('Max iterations reached for SC isotropic homogenization.'); end

                muhom_hf = C0(6,6)/2; khom_hf = C0(1,1) - 4/3 * muhom_hf;
                muhom_hf = max(muhom_hf, 1e-6); khom_hf = max(khom_hf, 1e-6); % Ensure positive
                denominator_E = (3*khom_hf + muhom_hf);
                if abs(denominator_E)<1e-10, Ehom_hf=0; else, Ehom_hf = 9*khom_hf*muhom_hf/denominator_E; end
                denominator_nu = (2*(3*khom_hf + muhom_hf));
                if abs(denominator_nu)<1e-10, vhom_hf=0.5; else, vhom_hf = (3*khom_hf-2*muhom_hf)/denominator_nu; end
                vhom_hf = max(-0.99, min(0.5, vhom_hf)); % Clamp Poisson's ratio
            end % End if/else for HF homogenization types

            %% ELASTICITY CEMENT PASTE (Mori-Tanaka)
            % --- Handle cases where input Chom_hf might be NaN ---
             if isnan(fcem_cp) || isnan(fhf_cp) || any(isnan(Chom_hf(:)))
                 Chom_cp = nan(6,6);
                 khom_cp = NaN; muhom_cp = NaN; vhom_cp = NaN; Ehom_cp = NaN;
                 disp('Skipping CP homogenization due to invalid inputs.');
             elseif fcem_cp < 1e-9 % Handle case of pure hydrate foam (fhf_cp approx 1)
                 Chom_cp = Chom_hf;
                 khom_cp = khom_hf; muhom_cp = muhom_hf; vhom_cp = vhom_hf; Ehom_cp = Ehom_hf;
                 disp('Cement paste is pure hydrate foam.');
             elseif fhf_cp < 1e-9 % Handle case of pure cement (fcem_cp approx 1)
                 Chom_cp = Cclin;
                 khom_cp = kclin; muhom_cp = muclin; vhom_cp = nuclin; Ehom_cp = Eclin;
                 disp('Cement paste is pure clinker.');
             else % --- Proceed with Mori-Tanaka ---
                 C0 = Chom_hf; % Matrix is the homogenized hydrate foam
                 P_sph = fun_P_sphere_iso(C0);
                 Ainf_cem = inv(I + P_sph * (Cclin - C0)); % Clinker inclusion
                 Ainf_hf = I; % Matrix phase
                 EEinfty_cp = inv(fcem_cp * Ainf_cem + fhf_cp * Ainf_hf);
                 A_cem = Ainf_cem * EEinfty_cp;
                 A_hf = Ainf_hf * EEinfty_cp;
                 Chom_cp = fcem_cp * Cclin * A_cem + fhf_cp * Chom_hf * A_hf; % MT formula

                 muhom_cp = Chom_cp(6,6)/2; khom_cp = Chom_cp(1,1) - 4/3 * muhom_cp;
                 muhom_cp = max(muhom_cp, 1e-6); khom_cp = max(khom_cp, 1e-6); % Ensure positive
                 denominator_E = (3*khom_cp + muhom_cp);
                 if abs(denominator_E)<1e-10, Ehom_cp=0; else, Ehom_cp = 9*khom_cp*muhom_cp/denominator_E; end
                 denominator_nu = (2*(3*khom_cp + muhom_cp));
                 if abs(denominator_nu)<1e-10, vhom_cp=0.5; else, vhom_cp = (3*khom_cp-2*muhom_cp)/denominator_nu; end
                 vhom_cp = max(-0.99, min(0.5, vhom_cp)); % Clamp Poisson's ratio
             end % End if/else for CP homogenization


            %% SECOND ORDER STRENGTH CALCULATIONS (Difference Quotients)
            % --- Initialize results ---
            Chom_cp_dev_transiso = nan(6,6);
            Chom_cp_vol_transiso = nan(6,6);
            diffQ_vol_cp = nan(6,6);
            diffQ_dev_cp = nan(6,6);

            % --- Check if prerequisites are valid before attempting calculation ---
            if ~any(isnan(Chom_cp(:))) && ~any(isnan(Chom_hf(:))) && ...
               ~isnan(fpor_hf) && ~isnan(fhyd_hf) && ~isnan(fcem_cp) && ~isnan(fhf_cp)

                % --- DEVIATORIC Perturbation ---
                disp('Starting Deviatoric Perturbation Calculation...');
                try
                    % 1. Hydrate Foam (Anisotropic SC...) - Calculation as before
                    C0 = Chom_hf; Chyd_diff = 3*khyd*J + 2*muhyddiff*K;
                    abweichung=1; countermark=0;
                    while abweichung > tol_aniso && countermark < maxcounter
                       countermark = countermark+1;
                       P_p=P_isotrans_sph(C0); Ainf_p=inv(I+P_p*(Cp-C0));
                       sumAinf_hyd=zeros(6,6);
                       for i=1:length(stroud_azi)
                          azi=stroud_azi(i); zeni=stroud_zeni(i); Q4=fun_Q4_bp(azi,zeni); Q4t=transpose(Q4);
                          C0_aniso=Q4*C0*Q4t; P_hyd_i_e3=fun_P_ellipsoid_transiso(C0_aniso,1,1e-20);
                          P_hyd_i=(Q4t*P_hyd_i_e3*Q4); Ainf_hyd_i=inv(I+P_hyd_i*(Chyd-C0));
                          sumAinf_hyd=sumAinf_hyd+Ainf_hyd_i;
                       end
                       num_stroud = length(stroud_azi); sumAinf_hyd=(1/num_stroud)*sumAinf_hyd;
                       P_hyd_e3=fun_P_ellipsoid_transiso(C0,1,1e-20);
                       Ainf_hyd_diff=inv(I+P_hyd_e3*(Chyd_diff-C0)); Ainf_hyd_minus=inv(I+P_hyd_e3*(Chyd-C0));
                       M_inv_orig = (fpor_hf*Ainf_p + fhyd_hf*sumAinf_hyd + fhyddiff*Ainf_hyd_diff - fhyddiff*Ainf_hyd_minus);
                       EEinfty_hf = inv(M_inv_orig); A_p=Ainf_p*EEinfty_hf; A_hyd=sumAinf_hyd*EEinfty_hf;
                       A_hyd_diff=Ainf_hyd_diff*EEinfty_hf; A_hyd_minus=Ainf_hyd_minus*EEinfty_hf;
                       Chom_hf_dev_transiso_orig = (fpor_hf*Cp*A_p + fhyd_hf*Chyd*A_hyd + fhyddiff*Chyd_diff*A_hyd_diff - fhyddiff*Chyd*A_hyd_minus);
                       C0_old=C0; C0=Chom_hf_dev_transiso_orig;
                       norm_C0 = norm(C0); if norm_C0 < 1e-12, abweichung=norm(C0-C0_old); else, abweichung=norm(C0-C0_old)/norm_C0; end
                    end
                    if countermark==maxcounter, warning('Max iterations reached for SC aniso (Dev)'); end
                    Chom_hf_dev_transiso = C0;

                    % 2. Cement Paste (Anisotropic MT...) - Calculation as before
                    C0 = Chom_hf_dev_transiso; P_sph=P_isotrans_sph(C0);
                    Ainf_cem=inv(I+P_sph*(Cclin-C0)); Ainf_hf=I;
                    EEinfty_cp_dev=inv(fcem_cp*Ainf_cem + fhf_cp*Ainf_hf);
                    A_cem_dev=Ainf_cem*EEinfty_cp_dev; A_hf_dev=Ainf_hf*EEinfty_cp_dev;
                    Chom_cp_dev_transiso=fcem_cp*Cclin*A_cem_dev + fhf_cp*C0*A_hf_dev; % MT result
                catch ME_dev
                   warning('Error during Deviatoric perturbation calculation: %s. Skipping.', ME_dev.message);
                   Chom_cp_dev_transiso = nan(6,6); % Ensure result is NaN on error
                end

                % --- VOLUMETRIC Perturbation ---
                disp('Starting Volumetric Perturbation Calculation...');
                try
                    % 1. Hydrate Foam (Anisotropic SC...) - Calculation as before
                     C0 = Chom_hf; Chyd_diff = 3*khyddiff*J + 2*muhyd*K; % Vol perturbation
                     abweichung=1; countermark=0;
                     while abweichung > tol_aniso && countermark < maxcounter
                       countermark = countermark+1;
                       P_p=P_isotrans_sph(C0); Ainf_p=inv(I+P_p*(Cp-C0));
                       sumAinf_hyd=zeros(6,6);
                       for i=1:length(stroud_azi)
                          azi=stroud_azi(i); zeni=stroud_zeni(i); Q4=fun_Q4_bp(azi,zeni); Q4t=transpose(Q4);
                          C0_aniso=Q4*C0*Q4t; P_hyd_i_e3=fun_P_ellipsoid_transiso(C0_aniso,1,1e-20);
                          P_hyd_i=(Q4t*P_hyd_i_e3*Q4); Ainf_hyd_i=inv(I+P_hyd_i*(Chyd-C0));
                          sumAinf_hyd=sumAinf_hyd+Ainf_hyd_i;
                       end
                       num_stroud = length(stroud_azi); sumAinf_hyd=(1/num_stroud)*sumAinf_hyd;
                       P_hyd_e3=fun_P_ellipsoid_transiso(C0,1,1e-20);
                       Ainf_hyd_diff=inv(I+P_hyd_e3*(Chyd_diff-C0)); Ainf_hyd_minus=inv(I+P_hyd_e3*(Chyd-C0));
                       M_inv_orig = (fpor_hf*Ainf_p + fhyd_hf*sumAinf_hyd + fhyddiff*Ainf_hyd_diff - fhyddiff*Ainf_hyd_minus);
                       EEinfty_hf = inv(M_inv_orig); A_p=Ainf_p*EEinfty_hf; A_hyd=sumAinf_hyd*EEinfty_hf;
                       A_hyd_diff=Ainf_hyd_diff*EEinfty_hf; A_hyd_minus=Ainf_hyd_minus*EEinfty_hf;
                       Chom_hf_vol_transiso_orig = (fpor_hf*Cp*A_p + fhyd_hf*Chyd*A_hyd + fhyddiff*Chyd_diff*A_hyd_diff - fhyddiff*Chyd*A_hyd_minus); % Use VOL Chyd_diff
                       C0_old=C0; C0=Chom_hf_vol_transiso_orig;
                       norm_C0 = norm(C0); if norm_C0 < 1e-12, abweichung=norm(C0-C0_old); else, abweichung=norm(C0-C0_old)/norm_C0; end
                    end
                    if countermark==maxcounter, warning('Max iterations reached for SC aniso (Vol)'); end
                    Chom_hf_vol_transiso = C0;

                    % 2. Cement Paste (Anisotropic MT...) - Calculation as before
                    C0 = Chom_hf_vol_transiso; P_sph=P_isotrans_sph(C0);
                    Ainf_cem=inv(I+P_sph*(Cclin-C0)); Ainf_hf=I;
                    EEinfty_cp_vol=inv(fcem_cp*Ainf_cem + fhf_cp*Ainf_hf);
                    A_cem_vol=Ainf_cem*EEinfty_cp_vol; A_hf_vol=Ainf_hf*EEinfty_cp_vol;
                    Chom_cp_vol_transiso=fcem_cp*Cclin*A_cem_vol + fhf_cp*C0*A_hf_vol; % MT result
                catch ME_vol
                    warning('Error during Volumetric perturbation calculation: %s. Skipping.', ME_vol.message);
                    Chom_cp_vol_transiso = nan(6,6); % Ensure result is NaN on error
                end

                %% Difference Quotients Calculation
                if ~any(isnan(Chom_cp_vol_transiso(:)))
                    diffQ_vol_cp = (1 / (fapp * khyd)) * (Chom_cp_vol_transiso - Chom_cp);
                else
                    warning('Cannot calculate Volumetric Diff Quotient due to NaNs in stiffness matrices.');
                    diffQ_vol_cp = nan(6,6);
                end
                if ~any(isnan(Chom_cp_dev_transiso(:)))
                    diffQ_dev_cp = (1 / (fapp * muhyd)) * (Chom_cp_dev_transiso - Chom_cp);
                else
                    warning('Cannot calculate Deviatoric Diff Quotient due to NaNs in stiffness matrices.');
                    diffQ_dev_cp = nan(6,6);
                end

            else
                 warning('Skipping Difference Quotient calculation due to invalid inputs.');
            end % End check for valid inputs to DiffQ


            %% Save results in cell structure (Recommended Structure - using 2D indexing)
            output_cell{wcit, xiit}.wc = wc;
            output_cell{wcit, xiit}.xi = xi;
            % Hydrate Foam results (Simplified)
            output_cell{wcit, xiit}.calc_hf.E = Ehom_hf;
            output_cell{wcit, xiit}.calc_hf.nu = vhom_hf;
            % Hydrate Foam volume fractions
            output_cell{wcit, xiit}.vol_hf.hyd = fhyd_hf; % Saved under vol_hf.hyd
            output_cell{wcit, xiit}.vol_hf.por = fpor_hf; % Saved under vol_hf.por (Corrected)
            % Cement Paste results
            output_cell{wcit, xiit}.calc_cp.E = Ehom_cp;
            output_cell{wcit, xiit}.calc_cp.nu = vhom_cp;
            output_cell{wcit, xiit}.calc_cp.diffQvol = diffQ_vol_cp;
            output_cell{wcit, xiit}.calc_cp.diffQdev = diffQ_dev_cp;
            % Cement Paste volume fractions
            output_cell{wcit, xiit}.vol_cp.cem = fcem_cp; % Saved under vol_cp.cem
            output_cell{wcit, xiit}.vol_cp.hf = fhf_cp;  % Saved under vol_cp.hf
            % Precision and numerics
            output_cell{wcit, xiit}.precision = [tol_iso, tol_aniso, maxcounter];
            output_cell{wcit, xiit}.numerics = [fapp, fhyddiff];
            
            % --- Save intermediate results to file ---
            save(filename, 'output_cell'); % Save the 2D cell array
            disp(['Results saved to ', filename])

        else % --- calc_exist was true ---
            disp(['Skipping calculation for wc=',num2str(wc),' xi=',num2str(xi),' (already exists).']);
        end % End if calc_exist==0

    end % End xi loop
end % End wc loop

disp('--- Calculation Finished ---');
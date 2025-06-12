%% PRECALCULATION of the difference Quotients for second order stress concentration to hydrates
% --- VERSION WITH EXTENDED POWERS' MODEL ---
clear all; clc; close all

% --- Add path for helper functions ---
% Ensure functions like fun_kmu_from_Enu, fun_CfromEnu, fun_P_sphere_iso, etc. are in this folder.
addpath('../functions');
% --- Load Stroud points for anisotropic calculations ---
stroud_points;

%% Input Parameters (Updated for Extended Model)
% --- Grid to calculate on ---
wc_list = [0.50, 0.55, 0.60];
xi_list = 0.05:0.05:1.0;

% --- NEW: Extended Powers' Model Parameters ---
ac_ratio = 4.8;               % [-] Aggregate-to-cement mass ratio. Assumed constant for this example.

% Expanded list for water absorption. Includes:
% - 0.0: Saturated (SSD) or non-porous aggregate
% - 0.005: Low-absorption aggregate
% - 0.010: Medium-absorption (your previous value)
% - 0.020: High-absorption aggregate
wa0_a_list = [0.0, 0.005, 0.010, 0.020];

% Expanded list for void-filling extent. Includes:
% - 0.0: No internal curing (Normal Powers' Model)
% - 0.3: Low efficiency internal curing
% - 0.6: Medium efficiency (your previous value)
% - 0.9: High efficiency internal curing
alpha_list = [0.0, 0.3, 0.6, 0.9];

% --- Densities ---
rhoH2O = 1000;      % [kg/m^3] Water density
rhohyd = 2073;      % [kg/m^3] Hydrate density
rhocem = 3150;      % [kg/m^3] Cement density

% --- Phase Stiffness Properties ---
% Cement Clinker
Eclin = 139.9; nuclin = 0.3; [kclin, muclin] = fun_kmu_from_Enu(Eclin, nuclin); Cclin = fun_CfromEnu(Eclin, nuclin);
% Hydrates
Ehyd = 29.15786664; nuhyd = 0.24; [khyd, muhyd] = fun_kmu_from_Enu(Ehyd, nuhyd); Chyd = fun_CfromEnu(Ehyd, nuhyd);
% Pores (filled with water)
k_H2O = 2.3; mu_H2O = 0; Cp = fun_Cfromkmu(k_H2O, mu_H2O);
% Cp = zeros(6,6)

% --- Homogenization Precision ---
tol_iso = 1e-10; tol_aniso = 1e-8; maxcounter = 5;

% --- Difference Quotient Calculation Numerics ---
fapp = 0.001;
khyddiff = (1 + fapp) * khyd;
muhyddiff = (1 + fapp) * muhyd;
fhyddiff = fapp / 5;

% --- Output ---
filename = 'precalc_extended_model_full.mat'; % Use a new filename

%% Initialization
% --- NEW: Update data structure to store new parameters ---
calc_struct=struct('E',NaN,'nu',NaN,'diffQvol',NaN,'diffQdev',NaN);
precalc_struct=struct('wc',NaN,'xi',NaN, ...
                      'wa0_a', NaN, 'alpha', NaN, 'wceff', NaN, ... % NEW: Fields for extended model
                      'calc_cp',calc_struct, ...
                      'calc_hf',struct('E',NaN,'nu',NaN), ...
                      'precision',NaN(1,3), 'numerics',NaN(1,2), ...
                      'vol_cp',struct('cem',NaN,'hf',NaN), ...
                      'vol_hf',struct('hyd',NaN,'por',NaN));

% --- Load existing data or initialize output structure (Now 4D) ---
output_dims = [length(wc_list), length(xi_list), length(wa0_a_list), length(alpha_list)];
if exist(filename, 'file') == 2
    load(filename)
    disp(['Warning: ',filename,' exists and is loaded.'])
    current_size = size(output_cell);
    if ~isequal(current_size, output_dims)
       disp('Warning: Loaded data size does not match current parameter lists. Re-initializing.');
       output_cell=repmat({precalc_struct}, output_dims);
    end
else
    disp([filename,' does not exist. Initializing new 4D output structure.'])
    output_cell=repmat({precalc_struct}, output_dims);
end

% --- Define standard identity tensors ---
if ~exist('I','var'), I=eye(6); J=zeros(6,6); J(1:3,1:3)=1/3; K=I-J; disp('Defined I, J, K.'); end


%% Main Calculation Loops (Now 4 Nested Loops)
% --- NEW: Outer loops for wa0_a and alpha ---
for wait = 1:length(wa0_a_list)
    wa0_a = wa0_a_list(wait);
    for alphait = 1:length(alpha_list)
        alpha = alpha_list(alphait);

        for wcit = 1:length(wc_list)
            wc = wc_list(wcit);
            for xiit = 1:length(xi_list)
                xi = xi_list(xiit);

                % Define the 4D index for the output cell
                idx = {wcit, xiit, wait, alphait};

                % Optional: Physical limit check
                xi_ult_physical = wc / 0.42;
                if xi > xi_ult_physical + 1e-6
                     disp(['Skipping for wc=',num2str(wc),' xi=',num2str(xi), ' (exceeds physical limit)']);
                     continue;
                end

                % Check if results already exist
                calc_exist = false;
                try
                    if ~isnan(output_cell{idx{:}}.calc_cp.E)
                        calc_exist = true;
                    end
                catch
                    calc_exist = false;
                end

                if ~calc_exist
                    disp(['Calculating for wc=',num2str(wc),' xi=',num2str(xi), ...
                          ' wa0/a=',num2str(wa0_a),' alpha=',num2str(alpha)]);

                    %% --- NEW: Extended Powers' Model Volume Fraction Calculation (Explicit version) ---
                    % This block implements the model by explicitly calculating each term from the journal.
                    
                    % 1. Calculate the initial effective w/c ratio in the paste
                    wc0_eff = wc - wa0_a * ac_ratio;
                    if wc0_eff <= 0
                        warning('Initial effective w/c ratio (wc0_eff=%.3f) is non-positive. Skipping.', wc0_eff);
                        continue; % Skip this combination of parameters
                    end
                    
                    % 2. Calculate the evolving effective w/c ratio for the current xi
                    k_factor = ((rhoH2O / rhocem) + wc0_eff)* ((3.31*alpha) / (20 + 63 * wc0_eff));
                    wceff = wc0_eff + k_factor * xi;
                    
                    % 3. Calculate the volume fractions of the solid phases based on 'wceff'
                    denominator = (1 + rhocem / rhoH2O * wceff);
                    fcem_PA = (1 - xi) / denominator;
                    fhyd_PA = (1.42 * rhocem * xi) / (rhohyd * denominator);

                    % 4. Calculate the volume fractions of the two types of pore space
                    % f_cap^cp: Volume of capillary water (water not consumed by hydration)
                    fH2O_capillary = max(0, (rhocem * (wceff - 0.42 * xi)) / (rhoH2O * denominator));

                    % f_void^cp: Volume of shrinkage-induced voids (from Journal Eq. 8, or by conservation)
                    f_voids = max(0, 1 - fcem_PA - fhyd_PA - fH2O_capillary);

                    % 5. Implement Journal Equations (44) and (45) to get the final water and air content
                    fH2O_PA = fH2O_capillary + alpha * f_voids; % Eq. (44)
                    fair_PA = (1 - alpha) * f_voids;           % Eq. (45)
                    
                    %% --- END of new volume fraction calculation ---

                    % --- Calculate CP and HF volume fractions ---
                    fcem_cp = fcem_PA;
                    fhyd_cp = fhyd_PA;
                    fhf_cp = fhyd_cp + fH2O_PA + fair_PA;
                    fhf_cp = max(0, fhf_cp);

                    if abs(fhf_cp) < 1e-10
                        fpor_hf = NaN; fhyd_hf = NaN;
                    else
                        fpor_hf = (fH2O_PA + fair_PA) / fhf_cp;
                        fhyd_hf = fhyd_cp / fhf_cp;
                    end
                    
                    % --- The rest of the homogenization code remains UNCHANGED ---
                    
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
                    
                    %% Save results in the 4D cell structure
                    output_cell{idx{:}}.wc = wc;
                    output_cell{idx{:}}.xi = xi;
                    output_cell{idx{:}}.wa0_a = wa0_a; output_cell{idx{:}}.alpha = alpha; output_cell{idx{:}}.wceff = wceff;
                    output_cell{idx{:}}.calc_hf.E = Ehom_hf; output_cell{idx{:}}.calc_hf.nu = vhom_hf;
                    output_cell{idx{:}}.vol_hf.hyd = fhyd_hf; output_cell{idx{:}}.vol_hf.por = fpor_hf;
                    output_cell{idx{:}}.calc_cp.E = Ehom_cp; output_cell{idx{:}}.calc_cp.nu = vhom_cp;
                    output_cell{idx{:}}.calc_cp.diffQvol = diffQ_vol_cp;
                    output_cell{idx{:}}.calc_cp.diffQdev = diffQ_dev_cp;
                    output_cell{idx{:}}.vol_cp.cem = fcem_cp; output_cell{idx{:}}.vol_cp.hf = fhf_cp;
                    output_cell{idx{:}}.precision = [tol_iso, tol_aniso, maxcounter];
                    output_cell{idx{:}}.numerics = [fapp, fhyddiff];

                    save(filename, 'output_cell');
                    
                else
                    disp(['Skipping calculation for wc=',num2str(wc),' xi=',num2str(xi), ...
                          ' wa0/a=',num2str(wa0_a),' alpha=',num2str(alpha),' (already exists).']);
                end
            end % End xi loop
        end % End wc loop
    end % End alpha loop
end % End wa0_a loop

disp('--- Calculation Finished ---');
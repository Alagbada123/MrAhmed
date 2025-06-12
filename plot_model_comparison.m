%% SCRIPT TO PLOT AND COMPARE NORMAL VS. EXTENDED POWERS' MODEL RESULTS
clear all; clc; close all;

% --- Configuration ---
filename = 'precalc_extended_model_with_Fpor_FULL.mat'; 
load(filename);

% --- Define the parameters used in the calculation ---
wc_list = [0.50, 0.55];
xi_list = 0.05:0.05:1.0;
Fpor_list = 1.0:0.1:1.7;
wa0_a_list = [0.0, 0.01];
alpha_list = [0.0, 0.6];
ac_ratio = 4.9;

% --- Material Densities ---
rhoH2O = 1000; rhohyd = 2073; rhocem = 3150;

%% PLOT 1: REPLICATE JOURNAL FIGURE 12
disp('Generating Plot 1: Replication of Journal Figure 12...');

% --- Define cases for plotting ---
wc_normal = 0.55; wa0a_normal = 0.0; alpha_normal = 0.0;
wc_extended = 0.55; wa0a_extended = 0.01; alpha_extended = 0.6;
Fpor_plot = 1.3;
if ~ismember(Fpor_plot, Fpor_list), error('Fpor_plot=%.2f is not in your Fpor_list.', Fpor_plot); end
xi_plot_range = [0, xi_list];

% --- Calculate volumes for both cases ---
[vol_cem_norm, vol_hyd_norm, vol_wat_norm, vol_air_norm] = calculate_volumes(xi_plot_range, wc_normal, wa0a_normal, alpha_normal, ac_ratio, Fpor_plot, rhocem, rhohyd, rhoH2O);
[vol_cem_ext, vol_hyd_ext, vol_wat_ext, vol_air_ext] = calculate_volumes(xi_plot_range, wc_extended, wa0a_extended, alpha_extended, ac_ratio, Fpor_plot, rhocem, rhohyd, rhoH2O);

% --- Create Figure 1 ---
figure('Name', 'Figure 12 Replication', 'Position', [100, 100, 1000, 450]);
% Plot (a): Normal Model
subplot(1, 2, 1);
Y_normal = [vol_cem_norm, vol_hyd_norm, vol_wat_norm, vol_air_norm];
h_area_norm = area(xi_plot_range, Y_normal, 'LineWidth', 1.5);
h_area_norm(1).FaceColor = [0.8 0.8 0.8]; h_area_norm(2).FaceColor = [0.9 0.9 0.9];
h_area_norm(3).FaceColor = [1.0 1.0 1.0]; h_area_norm(4).FaceColor = [1.0 1.0 1.0];
h_area_norm(1).EdgeColor = 'k'; h_area_norm(2).EdgeColor = 'k';
h_area_norm(3).EdgeColor = 'k'; h_area_norm(4).EdgeColor = [0.5 0.5 0.5]; h_area_norm(4).LineStyle = ':';
title('cement paste: w_{cp}(0)/c = 0.50, \alpha = 0.0', 'FontSize', 12);
xlabel('hydration degree [-]', 'FontSize', 12); ylabel('volume fractions [-]', 'FontSize', 12);
xlim([0 1]); ylim([0 1]); grid on; box on;
text(0.4, 0.1, 'cement', 'FontSize', 14, 'BackgroundColor', h_area_norm(1).FaceColor,'HorizontalAlignment','center');
text(0.4, 0.4, 'hydrates', 'FontSize', 14, 'BackgroundColor', h_area_norm(2).FaceColor,'HorizontalAlignment','center');
text(0.4, 0.8, 'water', 'FontSize', 14, 'BackgroundColor', h_area_norm(3).FaceColor,'HorizontalAlignment','center');
text(0.8, 0.9, 'air', 'FontSize', 14);

% Plot (b): Extended Model
subplot(1, 2, 2);
Y_extended = [vol_cem_ext, vol_hyd_ext, vol_wat_ext, vol_air_ext];
h_area_ext = area(xi_plot_range, Y_extended, 'LineWidth', 1.5);
h_area_ext(1).FaceColor = [0.8 0.8 0.8]; h_area_ext(2).FaceColor = [0.9 0.9 0.9];
h_area_ext(3).FaceColor = [1.0 1.0 1.0]; h_area_ext(4).FaceColor = [1.0 1.0 1.0];
h_area_ext(1).EdgeColor = 'k'; h_area_ext(2).EdgeColor = 'k';
h_area_ext(3).EdgeColor = 'k'; h_area_ext(4).EdgeColor = [0.5 0.5 0.5]; h_area_ext(4).LineStyle = ':';
title('cement paste: w_{cp}(0)/c = 0.4703, \alpha = 0.603', 'FontSize', 12);
xlabel('hydration degree [-]', 'FontSize', 12); ylabel('volume fractions [-]', 'FontSize', 12);
xlim([0 1]); ylim([0 1]); grid on; box on;
text(0.4, 0.1, 'cement', 'FontSize', 14, 'BackgroundColor', h_area_ext(1).FaceColor,'HorizontalAlignment','center');
text(0.4, 0.4, 'hydrates', 'FontSize', 14, 'BackgroundColor', h_area_ext(2).FaceColor,'HorizontalAlignment','center');
text(0.4, 0.8, 'water', 'FontSize', 14, 'BackgroundColor', h_area_ext(3).FaceColor,'HorizontalAlignment','center');
text(0.8, 0.9, 'air', 'FontSize', 14);

%% PLOT 2: DIAGNOSTIC COMPARISON PLOTS
disp('Generating Plot 2: Diagnostic Comparison Plots...');

% --- Find indices and extract data ---
idx_wc_norm = find(wc_list == wc_normal);
idx_wa0a_norm = find(wa0_a_list == wa0a_normal);
idx_alpha_norm = find(alpha_list == alpha_normal);
idx_Fpor_plot = find(Fpor_list == Fpor_plot);
idx_wc_ext = find(wc_list == wc_extended);
idx_wa0a_ext = find(wa0_a_list == wa0a_extended);
idx_alpha_ext = find(alpha_list == alpha_extended);

wceff_norm = cellfun(@(s) s.wceff, output_cell(idx_wc_norm, :, idx_Fpor_plot, idx_wa0a_norm, idx_alpha_norm));
E_cp_norm = cellfun(@(s) s.calc_cp.E, output_cell(idx_wc_norm, :, idx_Fpor_plot, idx_wa0a_norm, idx_alpha_norm));
wceff_ext = cellfun(@(s) s.wceff, output_cell(idx_wc_ext, :, idx_Fpor_plot, idx_wa0a_ext, idx_alpha_ext));
E_cp_ext = cellfun(@(s) s.calc_cp.E, output_cell(idx_wc_ext, :, idx_Fpor_plot, idx_wa0a_ext, idx_alpha_ext));

% --- Calculate Saturation Ratio (S) ---
[~, vh_n, vw_n, va_n] = calculate_volumes(xi_list, wc_normal, wa0a_normal, alpha_normal, ac_ratio, Fpor_plot, rhocem, rhohyd, rhoH2O);
gel_pores_norm = 0.28 * vh_n; total_pores_norm = vw_n + va_n + gel_pores_norm; S_norm = vw_n ./ total_pores_norm;
[~, vh_e, vw_e, va_e] = calculate_volumes(xi_list, wc_extended, wa0a_extended, alpha_extended, ac_ratio, Fpor_plot, rhocem, rhohyd, rhoH2O);
gel_pores_ext = 0.28 * vh_e; total_pores_ext = vw_e + va_e + gel_pores_ext; S_ext = vw_e ./ total_pores_ext;

% --- Create Figure 2 ---
figure('Name', 'Diagnostic Model Comparison', 'Position', [150, 150, 1400, 400]);
subplot(1, 3, 1);
plot(xi_list, wceff_norm, 'r-o', 'LineWidth', 2, 'MarkerFaceColor', 'r', 'DisplayName', 'Normal Model'); hold on;
plot(xi_list, wceff_ext, 'b-s', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'DisplayName', 'Extended Model'); hold off;
title('Effective w/c Ratio Evolution', 'FontSize', 14); xlabel('Hydration Degree \xi [-]', 'FontSize', 12);
ylabel('w_{cp}(\xi)/c [-]', 'FontSize', 12); legend('show', 'Location', 'best'); grid on; box on;
subplot(1, 3, 2);
plot(xi_list, E_cp_norm, 'r-o', 'LineWidth', 2, 'MarkerFaceColor', 'r', 'DisplayName', 'Normal Model'); hold on;
plot(xi_list, E_cp_ext, 'b-s', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'DisplayName', 'Extended Model'); hold off;
title('Cement Paste Stiffness (E_{cp})', 'FontSize', 14); xlabel('Hydration Degree \xi [-]', 'FontSize', 12);
ylabel('Young''s Modulus [GPa]', 'FontSize', 12); legend('show', 'Location', 'best'); grid on; box on;
subplot(1, 3, 3);
plot(xi_list, S_norm, 'r-o', 'LineWidth', 2, 'MarkerFaceColor', 'r', 'DisplayName', 'Self-Desiccation'); hold on;
plot(xi_list, S_ext, 'b-s', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'DisplayName', 'Internal Curing'); hold off;
title('Internal Saturation Ratio', 'FontSize', 14); xlabel('Hydration Degree \xi [-]', 'FontSize', 12);
ylabel('Saturation Ratio S [-]', 'FontSize', 12); ylim([0.7 1.05]); grid on; box on;
legend('show', 'Location', 'best');


%% PLOT 3: NEW INSIGHTFUL PLOT - COMPOSITION OF PORE SPACE
disp('Generating Plot 3: Composition of Pore Space...');

figure('Name', 'Insightful Plot: Pore Space Composition', 'Position', [200, 200, 1000, 450]);

% Data for this plot comes from the vol_wat and vol_air arrays already calculated
porosity_data_norm = [vol_wat_norm, vol_air_norm];
porosity_data_ext = [vol_wat_ext, vol_air_ext];

% Plot (a): Normal Model Pore Space
subplot(1, 2, 1);
h_area_por_norm = area(xi_plot_range, porosity_data_norm, 'LineWidth', 1.5, 'EdgeColor', 'k');
h_area_por_norm(1).FaceColor = [0.6 0.8 1.0]; % Water-filled (light blue)
h_area_por_norm(2).FaceColor = [0.9 0.9 0.9]; % Air-filled (light gray)
title('Pore Composition: Normal Model', 'FontSize', 14);
xlabel('Hydration Degree \xi [-]', 'FontSize', 12);
ylabel('Porosity Volume Fraction [-]', 'FontSize', 12);
legend({'Water-filled', 'Air-filled'}, 'Location', 'best');
grid on; box on;
max_y = max(sum(porosity_data_norm, 2)); % Find max porosity for consistent y-axis

% Plot (b): Extended Model Pore Space
subplot(1, 2, 2);
h_area_por_ext = area(xi_plot_range, porosity_data_ext, 'LineWidth', 1.5, 'EdgeColor', 'k');
h_area_por_ext(1).FaceColor = [0.6 0.8 1.0]; % Water-filled (light blue)
h_area_por_ext(2).FaceColor = [0.9 0.9 0.9]; % Air-filled (light gray)
title('Pore Composition: Extended Model', 'FontSize', 14);
xlabel('Hydration Degree \xi [-]', 'FontSize', 12);
ylabel('Porosity Volume Fraction [-]', 'FontSize', 12);
legend({'Water-filled', 'Air-filled'}, 'Location', 'best');
grid on; box on;
max_y = max([max_y, max(sum(porosity_data_ext, 2))]); % Get overall max porosity

% Set uniform y-axis for both subplots for fair comparison
subplot(1, 2, 1); ylim([0 max_y*1.1]);
subplot(1, 2, 2); ylim([0 max_y*1.1]);


%% --- LOCAL FUNCTION DEFINITION ---
function [vol_cem, vol_hyd, vol_water, vol_air] = calculate_volumes(xi_range, wc, wa0a, alpha, ac, Fpor, rhoC, rhoH, rhoW)
    n_points = length(xi_range);
    vol_cem = zeros(n_points, 1);
    vol_hyd = zeros(n_points, 1);
    vol_water = zeros(n_points, 1);
    vol_air = zeros(n_points, 1);
    
    for i = 1:n_points
        xi = xi_range(i);
        
        % Extended model logic
        wc0_eff = wc - wa0a * ac;
        if wc0_eff <= 0, wc0_eff = 1e-6; end
        
        k_factor = (1.051 + 3.31 * wc0_eff) / (20 + 63 * wc0_eff);
        wceff = wc0_eff + k_factor * alpha * xi;

        % Unscaled volumes
        denominator = (1 + rhoC/rhoW * wceff);
        fcem_unscaled = (1 - xi) / denominator;
        fhyd_unscaled = (1.42 * rhoC * xi) / (rhoH * denominator);
        fH2O_cap_unscaled = max(0, (rhoC*(wceff - 0.42*xi))/(rhoW*denominator));
        f_voids_unscaled = max(0, 1 - fcem_unscaled - fhyd_unscaled - fH2O_cap_unscaled);

        % Apply Fpor scaling
        f_solids = fcem_unscaled + fhyd_unscaled;
        f_porosity_unscaled = fH2O_cap_unscaled + f_voids_unscaled;
        f_porosity_scaled = f_porosity_unscaled * Fpor;
        
        if f_solids < 1e-9, renorm_factor = 0; else, renorm_factor = (1 - f_porosity_scaled) / f_solids; end
        if renorm_factor < 0, continue; end
        
        fcem_total = fcem_unscaled * renorm_factor;
        fhyd_total = fhyd_unscaled * renorm_factor;
        
        if f_porosity_unscaled < 1e-9, por_scale_factor = 0; else, por_scale_factor = f_porosity_scaled / f_porosity_unscaled; end
        f_voids_scaled = f_voids_unscaled * por_scale_factor;
        fH2O_cap_scaled = fH2O_cap_unscaled * por_scale_factor;
        
        % Final volumes using Eqs. (44) and (45)
        vol_cem(i) = fcem_total;
        vol_hyd(i) = fhyd_total;
        vol_water(i) = fH2O_cap_scaled + alpha * f_voids_scaled;
        vol_air(i) = (1 - alpha) * f_voids_scaled;
    end
end
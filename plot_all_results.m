%% UNIFIED SCRIPT TO PLOT AND COMPARE MODEL RESULTS
% This script can analyze data from calculations WITH or WITHOUT an Fpor loop.
% It will automatically detect the data structure and adapt.
clear all; clc; close all;

%% --- 1. USER INTERACTION: SELECT FILE AND VERIFY PARAMETERS ---

% --- Let user select the .mat file to analyze ---
[filename, pathname] = uigetfile('*.mat', 'Select the output .mat file');
if isequal(filename,0)
    disp('User selected Cancel');
    return; % Stop the script if no file is chosen
else
    disp(['User selected: ', fullfile(pathname, filename)]);
    load(fullfile(pathname, filename));
end

% --- USER: VERIFY THESE LISTS MATCH YOUR CALCULATION ---
% IMPORTANT: These lists must exactly match the ones used to generate the
% .mat file you just selected.
wc_list = [0.50, 0.55, 0.60];
xi_list = 0.05:0.05:1.0;
Fpor_list = 1.0:0.1:1.7; % This is only used if the 5D file is loaded
wa0_a_list = [0.0, 0.005, 0.01];
alpha_list = [0.0, 0.6];
ac_ratio = 4.9;

% --- Material Densities ---
rhoH2O = 1000; rhohyd = 2073; rhocem = 3150;

%% --- 2. DETECT DATA STRUCTURE AND EXTRACT DATA ---

% Define the specific cases we want to plot for comparison
wc_normal = 0.55; wa0a_normal = 0.0; alpha_normal = 0.0;
wc_extended = 0.55; wa0a_extended = 0.01; alpha_extended = 0.6;
Fpor_plot = 1.0; % We will always compare cases at Fpor=1.0

% Check if the required plot parameters exist in the lists
if ~ismember(wc_normal, wc_list) || ~ismember(wc_extended, wc_list) || ...
   ~ismember(wa0a_normal, wa0_a_list) || ~ismember(wa0a_extended, wa0_a_list) || ...
   ~ismember(alpha_normal, alpha_list) || ~ismember(alpha_extended, alpha_list)
    error('One of the plot parameters (wc, wa0a, alpha) is not in the defined lists.');
end

% --- Find indices for the plot cases ---
idx_wc_norm = find(wc_list == wc_normal);
idx_wa0a_norm = find(wa0_a_list == wa0a_normal);
idx_alpha_norm = find(alpha_list == alpha_normal);

idx_wc_ext = find(wc_list == wc_extended);
idx_wa0a_ext = find(wa0_a_list == wa0a_extended);
idx_alpha_ext = find(alpha_list == alpha_extended);

% --- Intelligent Data Extraction based on dimensions ---
if ndims(output_cell) == 5 % Case WITH Fpor
    disp('Detected 5D data structure (with Fpor).');
    if ~ismember(Fpor_plot, Fpor_list)
        error('Fpor_plot=%.2f is not in the 5D data''s Fpor_list.', Fpor_plot);
    end
    idx_Fpor_plot = find(Fpor_list == Fpor_plot);
    
    % Extract using 5D indexing: output_cell(wc, xi, Fpor, wa0a, alpha)
    wceff_norm = cellfun(@(s) s.wceff, output_cell(idx_wc_norm, :, idx_Fpor_plot, idx_wa0a_norm, idx_alpha_norm));
    E_cp_norm = cellfun(@(s) s.calc_cp.E, output_cell(idx_wc_norm, :, idx_Fpor_plot, idx_wa0a_norm, idx_alpha_norm));
    wceff_ext = cellfun(@(s) s.wceff, output_cell(idx_wc_ext, :, idx_Fpor_plot, idx_wa0a_ext, idx_alpha_ext));
    E_cp_ext = cellfun(@(s) s.calc_cp.E, output_cell(idx_wc_ext, :, idx_Fpor_plot, idx_wa0a_ext, idx_alpha_ext));
    
elseif ndims(output_cell) == 4 % Case WITHOUT Fpor
    disp('Detected 4D data structure (without Fpor).');
    
    % Extract using 4D indexing: output_cell(wc, xi, wa0a, alpha)
    wceff_norm = cellfun(@(s) s.wceff, output_cell(idx_wc_norm, :, idx_wa0a_norm, idx_alpha_norm));
    E_cp_norm = cellfun(@(s) s.calc_cp.E, output_cell(idx_wc_norm, :, idx_wa0a_norm, idx_alpha_norm));
    wceff_ext = cellfun(@(s) s.wceff, output_cell(idx_wc_ext, :, idx_wa0a_ext, idx_alpha_ext));
    E_cp_ext = cellfun(@(s) s.calc_cp.E, output_cell(idx_wc_ext, :, idx_wa0a_ext, idx_alpha_ext));

else
    error('Unsupported data structure. Expected 4 or 5 dimensions, but found %d.', ndims(output_cell));
end


%% --- 3. GENERATE ALL PLOTS ---

% The rest of the script is now independent of the source file structure,
% as it uses the correctly extracted data vectors.

% --- PLOT 1: REPLICATE JOURNAL FIGURE 12 ---
disp('Generating Plot 1: Replication of Journal Figure 12...');
xi_plot_range = [0, xi_list];
[vol_cem_norm, vol_hyd_norm, vol_wat_norm, vol_air_norm] = calculate_volumes(xi_plot_range, wc_normal, wa0a_normal, alpha_normal, ac_ratio, Fpor_plot, rhocem, rhohyd, rhoH2O);
[vol_cem_ext, vol_hyd_ext, vol_wat_ext, vol_air_ext] = calculate_volumes(xi_plot_range, wc_extended, wa0a_extended, alpha_extended, ac_ratio, Fpor_plot, rhocem, rhohyd, rhoH2O);
figure('Name', 'Figure 12 Replication', 'Position', [100, 100, 1000, 450]);
% Plot (a)
subplot(1, 2, 1);
Y_normal = [vol_cem_norm, vol_hyd_norm, vol_wat_norm, vol_air_norm];
h_area_norm = area(xi_plot_range, Y_normal, 'LineWidth', 1.5);
h_area_norm(1).FaceColor=[.8 .8 .8];h_area_norm(2).FaceColor=[.9 .9 .9];h_area_norm(3).FaceColor=[1 1 1];h_area_norm(4).FaceColor=[1 1 1];
h_area_norm(1).EdgeColor='k';h_area_norm(2).EdgeColor='k';h_area_norm(3).EdgeColor='k';h_area_norm(4).EdgeColor=[.5 .5 .5];h_area_norm(4).LineStyle=':';
title('cement paste: w_{cp}(0)/c = 0.50, \alpha = 0.0', 'FontSize', 12);
xlabel('hydration degree [-]'); ylabel('volume fractions [-]'); xlim([0 1]); ylim([0 1]); grid on; box on;
text(0.4,0.1,'cement','FontSize',14,'BackgroundColor',h_area_norm(1).FaceColor,'HorizontalAlignment','center');
text(0.4,0.4,'hydrates','FontSize',14,'BackgroundColor',h_area_norm(2).FaceColor,'HorizontalAlignment','center');
text(0.4,0.8,'water','FontSize',14,'BackgroundColor',h_area_norm(3).FaceColor,'HorizontalAlignment','center');
text(0.8,0.9,'air','FontSize',14);
% Plot (b)
subplot(1, 2, 2);
Y_extended = [vol_cem_ext, vol_hyd_ext, vol_wat_ext, vol_air_ext];
h_area_ext = area(xi_plot_range, Y_extended, 'LineWidth', 1.5);
h_area_ext(1).FaceColor=[.8 .8 .8];h_area_ext(2).FaceColor=[.9 .9 .9];h_area_ext(3).FaceColor=[1 1 1];h_area_ext(4).FaceColor=[1 1 1];
h_area_ext(1).EdgeColor='k';h_area_ext(2).EdgeColor='k';h_area_ext(3).EdgeColor='k';h_area_ext(4).EdgeColor=[.5 .5 .5];h_area_ext(4).LineStyle=':';
title('cement paste: w_{cp}(0)/c = 0.4703, \alpha = 0.603', 'FontSize', 12);
xlabel('hydration degree [-]'); ylabel('volume fractions [-]'); xlim([0 1]); ylim([0 1]); grid on; box on;
text(0.4,0.1,'cement','FontSize',14,'BackgroundColor',h_area_ext(1).FaceColor,'HorizontalAlignment','center');
text(0.4,0.4,'hydrates','FontSize',14,'BackgroundColor',h_area_ext(2).FaceColor,'HorizontalAlignment','center');
text(0.4,0.8,'water','FontSize',14,'BackgroundColor',h_area_ext(3).FaceColor,'HorizontalAlignment','center');
text(0.8,0.9,'air','FontSize',14);

% --- PLOT 2: DIAGNOSTIC COMPARISON PLOTS ---
disp('Generating Plot 2: Diagnostic Comparison Plots...');
[~, vh_n, vw_n, va_n] = calculate_volumes(xi_list, wc_normal, wa0a_normal, alpha_normal, ac_ratio, Fpor_plot, rhocem, rhohyd, rhoH2O);
gel_pores_norm = 0.28 * vh_n; total_pores_norm = vw_n + va_n + gel_pores_norm; S_norm = vw_n ./ total_pores_norm;
[~, vh_e, vw_e, va_e] = calculate_volumes(xi_list, wc_extended, wa0a_extended, alpha_extended, ac_ratio, Fpor_plot, rhocem, rhohyd, rhoH2O);
gel_pores_ext = 0.28 * vh_e; total_pores_ext = vw_e + va_e + gel_pores_ext; S_ext = vw_e ./ total_pores_ext;
figure('Name', 'Diagnostic Model Comparison', 'Position', [150, 150, 1400, 400]);
subplot(1,3,1); plot(xi_list,wceff_norm,'r-o','LineWidth',2,'MarkerFaceColor','r','DisplayName','Normal Model'); hold on; plot(xi_list,wceff_ext,'b-s','LineWidth',2,'MarkerFaceColor','b','DisplayName','Extended Model'); hold off; title('Effective w/c Ratio Evolution'); xlabel('Hydration Degree \xi [-]'); ylabel('w_{cp}(\xi)/c [-]'); legend('show','Location','best'); grid on; box on;
subplot(1,3,2); plot(xi_list,E_cp_norm,'r-o','LineWidth',2,'MarkerFaceColor','r','DisplayName','Normal Model'); hold on; plot(xi_list,E_cp_ext,'b-s','LineWidth',2,'MarkerFaceColor','b','DisplayName','Extended Model'); hold off; title('Cement Paste Stiffness (E_{cp})'); xlabel('Hydration Degree \xi [-]'); ylabel('Young''s Modulus [GPa]'); legend('show','Location','best'); grid on; box on;
subplot(1,3,3); plot(xi_list,S_norm,'r-o','LineWidth',2,'MarkerFaceColor','r','DisplayName','Self-Desiccation'); hold on; plot(xi_list,S_ext,'b-s','LineWidth',2,'MarkerFaceColor','b','DisplayName','Internal Curing'); hold off; title('Internal Saturation Ratio'); xlabel('Hydration Degree \xi [-]'); ylabel('Saturation Ratio S [-]'); ylim([0.7 1.05]); grid on; box on; legend('show','Location','best');

% --- PLOT 3: PORE SPACE COMPOSITION ---
disp('Generating Plot 3: Composition of Pore Space...');
figure('Name', 'Insightful Plot: Pore Space Composition', 'Position', [200, 200, 1000, 450]);
porosity_data_norm = [vol_wat_norm, vol_air_norm];
porosity_data_ext = [vol_wat_ext, vol_air_ext];
% Plot (a)
subplot(1, 2, 1);
h_area_por_norm = area(xi_plot_range, porosity_data_norm, 'LineWidth', 1.5, 'EdgeColor', 'k');
h_area_por_norm(1).FaceColor=[0.6 0.8 1.0]; h_area_por_norm(2).FaceColor=[0.9 0.9 0.9];
title('Pore Composition: Normal Model'); xlabel('Hydration Degree \xi [-]'); ylabel('Porosity Volume Fraction [-]');
legend({'Water-filled', 'Air-filled'},'Location','best'); grid on; box on;
max_y = max(sum(porosity_data_norm,2));
% Plot (b)
subplot(1, 2, 2);
h_area_por_ext = area(xi_plot_range, porosity_data_ext, 'LineWidth', 1.5, 'EdgeColor', 'k');
h_area_por_ext(1).FaceColor=[0.6 0.8 1.0]; h_area_por_ext(2).FaceColor=[0.9 0.9 0.9];
title('Pore Composition: Extended Model'); xlabel('Hydration Degree \xi [-]'); ylabel('Porosity Volume Fraction [-]');
legend({'Water-filled', 'Air-filled'},'Location','best'); grid on; box on;
max_y = max([max_y, max(sum(porosity_data_ext,2))]);
% Set uniform y-axis
subplot(1, 2, 1); ylim([0 max_y*1.1]);
subplot(1, 2, 2); ylim([0 max_y*1.1]);


%% --- LOCAL FUNCTION DEFINITION ---
function [vol_cem, vol_hyd, vol_water, vol_air] = calculate_volumes(xi_range, wc, wa0a, alpha, ac, Fpor, rhoC, rhoH, rhoW)
    n_points = length(xi_range);
    vol_cem=zeros(n_points,1); vol_hyd=zeros(n_points,1); vol_water=zeros(n_points,1); vol_air=zeros(n_points,1);
    for i = 1:n_points
        xi = xi_range(i);
        wc0_eff = wc - wa0a * ac;
        if wc0_eff <= 0, wc0_eff = 1e-6; end
        k_factor = (1.051 + 3.31 * wc0_eff) / (20 + 63 * wc0_eff);
        wceff = wc0_eff + k_factor * alpha * xi;
        denominator = (1 + rhoC/rhoW * wceff);
        fcem_unscaled = (1 - xi) / denominator;
        fhyd_unscaled = (1.42 * rhoC * xi) / (rhoH * denominator);
        fH2O_cap_unscaled = max(0, (rhoC*(wceff-0.42*xi))/(rhoW*denominator));
        f_voids_unscaled = max(0, 1 - fcem_unscaled - fhyd_unscaled - fH2O_cap_unscaled);
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
        vol_cem(i) = fcem_total;
        vol_hyd(i) = fhyd_total;
        vol_water(i) = fH2O_cap_scaled + alpha * f_voids_scaled;
        vol_air(i) = (1 - alpha) * f_voids_scaled;
    end
end
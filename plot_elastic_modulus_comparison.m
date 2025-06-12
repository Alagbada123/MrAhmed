% plot_elastic_modulus_comparison.m
clear; close all; clc;

%% 1) Load modeled data
% elastic‐only results
S_od  = load('RC_final_Eod141312_sensFpor.mat','RC_sens');
S_ssd = load('RC_final_Essd141312_sensFpor.mat','RC_sens');
% extract Econc for wc=0.55, Fpor index=1
E_od  = cell2mat( squeeze( S_od.RC_sens(1,:,1) ) );
E_ssd = cell2mat( squeeze( S_ssd.RC_sens(1,:,1) ) );

%% 2) Define hydration degrees
xi_values = 0.05:0.05:1;
xi_plot   = [0, xi_values];
Eplot_od  = [0, E_od];
Eplot_ssd = [0, E_ssd];

%% 3) Load experimental data
run('create_experimental_data.m');    % defines experimental_data
row = 4;  % RCA sample row
exp_xi  = experimental_data.xi_marker_vals;
exp_E   = experimental_data.SSD.youngs_modulus(row,:);
exp_SE  = experimental_data.SSD.youngs_modulus_stddev(row,:);

%% 4) Plot
figure; hold on;
plot(xi_plot, Eplot_od,  '-o','LineWidth',1.5,'DisplayName','Modeled E\_OD');
plot(xi_plot, Eplot_ssd, '-s','LineWidth',1.5,'DisplayName','Modeled E\_SSD');
errorbar(exp_xi, exp_E, exp_SE, '^','MarkerSize',6,'LineWidth',1.2, ...
    'DisplayName','Experimental RCA');

% formatting
grid on; box on;
set(gca,'FontName','Times New Roman','FontSize',12,'LineWidth',1);
ax = gca;
ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
ax.XMinorGrid = 'on'; ax.YMinorGrid = 'on';
ax.GridLineStyle = ':'; ax.MinorGridLineStyle = ':';
ax.GridAlpha = 0.4;   ax.MinorGridAlpha = 0.4;
ax.XAxis.MinorTickValues = 0.1:0.2:0.9;
ax.YAxis.MinorTickValues = 5:10:45;

xlabel('Hydration Degree \xi [-]','FontSize',12);
ylabel('E_{conc} [GPa]','FontSize',12,'Interpreter','tex');
xlim([0 1]); ylim([0 50]);
legend('Location','southeast','FontSize',10);
hold off;

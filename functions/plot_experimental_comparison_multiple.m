% plot_experimental_comparison_multiple.m
function plot_experimental_comparison_multiple(analysisParams, E_results, material_names, experimental_data)
    % Create figure with specific size
    figure('Name', 'Elasticity of Concrete With and Without Coating for Aggregate', ...
           'Position', [100, 100, 1000, 600]);
    hold on;
    
    % Set up grid and appearance
    grid on;
    box on;
    set(gca, 'FontSize', 12, 'LineWidth', 1);
    
    % Define mapping between material cases and experimental data rows
    exp_data_mapping = struct(...
        'NAC_ITZ', 1, ...    % NAC is row 1
        'RCA', 4, ...        % RC-SSD is row 4
        'RCA_ITZ', 4);       % RC-SSD is row 4 (same as RCA)
    
    % Plot numerical results with specific styles
    plot(analysisParams.xiValues, E_results(1, :), 'b-', 'LineWidth', 2, ...
        'DisplayName', 'Without Coating');
    plot(analysisParams.xiValues, E_results(2, :), 'r--', 'LineWidth', 2, ...
        'DisplayName', 'With Coating');
    
    % Define colors and markers for experimental data
    exp_colors = [...
        1.00, 0.65, 0.00;  % Orange for NAC
        0.49, 0.18, 0.56;  % Purple for RPB-OD
        0.47, 0.67, 0.19;  % Green for RSB-OD
        0.30, 0.75, 0.93;  % Light blue for RC-OD
        0.85, 0.33, 0.10;  % Red-orange for RCa-OD
        0.00, 0.45, 0.74;  % Blue for RCc-OD
        0.64, 0.08, 0.18]; % Dark red for RCh-OD
    
    markers = {'o', 's', 'd', 'v', '^', '>', '<'};
    
    % Plot all experimental datasets
    for i = 1:7
        errorbar(experimental_data.xi_marker_vals, ...
                experimental_data.OD.youngs_modulus(i,:), ...
                experimental_data.OD.standard_deviation(i,:), ...
                'Color', exp_colors(i,:), ...
                'Marker', markers{i}, ...
                'MarkerFaceColor', exp_colors(i,:), ...
                'MarkerSize', 6, ...
                'LineStyle', '-', ...
                'LineWidth', 1, ...
                'CapSize', 5, ...
                'DisplayName', experimental_data.OD.sample_descriptions{i});
    end
    
    % Set axis labels
    xlabel('Hydration Degree ξ [-]', 'FontSize', 12);
    ylabel('E-Modulus E_{conc} [GPa]', 'FontSize', 12);
    title('Elasticity of Concrete With and Without Coating for Aggregate', ...
          'FontSize', 14);
    
    % Set axis limits
    xlim([0 1]);
    ylim([0 50]);
    
    % Create legend
    legend('Location', 'eastoutside', 'FontSize', 10);
    
    % Adjust plot margins to prevent cutting off
    set(gca, 'Position', [0.1 0.1 0.7 0.8]);
    
    hold off;
end
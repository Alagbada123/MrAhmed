function plot_experimental_comparison(analysisParams, E_with_ITZ, E_without_ITZ, experimental_data, plot_type)
    % plot_type: 'single' for one material, 'comparison' for two materials
    
    % Create two subplots for OD and SSD data
    figure('Position', [100, 100, 1200, 1200]);

    % Define markers and colors for experimental data
    markers = {'o', 's', 'd', '^', 'v', '>', '<'};
    colors = {[0 0 0], [0.5 0 0.5], [0.4660 0.6740 0.1880], [1 0 0], ...
              [0.8500 0.3250 0.0980], [0 0.7 1], [1 0.5 0]};

    %% First subplot: OD Data
    subplot(2,1,1);
    
    % Plot theoretical curves based on plot_type
    if strcmp(plot_type, 'comparison')
        plot(analysisParams.xiValues, E_without_ITZ, 'b-', 'LineWidth', 2, 'DisplayName', 'Without Coating');
        hold on;
        plot(analysisParams.xiValues, E_with_ITZ, 'r--', 'LineWidth', 2, 'DisplayName', 'With Coating');
    else
        plot(analysisParams.xiValues, E_with_ITZ, 'r--', 'LineWidth', 2, 'DisplayName', 'With Coating');
    end
    hold on;

    % Plot OD experimental data
    for i = 1:size(experimental_data.OD.youngs_modulus, 1)
        errorbar(experimental_data.xi_marker_vals, ...
            experimental_data.OD.youngs_modulus(i,:), ...
            experimental_data.OD.standard_deviation(i,:), ...
            markers{i}, 'Color', colors{i}, ...
            'MarkerFaceColor', colors{i}, ...
            'LineWidth', 1.5, ...
            'DisplayName', strrep(experimental_data.OD.sample_descriptions{i}, ' - ', ' '));
    end

    % Format OD plot
    xlabel('Hydration Degree \xi [-]', 'FontSize', 12);
    ylabel('E-Modulus E_{conc} [GPa]', 'FontSize', 12);
    title('Elasticity of Concrete (Oven-Dried Samples)', 'FontSize', 14);
    format_plot_properties();

    %% Second subplot: SSD Data
    subplot(2,1,2);
    
    % Plot theoretical curves based on plot_type
    if strcmp(plot_type, 'comparison')
        plot(analysisParams.xiValues, E_without_ITZ, 'b-', 'LineWidth', 2, 'DisplayName', 'Without Coating');
        hold on;
        plot(analysisParams.xiValues, E_with_ITZ, 'r--', 'LineWidth', 2, 'DisplayName', 'With Coating');
    else
        plot(analysisParams.xiValues, E_with_ITZ, 'r--', 'LineWidth', 2, 'DisplayName', 'With Coating');
    end
    hold on;

    % Plot SSD experimental data
    for i = 1:size(experimental_data.SSD.youngs_modulus, 1)
        errorbar(experimental_data.xi_marker_vals, ...
            experimental_data.SSD.youngs_modulus(i,:), ...
            experimental_data.SSD.standard_deviation(i,:), ...
            markers{i}, 'Color', colors{i}, ...
            'MarkerFaceColor', colors{i}, ...
            'LineWidth', 1.5, ...
            'DisplayName', strrep(experimental_data.SSD.sample_descriptions{i}, ' - ', ' '));
    end

    % Format SSD plot
    xlabel('Hydration Degree \xi [-]', 'FontSize', 12);
    ylabel('E-Modulus E_{conc} [GPa]', 'FontSize', 12);
    title('Elasticity of Concrete (Saturated Surface-Dried Samples)', 'FontSize', 14);
    format_plot_properties();

    % Adjust spacing between subplots
    set(gcf, 'Color', 'white');
    set(gcf, 'Position', [100, 100, 1200, 1200]);
end

function format_plot_properties()
    % Helper function to format plot properties
    grid on;
    legend('Location', 'southeast', 'FontSize', 10);
    set(gca, 'FontSize', 10);
    
    % Set axis limits and ticks
    xlim([0 1]);
    ylim([0 50]);
    xticks(0:0.1:1);
    yticks(0:5:50);
    
    % Make grid more visible
    set(gca, 'GridAlpha', 0.2);
    set(gca, 'MinorGridAlpha', 0.1);
    grid minor;
    
    % Adjust figure properties
    box on;
    set(gca, 'LineWidth', 1);
end
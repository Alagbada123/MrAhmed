function plot_basic_comparison(analysisParams, E_with_ITZ, E_without_ITZ, Aggregate_with_ITZ)
    figure('Position', [100, 100, 800, 600]);
    plot(analysisParams.xiValues, E_with_ITZ, 'r--', 'LineWidth', 2, 'DisplayName', 'With ITZ');
    hold on;
    plot(analysisParams.xiValues, E_without_ITZ, 'b-', 'LineWidth', 2, 'DisplayName', 'Without ITZ');
    hold off;

    % Format plot
    xlabel('Degree of Hydration \xi', 'FontSize', 14);
    ylabel('Young''s Modulus (GPa)', 'FontSize', 14);
    title('Concrete Stiffness Evolution: Effect of ITZ', 'FontSize', 16);
    grid on;
    legend('Location', 'southeast', 'FontSize', 12);
    set(gca, 'FontSize', 12);

    % Set axis limits explicitly
    xlim([0 1]);
    xticks(0:0.1:1);

    % Add mix design information
    text_str = sprintf('Mix Design:\nw/c = %.2f\nsand/c = %.3f\ngravel/c = %.3f', ...
        Aggregate_with_ITZ.mixDesign.waterCementRatio, ...
        Aggregate_with_ITZ.mixDesign.aggregateCementRatios.sand, ...
        Aggregate_with_ITZ.mixDesign.coatedAggregateCementRatios.gravel);
    annotation('textbox', [0.15, 0.7, 0.2, 0.2], 'String', text_str, ...
        'FitBoxToText', 'on', 'BackgroundColor', 'white');
end
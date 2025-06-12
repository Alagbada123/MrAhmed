function plotResults(xi_values, results)
    % Create plots of mechanical properties vs hydration degree
    figure('Name', 'Concrete Properties vs. Hydration Degree');
    
    subplot(3,1,1);
    plot(xi_values, results.bulkModulus);
    title('Bulk Modulus vs. Hydration Degree');
    xlabel('Hydration Degree');
    ylabel('Bulk Modulus [GPa]');
    
    subplot(3,1,2);
    plot(xi_values, results.shearModulus);
    title('Shear Modulus vs. Hydration Degree');
    xlabel('Hydration Degree');
    ylabel('Shear Modulus [GPa]');
    
    subplot(3,1,3);
    plot(xi_values, results.youngsModulus);
    title('Young''s Modulus vs. Hydration Degree');
    xlabel('Hydration Degree');
    ylabel('Young''s Modulus [GPa]');
end
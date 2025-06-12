function plot_strength_envelopes(SI_ult, phi_hyd, c_hyd, stepsangle)
    % Create figure with multiple subplots for different stress spaces
    figure('Position', [100, 100, 1200, 800]);
    
    %% 1. Principal Stress Space (σ1-σ3)
    subplot(2,2,1)
    
    % Generate points for strength envelope
    sigma3 = linspace(-SI_ult*1000, 0, 100); % Compressive stress (MPa)
    sigma1_MC = zeros(size(sigma3));
    
    % Mohr-Coulomb envelope
    for i = 1:length(sigma3)
        sigma1_MC(i) = sigma3(i)*(1-sin(phi_hyd))/(1+sin(phi_hyd)) + ...
                       2*c_hyd*1000*cos(phi_hyd)/(1+sin(phi_hyd));
    end
    
    % Plot envelope
    plot(sigma3, sigma1_MC, 'b-', 'LineWidth', 2);
    hold on;
    plot([0 -SI_ult*1000], [0 -SI_ult*1000], 'k--'); % Hydrostatic line
    
    % Format plot
    xlabel('\sigma_3 [MPa]');
    ylabel('\sigma_1 [MPa]');
    title('Principal Stress Space (\sigma_1-\sigma_3)');
    grid on;
    axis equal;
    
    %% 2. Meridian Plot (p-q space)
    subplot(2,2,2)
    
    % Generate points for p-q space
    p = linspace(-SI_ult*1000/3, 0, 100); % Mean stress
    q_MC = zeros(size(p));
    
    % Convert MC parameters to p-q space
    alpha_pq = 6*sin(phi_hyd)/(3-sin(phi_hyd));
    k_pq = 6*c_hyd*1000*cos(phi_hyd)/(3-sin(phi_hyd));
    
    % Calculate q values
    for i = 1:length(p)
        q_MC(i) = -alpha_pq*p(i) + k_pq;
    end
    
    % Plot envelope
    plot(p, q_MC, 'r-', 'LineWidth', 2);
    
    % Format plot
    xlabel('p [MPa]');
    ylabel('q [MPa]');
    title('p-q Space');
    grid on;
    
    %% 3. Octahedral Plane
    subplot(2,2,3)
    
    % Generate points for octahedral plane
    theta = linspace(0, 2*pi, 100);
    r_oct = zeros(size(theta));
    p_fixed = -SI_ult*1000/3; % Fixed mean stress level
    
    % Calculate octahedral radius
    for i = 1:length(theta)
        cos3theta = cos(3*theta(i));
        r_oct(i) = (2*sqrt(2)/3)*(-alpha_pq*p_fixed + k_pq)/(1 + alpha_pq*cos3theta);
    end
    
    % Plot envelope
    polarplot(theta, r_oct, 'g-', 'LineWidth', 2);
    
    % Format plot
    title(['Octahedral Plane at p = ' num2str(p_fixed, '%.1f') ' MPa']);
    
    %% 4. 3D Strength Surface
    subplot(2,2,4)
    
    % Generate 3D surface points
    [theta_3D, rho_3D] = meshgrid(linspace(0, 2*pi, 50), ...
                                 linspace(0, SI_ult*1000, 50));
    X = zeros(size(theta_3D));
    Y = zeros(size(theta_3D));
    Z = zeros(size(theta_3D));
    
    for i = 1:size(rho_3D, 1)
        for j = 1:size(theta_3D, 2)
            r = rho_3D(i,j);
            t = theta_3D(i,j);
            
            % Convert to Cartesian coordinates
            X(i,j) = r*cos(t);
            Y(i,j) = r*sin(t);
            Z(i,j) = -alpha_pq*(X(i,j) + Y(i,j))/sqrt(6) + k_pq;
        end
    end
    
    % Plot 3D surface
    surf(X, Y, Z);
    colormap('winter');
    shading interp;
    
    % Format plot
    xlabel('\sigma_1 [MPa]');
    ylabel('\sigma_2 [MPa]');
    zlabel('\sigma_3 [MPa]');
    title('3D Strength Surface');
    grid on;
    view(45, 30);
    
    % Add overall title
    sgtitle(['Concrete Strength Envelopes (f_c = ' num2str(SI_ult*1000, '%.1f') ' MPa)'], ...
            'FontSize', 14, 'FontWeight', 'bold');
    
    % Adjust layout
    set(gcf, 'Color', 'white');
end

% % implementation in the main code After strength calculations
% if choose_Eonly == 0
%     % Plot strength envelopes
%     plot_strength_envelopes(SI_ult, phi_hyd, c_hyd, stepsangle);
% end

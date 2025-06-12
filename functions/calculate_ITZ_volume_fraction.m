function V_ITZ = calculate_ITZ_volume_fraction(Va, Dmin, Dmax, t_ITZ, A)
    % Function to calculate the ITZ volume fraction using Fuller's gradation
    % Inputs:
    % Va    : Volume fraction of aggregate
    % Dmin  : Minimum particle diameter (mm)
    % Dmax  : Maximum particle diameter (mm)
    % t_ITZ : Thickness of the ITZ (mm)
    % A     : Coefficient for the cubic term in g (0, 2, or 3)
    %
    % Output:
    % V_ITZ : ITZ volume fraction
    % 
    % Reference: doi: 10.1007/s11431-011-4737-x
    
    % Validate A
    if ~ismember(A, [0, 2, 3])
        error('A must be 0, 2, or 3');
    end
    
    % Calculate moments of the aggregate size distribution
    R1 = (5 * (Dmin * Dmax^2.5 - Dmin^2.5 * Dmax)) / (6 * (Dmax^2.5 - Dmin^2.5));
    R2 = (5 * (Dmin^2 * Dmax^2.5 - Dmin^2.5 * Dmax^2)) / (4 * (Dmax^2.5 - Dmin^2.5));
    R3 = (5 * (Dmin^2.5 * Dmax^3 - Dmin^3 * Dmax^2.5)) / (8 * (Dmax^2.5 - Dmin^2.5));
    
    % Calculate the number density of particles per unit volume
    Nv = (3 * Va) / (4 * pi * R3);
    
    % Calculate the coefficients c, d, g
    c = (4 * R2) / (1 - Va);
    d = (4 * R1 / (1 - Va)) + ((8 * pi * Nv * R2^2) / (1 - Va)^2);
    g = (4 * R1 / (3 * (1 - Va))) + ((16 * pi * Nv * R2^2 * R1) / (3 * (1 - Va)^2)) + ...
        ((64 * A * pi^2 * Nv^2 * R2^3) / (27 * (1 - Va)^3));
    
    % Calculate the ITZ volume fraction
    V_ITZ = 1 - Va - (1 - Va) * exp(-pi * Nv * (c * t_ITZ + d * t_ITZ^2 + g * t_ITZ^3));
end

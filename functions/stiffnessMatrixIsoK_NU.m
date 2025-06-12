function C = stiffnessMatrixIsoK_NU(k, mu)
    % Function to calculate the stiffness matrix of an isotropic material
    % based on Young's modulus (E) and Poisson's ratio (NU)
    % 
    % Inputs:
    % k  - Bulk modulus
    % mu - Shear modulus
    %
    % Output:
    % C - 6x6 stiffness matrix for the isotropic material

    % Compute I tensor - unity tensor (6x6 identity matrix)
    I = eye(6);

    % Volumetric part of unity tensor
    J = (1/3) * [ones(3, 3), zeros(3, 3); zeros(3, 6)];

    % Deviatoric part of unity tensor
    K = I - J;

    % Define the stiffness tensor using the elastic moduli (k and MU)
    C = (3 * k * J) + (2 * mu * K);
end

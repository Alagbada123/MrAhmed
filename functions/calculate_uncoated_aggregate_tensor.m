function Ainf = calculate_uncoated_aggregate_tensor(Chom_cp, stiffness_uncoated, P_sph, I)
% CALCULATE_UNCOATED_AGGREGATE_TENSOR Calculates strain concentration tensor 
% for uncoated aggregate using Mori-Tanaka scheme
%
% Inputs:
%   Chom_cp          - Homogenized stiffness of cement paste (matrix)
%   stiffness_uncoated - Single structure from stiffness.uncoated array
%                      containing C, k, mu for one uncoated aggregate
%   P_sph            - Hill's polarization tensor for spherical inclusion
%   I                - Fourth-order identity tensor
%
% Output:
%   Ainf             - Strain concentration tensor for uncoated aggregate

    % Calculate strain concentration tensor using Mori-Tanaka scheme
    Ainf = inv(I + P_sph*(stiffness_uncoated.C - Chom_cp));
end
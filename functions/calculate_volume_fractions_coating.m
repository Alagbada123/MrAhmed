function fractions = calculate_volume_fractions_coating(xi, materialProperties, mixDesign, coating_factor, Fpor)
    % CALCULATE_VOLUME_FRACTIONS Calculates the volume fractions for each phase.
    %
    % Syntax:
    %   fractions = calculate_volume_fractions(xi, materialProperties, mixDesign, coating_factor, Fpor)
    %
    % Inputs:
    %   xi                - Hydration degree (nominally between 0 and 1).
    %   materialProperties - Struct containing material properties (e.g., densities).
    %   mixDesign         - Struct containing mix design parameters (e.g., water/cement ratio).
    %   coating_factor    - Factor determining ITZ thickness (default: 0.2).
    %   Fpor              - Porosity increase factor for ITZ (default: 1).
    %
    % Outputs:
    %   fractions         - Struct with nested fields for:
    %       .uncoated_aggregates
    %       .coated_aggregates
    %       .cement_paste
    %       .cement_paste_composition
    %       .hydrate_foam
    %       .ITZ_composition

    %% (0) Default coating factor & Fpor
    if nargin < 4
        coating_factor = 0.2;
    end
    if nargin < 5
        Fpor = 1; % Default porosity increase factor
    end
    
    if xi < 0 || xi > 1
        warning('Hydration degree xi = %.3f is outside [0,1]. Check if that is intended.', xi);
    end

    %% (1) Initialize aggregate fields
    uncoated_names = fieldnames(materialProperties.uncoated_aggregates);
    n_uncoated = length(uncoated_names);
    
    has_coated = isfield(materialProperties, 'coated_aggregates');
    if has_coated
        coated_names = fieldnames(materialProperties.coated_aggregates);
        n_coated = length(coated_names);
    else
        n_coated = 0;
    end
    
    %% (2) Calculate denominator term for volume fractions
    sum_agg_term = 0;
    
    % Uncoated aggregates
    for i = 1:n_uncoated
        agg_name = uncoated_names{i};
        sum_agg_term = sum_agg_term + ...
            mixDesign.aggregateCementRatios.(agg_name) / materialProperties.uncoated_aggregates.(agg_name).rho;
    end
    
    % Coated aggregates, if present
    if has_coated
        for i = 1:n_coated
            agg_name = coated_names{i};
            sum_agg_term = sum_agg_term + ...
                mixDesign.coatedAggregateCementRatios.(agg_name) / materialProperties.coated_aggregates.(agg_name).rho;
        end
    end
    
    denominator = 1 / materialProperties.clinker.rho + ...
                  mixDesign.waterCementRatio / mixDesign.waterDensity + ...
                  sum_agg_term;
    
    %% (3) Compute aggregate fractions
    f_uncoated = zeros(n_uncoated, 1);
    for i = 1:n_uncoated
        agg_name = uncoated_names{i};
        f_uncoated(i) = (mixDesign.aggregateCementRatios.(agg_name) / ...
            materialProperties.uncoated_aggregates.(agg_name).rho) / denominator;
    end
    
    f_coated = zeros(n_coated, 1);
    if has_coated
        for i = 1:n_coated
            agg_name = coated_names{i};
            f_coated(i) = (mixDesign.coatedAggregateCementRatios.(agg_name) / ...
                materialProperties.coated_aggregates.(agg_name).rho) / denominator;
        end
    end
    
    % Total aggregate fraction
    total_agg_fraction = sum(f_uncoated) + sum(f_coated);
    f_cp = 1 - total_agg_fraction;  % Cement paste fraction (before coating)
    
    if f_cp <= 0
        error('Cement paste fraction is negative or zero (%.3f). Check aggregate ratios.', f_cp);
    elseif f_cp < 0.1
        warning('Cement paste fraction is very low (%.3f).', f_cp);
    end
    
    % Allocate fraction for coating if coated aggregates exist
    if has_coated
        f_coating = coating_factor * f_cp;
        f_cp_remaining = (1 - coating_factor) * f_cp;
    else
        f_coating = 0;
        f_cp_remaining = f_cp;
    end
    
    %% (4) Assign fractions to output structure (aggregates)
    % Uncoated aggregates
    for i = 1:n_uncoated
        fractions.uncoated_aggregates.(uncoated_names{i}) = f_uncoated(i);
    end
    
    % Coated aggregates
    if has_coated
        for i = 1:n_coated
            agg_name = coated_names{i};
            fractions.coated_aggregates.(agg_name).core = f_coated(i);
            fractions.coated_aggregates.(agg_name).coating = f_coating / n_coated;
        end
    end
    
    % Cement paste fraction (after coating)
    fractions.cement_paste = f_cp_remaining;
    fractions.ITZ = f_coating;
    
    %% (5) Powers model composition within cement paste
    f_clin_cp = (1 - xi) / ...
        (1 + materialProperties.clinker.rho / mixDesign.waterDensity * mixDesign.waterCementRatio);
    
    f_water_cp = (materialProperties.clinker.rho * (mixDesign.waterCementRatio - 0.42 * xi)) / ...
                 (mixDesign.waterDensity * ...
                 (1 + materialProperties.clinker.rho / mixDesign.waterDensity * mixDesign.waterCementRatio));
    
    f_hyd_cp = (1.42 * materialProperties.clinker.rho * xi) / ...
               (materialProperties.hydrate.rho * ...
               (1 + materialProperties.clinker.rho / mixDesign.waterDensity * mixDesign.waterCementRatio));
    
    f_air_cp = 1 - f_water_cp - f_hyd_cp - f_clin_cp;
    
    fractions.cement_paste_composition = struct(...
        'clinker', f_clin_cp, ...
        'water', f_water_cp, ...
        'hydrates', f_hyd_cp, ...
        'air', f_air_cp ...
    );
    
    %% (6) Hydrate foam fractions within cement paste
    f_hf_total = 1 - f_clin_cp; % Total hydrate foam volume in cement paste
    f_hyd_hf = f_hyd_cp / f_hf_total; % Fraction of foam that is hydrates
    f_por_hf = (f_water_cp + f_air_cp) / f_hf_total; % Fraction of foam that is porosity

    fractions.hydrate_foam = struct(...
        'hydrates', f_hyd_hf, ...
        'porosity', f_por_hf, ...
        'total_in_paste', f_hf_total ...
    );

    %% (7) ITZ Volume Fraction and Hydrate Foam Adjustment
    if has_coated
        % ITZ composition calculations
        f_clin_itz = ((1 - f_por_hf * f_hf_total * Fpor) / ...
                     (1 - f_por_hf * f_hf_total)) * f_clin_cp;

        f_hf_itz = 1 - f_clin_itz; % ITZ total hydrate foam
        f_por_itz = (f_por_hf * f_hf_total * Fpor) / f_hf_itz; % ITZ porosity fraction
        f_hyd_itz = 1 - f_por_itz; % ITZ hydrate fraction

        % Store ITZ-specific values in the output structure
        fractions.ITZ_composition = struct(... 
            'clinker', f_clin_itz, ...
            'total_hydrate_foam', f_hf_itz ...
        );

        % Store ITZ hydrate foam breakdown
        fractions.ITZ_hydrate_foam = struct(... 
            'hydrates', f_hyd_itz, ...
            'porosity', f_por_itz, ...
            'total_hydrate_foam', f_hf_itz ...
        );
    else
        % No ITZ adjustment if no coated aggregates
        fractions.ITZ_composition = struct(... 
            'clinker', NaN, ...
            'total_hydrate_foam', NaN, ...
            'volume_fraction', 0 ...
        );

        fractions.ITZ_hydrate_foam = struct(... 
            'hydrates', NaN, ...
            'porosity', NaN, ...
            'total_hydrate_foam', NaN ...
        );
    end


    %% (8) Verification checks
    % Verify hydrate foam fractions sum to 1
    testSumHF = f_hyd_hf + f_por_hf;
    assert(abs(testSumHF - 1) < 1e-10, ...
        'Hydrate foam fractions do not sum to 1 (sum=%.4g).', testSumHF);
    
    % Verify ITZ fractions sum to 1 if present
    if has_coated
        testSumITZ = f_hyd_itz + f_por_itz;
        assert(abs(testSumITZ - 1) < 1e-10, ...
            'ITZ hydrate foam fractions do not sum to 1 (sum=%.4g).', testSumITZ);
    end
end
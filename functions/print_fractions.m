function print_fractions(fractions)
    fprintf('\n=== Volume Fractions Summary ===\n\n');
    
    % Print uncoated aggregate fractions
    fprintf('Uncoated Aggregates:\n');
    uncoated_names = fieldnames(fractions.uncoated_aggregates);
    for i = 1:length(uncoated_names)
        fprintf('  %s: %.3f\n', uncoated_names{i}, fractions.uncoated_aggregates.(uncoated_names{i}));
    end
    fprintf('\n');
    
    % Print coated aggregate fractions if present
    if isfield(fractions, 'coated_aggregates')
        fprintf('Coated Aggregates:\n');
        coated_names = fieldnames(fractions.coated_aggregates);
        for i = 1:length(coated_names)
            fprintf('  %s:\n', coated_names{i});
            fprintf('    Core: %.3f\n', fractions.coated_aggregates.(coated_names{i}).core);
            fprintf('    Coating: %.3f\n', fractions.coated_aggregates.(coated_names{i}).coating);
        end
        fprintf('\n');
    end
    
    % Print cement paste fraction
    fprintf('Cement Paste: %.3f\n\n', fractions.cement_paste);
    
    % Print cement paste composition
    fprintf('Cement Paste Composition:\n');
    fprintf('  Clinker: %.3f\n', fractions.cement_paste_composition.clinker);
    fprintf('  Water: %.3f\n', fractions.cement_paste_composition.water);
    fprintf('  Hydrates: %.3f\n', fractions.cement_paste_composition.hydrates);
    fprintf('  Air: %.3f\n\n', fractions.cement_paste_composition.air);
    
    % Print hydrate foam composition
    fprintf('Hydrate Foam Composition:\n');
    fprintf('  Hydrates: %.3f\n', fractions.hydrate_foam.hydrates);
    fprintf('  Porosity: %.3f\n', fractions.hydrate_foam.porosity);
    fprintf('  Total in paste: %.3f\n\n', fractions.hydrate_foam.total_in_paste);
    
    % Print total volume check
    total_uncoated = sum(struct2array(fractions.uncoated_aggregates));
    total_coated = 0;
    if isfield(fractions, 'coated_aggregates')
        coated_names = fieldnames(fractions.coated_aggregates);
        for i = 1:length(coated_names)
            total_coated = total_coated + ...
                fractions.coated_aggregates.(coated_names{i}).core + ...
                fractions.coated_aggregates.(coated_names{i}).coating;
        end
    end
    total_volume = total_uncoated + total_coated + fractions.cement_paste;
    
    fprintf('Volume Check:\n');
    fprintf('  Total uncoated aggregates: %.3f\n', total_uncoated);
    fprintf('  Total coated aggregates (with coating): %.3f\n', total_coated);
    fprintf('  Cement paste: %.3f\n', fractions.cement_paste);
    fprintf('  Total volume: %.3f\n', total_volume);
end
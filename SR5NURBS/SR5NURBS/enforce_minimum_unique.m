function population = enforce_minimum_unique(population, source, M, V, minimum_unique)
%ENFORCE_MINIMUM_UNIQUE Replace duplicate copies with the best missing rows.

if nargin < 5 || isempty(minimum_unique) || isempty(population)
    return
end

minimum_unique = max(1, round(minimum_unique));
source_decisions = source(:, 1:V);
available_unique = size(unique(source_decisions, 'rows'), 1);
target_unique = min([minimum_unique, size(population, 1), available_unique]);
if size(unique(population(:, 1:V), 'rows'), 1) >= target_unique
    return
end

rank_column = V + M + 1;
crowding_column = rank_column + 1;
if size(source, 2) >= crowding_column
    ranking_key = [source(:, rank_column), -source(:, crowding_column)];
    [~, source_order] = sortrows(ranking_key, [1 2]);
else
    source_order = (1:size(source, 1))';
end

while size(unique(population(:, 1:V), 'rows'), 1) < target_unique
    selected_decisions = population(:, 1:V);
    replacement_row = find_duplicate_copy(selected_decisions);
    if isempty(replacement_row)
        break
    end

    candidate_row = [];
    for index = reshape(source_order, 1, [])
        if ~ismember(source(index, 1:V), selected_decisions, 'rows')
            candidate_row = index;
            break
        end
    end
    if isempty(candidate_row)
        break
    end
    population(replacement_row, :) = source(candidate_row, :);
end
end

function replacement_row = find_duplicate_copy(decisions)
[~, ~, group] = unique(decisions, 'rows');
counts = accumarray(group, 1);
duplicate_rows = find(counts(group) > 1);
if isempty(duplicate_rows)
    replacement_row = [];
else
    replacement_row = duplicate_rows(end);
end
end

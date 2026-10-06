function [chromosome, iteration_history, completed_generations] = ...
        load_nsga_checkpoint(checkpoint_file, population_size, ...
        target_generations, expected_columns)
%LOAD_NSGA_CHECKPOINT Load the most recent complete NSGA-II generation.

if ~isfile(checkpoint_file)
    error('NSGA2:CheckpointNotFound', ...
        'Checkpoint file not found: %s', checkpoint_file);
end

iteration_history = load(checkpoint_file);
if ~isnumeric(iteration_history) || isempty(iteration_history)
    error('NSGA2:InvalidCheckpoint', ...
        'Checkpoint must contain a nonempty numeric matrix: %s', ...
        checkpoint_file);
end
if size(iteration_history, 2) ~= expected_columns
    error('NSGA2:CheckpointColumnMismatch', ...
        ['Checkpoint has %d columns, but the current configuration ' ...
        'requires %d columns.'], size(iteration_history, 2), ...
        expected_columns);
end
if any(~isfinite(iteration_history(:)))
    error('NSGA2:InvalidCheckpoint', ...
        'Checkpoint contains NaN or Inf values: %s', checkpoint_file);
end

completed_generations = floor(size(iteration_history, 1) / population_size);
if completed_generations < 1
    error('NSGA2:IncompleteCheckpoint', ...
        'Checkpoint does not contain one complete generation.');
end
if completed_generations >= target_generations
    error('NSGA2:CheckpointAlreadyComplete', ...
        ['Checkpoint already contains %d generations; the requested ' ...
        'target is %d.'], completed_generations, target_generations);
end

complete_rows = completed_generations * population_size;
iteration_history = iteration_history(1:complete_rows, :);
chromosome = iteration_history(complete_rows - population_size + 1:complete_rows, :);
end

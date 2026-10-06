function f = genetic_operator(parent_chromosome, M, V, mu, mum, ...
        l_limit, u_limit, mutation_probability, minimum_unique, num_weights)

%% function f  = genetic_operator(parent_chromosome, M, V, mu, mum, l_limit, u_limit)
% 
% This function is utilized to produce offsprings from parent chromosomes.
% The genetic operators corssover and mutation which are carried out with
% slight modifications from the original design. For more information read
% the document enclosed. 
%
% parent_chromosome - the set of selected chromosomes.
% M - number of objective functions
% V - number of decision varaiables
% mu - distribution index for crossover (read the enlcosed pdf file)
% mum - distribution index for mutation (read the enclosed pdf file)
% l_limit - a vector of lower limit for the corresponding decsion variables
% u_limit - a vector of upper limit for the corresponding decsion variables
%
% The genetic operation is performed only on the decision variables, that
% is the first V elements in the chromosome vector. 

%  Copyright (c) 2009, Aravind Seshadri
%  All rights reserved.
%
%  Redistribution and use in source and binary forms, with or without 
%  modification, are permitted provided that the following conditions are 
%  met:
%
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%      
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
%  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
%  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
%  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
%  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%  POSSIBILITY OF SUCH DAMAGE.

[N,m] = size(parent_chromosome);
if nargin < 8 || isempty(mutation_probability)
    mutation_probability = 0.10;
end
if nargin < 9 || isempty(minimum_unique)
    minimum_unique = 1;
end
if nargin < 10 || isempty(num_weights)
    num_weights = V;
end
mutation_probability = min(max(real(mutation_probability), 0), 1);
minimum_unique = max(1, round(minimum_unique));

clear m
p = 1;
% Flags used to set if crossover and mutation were actually performed. 
was_crossover = 0;
was_mutation = 0;


for i = 1 : N
    % Perform crossover unless this reproduction event is selected for
    % mutation. The mutation probability is configured by trajectory_config.
    if rand(1) >= mutation_probability
        % Initialize the children to be null vector.
        child_1 = [];
        child_2 = [];
        % Select the first parent
        parent_1 = round(N*rand(1));
        if parent_1 < 1
            parent_1 = 1;
        end
        % Select the second parent
        parent_2 = round(N*rand(1));
        if parent_2 < 1
            parent_2 = 1;
        end
        % Prefer two different parents, but do not wait forever when the
        % tournament pool has converged to one repeated chromosome.
        parent_attempt = 0;
        max_parent_attempts = max(20, 5 * N);
        while isequal(parent_chromosome(parent_1,:),parent_chromosome(parent_2,:)) ...
                && parent_attempt < max_parent_attempts
            parent_2 = round(N*rand(1));
            if parent_2 < 1
                parent_2 = 1;
            end
            parent_attempt = parent_attempt + 1;
        end
        if isequal(parent_chromosome(parent_1,:),parent_chromosome(parent_2,:))
            for distinct_parent = 1:N
                if ~isequal(parent_chromosome(parent_1,:), ...
                        parent_chromosome(distinct_parent,:))
                    parent_2 = distinct_parent;
                    break
                end
            end
        end
        % Get the chromosome information for each randomnly selected
        % parents
        parent_1 = parent_chromosome(parent_1,:);
        parent_2 = parent_chromosome(parent_2,:);
        % Perform corssover for each decision variable in the chromosome.
        for j = 1 : V
            % SBX (Simulated Binary Crossover).
            % For more information about SBX refer the enclosed pdf file.
            % Generate a random number
            u(j) = rand(1);
            if u(j) <= 0.5
                bq(j) = (2*u(j))^(1/(mu+1));
            else
                bq(j) = (1/(2*(1 - u(j))))^(1/(mu+1));
            end
            % Generate the jth element of first child
            precision = 3;
            if j <= num_weights
                precision = 2;
            end
            child_1(j) = round(0.5*(((1 + bq(j))*parent_1(j)) + (1 - bq(j))*parent_2(j)), precision);
            % each 2 of 50 parents generate 2 children (three variables)
            % Generate the jth element of second child
            child_2(j) = round(0.5*(((1 - bq(j))*parent_1(j)) + (1 + bq(j))*parent_2(j)), precision);
            % Make sure that the generated element is within the specified
            % decision space else set it to the appropriate extrema.
            if child_1(j) > u_limit(j)
                child_1(j) = u_limit(j);
            elseif child_1(j) < l_limit(j)
                child_1(j) = l_limit(j);
            end
            if child_2(j) > u_limit(j)
                child_2(j) = u_limit(j);
            elseif child_2(j) < l_limit(j)
                child_2(j) = l_limit(j);
            end
        end
        % Evaluate the objective function for the offsprings and as before
        % concatenate the offspring chromosome with objective value.
        child_1(:,V + 1: M + V) = evaluate_objective(child_1, M, V);
        child_2(:,V + 1: M + V) = evaluate_objective(child_2, M, V);
        % Set the crossover flag. When crossover is performed two children
        % are generate, while when mutation is performed only only child is
        % generated.
        was_crossover = 1;
        was_mutation = 0;
    % Otherwise perform polynomial mutation.
    else
        % Select at random the parent.
        parent_3 = round(N*rand(1));
        if parent_3 < 1
            parent_3 = 1;
        end
        child_3 = polynomial_mutation(parent_chromosome(parent_3, 1:V), ...
            V, mum, l_limit, u_limit, num_weights);
        % Evaluate the objective function for the offspring and as before
        % concatenate the offspring chromosome with objective value.    
        child_3(:,V + 1: M + V) = evaluate_objective(child_3, M, V);
        % Set the mutation flag
        was_mutation = 1;
        was_crossover = 0;
    end
    % Keep proper count and appropriately fill the child variable with all
    % the generated children for the particular generation.
    if was_crossover
        child(p,:) = child_1;
        child(p+1,:) = child_2;
        was_crossover = 0;
        p = p + 2;
    elseif was_mutation
        child(p,:) = child_3(1,1 : M + V);
        was_mutation = 0;
        p = p + 1;
    end
end
child = ensure_unique_offspring(child, parent_chromosome, M, V, mum, ...
    l_limit, u_limit, minimum_unique, num_weights);
f = child;
end

function child = ensure_unique_offspring(child, parents, M, V, mum, ...
        l_limit, u_limit, minimum_unique, num_weights)
target_unique = minimum_unique;
attempt = 0;
max_attempts = max(100, 20 * target_unique);

while size(unique(child(:, 1:V), 'rows'), 1) < target_unique && ...
        attempt < max_attempts
    parent_index = randi(size(parents, 1));
    candidate = polynomial_mutation(parents(parent_index, 1:V), ...
        V, mum, l_limit, u_limit, num_weights);
    if ~ismember(candidate, child(:, 1:V), 'rows')
        candidate(:, V + 1:V + M) = evaluate_objective(candidate, M, V);
        child(end + 1, :) = candidate; %#ok<AGROW>
    end
    attempt = attempt + 1;
end

% The polynomial mutation normally supplies enough diversity. Use bounded
% random immigrants only as a deterministic fallback for a collapsed pool.
attempt = 0;
while size(unique(child(:, 1:V), 'rows'), 1) < target_unique && ...
        attempt < max_attempts
    candidate = l_limit + (u_limit - l_limit) .* rand(1, V);
    candidate = round_and_bound(candidate, V, l_limit, u_limit, num_weights);
    if ~ismember(candidate, child(:, 1:V), 'rows')
        candidate(:, V + 1:V + M) = evaluate_objective(candidate, M, V);
        child(end + 1, :) = candidate; %#ok<AGROW>
    end
    attempt = attempt + 1;
end

if size(unique(child(:, 1:V), 'rows'), 1) < target_unique
    error('genetic_operator:InsufficientDiversity', ...
        'Unable to generate %d unique offspring chromosomes.', target_unique);
end
end

function child = polynomial_mutation(parent, V, mum, l_limit, u_limit, num_weights)
child = parent(1:V);
for j = 1:V
    random_value = rand(1);
    if random_value < 0.5
        delta = (2 * random_value)^(1 / (mum + 1)) - 1;
    else
        delta = 1 - (2 * (1 - random_value))^(1 / (mum + 1));
    end
    child(j) = child(j) + delta;
end
child = round_and_bound(child, V, l_limit, u_limit, num_weights);
end

function decision = round_and_bound(decision, V, l_limit, u_limit, num_weights)
for j = 1:V
    precision = 3;
    if j <= num_weights
        precision = 2;
    end
    decision(j) = round(decision(j), precision);
    decision(j) = min(max(decision(j), l_limit(j)), u_limit(j));
end
end

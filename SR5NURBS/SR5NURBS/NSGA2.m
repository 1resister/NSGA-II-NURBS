function NSGA2(pop, gen, M, V, min_range, max_range, resume_from_checkpoint)

%% function nsga_2(pop,gen)
% is a multi-objective optimization function where the input arguments are 
% pop - Population size
% gen - Total number of generations
% 
% This functions is based on evolutionary algorithm for finding the optimal
% solution for multiple objective i.e. pareto front for the objectives. 
% Initially enter only the population size and the stoping criteria or
% the total number of generations after which the algorithm will
% automatically stopped. 
%
% You will be asked to enter the number of objective functions, the number
% of decision variables and the range space for the decision variables.
% Also you will have to define your own objective funciton by editing the
% evaluate_objective() function. A sample objective function is described
% in evaluate_objective.m. Kindly make sure that the objective function
% which you define match the number of objectives that you have entered as
% well as the number of decision variables that you have entered. The
% decision variable space is continuous for this function, but the
% objective space may or may not be continuous.
%
% Original algorithm NSGA-II was developed by researchers in Kanpur Genetic
% Algorithm Labarotary and kindly visit their website for more information
% http://www.iitk.ac.in/kangal/


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

%% Simple error checking
% Number of Arguments
% Check for the number of arguments. The two input arguments are necessary
% to run this function.
if nargin < 2
    error('NSGA-II: Please enter the population size and number of generations as input arguments.');
end
if nargin < 7
    resume_from_checkpoint = false;
end
% Both the input arguments need to of integer data type
if isnumeric(pop) == 0 || isnumeric(gen) == 0
    error('Both input arguments pop and gen should be integer datatype');
end
% Minimum population size has to be 20 individuals
if pop < 20
    error('Minimum population for running this function is 20');
end
if gen < 5
    error('Minimum number of generations is 5');
end
% Make sure pop and gen are integers
pop = round(pop);
gen = round(gen);
cfg = trajectory_config(fileparts(mfilename('fullpath')));

%% Initialize the population
% Population is initialized with random values which are within the
% specified range. Each chromosome consists of the decision variables. Also
% the value of the objective functions, rank and crowding distance
% information is also added to the chromosome vector but only the elements
% of the vector which has the decision variables are operated upon to
% perform the genetic operations like corssover and mutation.

% the chromosome is composed of the variables and function values
% chromosome = [x(1) x(2) f_1 f_2];
if resume_from_checkpoint
    [chromosome, Iteration, completed_generations] = load_nsga_checkpoint( ...
        cfg.files.chromosome_iter, pop, gen, V + M + 2);
    start_generation = completed_generations + 1;
    fprintf('Resuming NSGA-II from generation %d; target generation %d.\n', ...
        completed_generations, gen);
else
    chromosome = initialize_variables(pop, M, V, min_range, max_range, ...
        cfg.nsga.num_weights, cfg.nsga.reference_decision);

    %% Sort the initialized population
    % Sort the population using non-domination-sort. This returns two columns
    % for each individual which are the rank and the crowding distance
    % corresponding to their position in the front they belong. At this stage
    % the rank and the crowding distance for each chromosome is added to the
    % chromosome vector for easy of computation.
    chromosome = non_domination_sort_mod(chromosome, M, V);
    Iteration = zeros(0, size(chromosome, 2));
    start_generation = 1;

    %% save the initial chromosome
    disp('==================================================================');
    disp('Initial ranked chromosomes generated ...');
    disp('==================================================================');
    save(cfg.files.chromosome_initial, 'chromosome', '-ASCII');
end

%% Start the evolution process
% The following are performed in each generation
% * Select the parents which are fit for reproduction
% * Perfrom crossover and Mutation operator on the selected parents
% * Perform Selection from the parents and the offsprings
% * Replace the unfit individuals with the fit individuals to maintain a
%   constant population size.

for i = start_generation : gen
    generation_timer = tic;
    % Select the parents
    % Parents are selected for reproduction to generate offspring. The
    % original NSGA-II uses a binary tournament selection based on the
    % crowded-comparision operator. The arguments are 
    % pool - size of the mating pool. It is common to have this to be half the
    %        population size.
    % tour - Tournament size. Original NSGA-II uses a binary tournament
    %        selection, but to see the effect of tournament size this is kept
    %        arbitary, to be choosen by the user.
     pool = round(pop/2);
     tour = 2;
    % Selection process
    % A binary tournament selection is employed in NSGA-II. In a binary
    % tournament selection process two individuals are selected at random
    % and their fitness is compared. The individual with better fitness is
    % selcted as a parent. Tournament selection is carried out until the
    % pool size is filled. Basically a pool size is the number of parents
    % to be selected. The input arguments to the function
    % tournament_selection are chromosome, pool, tour. The function uses
    % only the information from last two elements in the chromosome vector.
    % The last element has the crowding distance information while the
    % penultimate element has the rank information. Selection is based on
    % rank and if individuals with same rank are encountered, crowding
    % distance is compared. A lower rank and higher crowding distance is
    % the selection criteria.
    parent_chromosome = tournament_selection(chromosome, pool, tour);
    parent_chromosome = enforce_minimum_unique(parent_chromosome, ...
        chromosome, M, V, cfg.nsga.min_unique_chromosomes);
    parent_unique = size(unique(parent_chromosome(:, 1:V), 'rows'), 1);
    required_parent_unique = min(cfg.nsga.min_unique_chromosomes, ...
        size(parent_chromosome, 1));
    if parent_unique < required_parent_unique
        error('NSGA2:InsufficientParentDiversity', ...
            'Parent pool contains %d unique chromosomes; %d are required.', ...
            parent_unique, required_parent_unique);
    end

    % Perfrom crossover and Mutation operator
    % The original NSGA-II algorithm uses Simulated Binary Crossover (SBX) and
    % Polynomial  mutation. Crossover probability pc = 0.9 and mutation
    % probability is pm = 1/n, where n is the number of decision variables.
    % Both real-coded GA and binary-coded GA are implemented in the original
    % algorithm, while in this program only the real-coded GA is considered.
    % The distribution indeices for crossover and mutation operators as mu = 20
    % and mum = 20 respectively.
    mu = 20;
    mum = 20;
    offspring_chromosome = genetic_operator(parent_chromosome, ...
        M, V, mu, mum, min_range, max_range, ...
        cfg.nsga.mutation_probability, cfg.nsga.min_unique_chromosomes, ...
        cfg.nsga.num_weights);

    % Intermediate population
    % Intermediate population is the combined population of parents and
    % offsprings of the current generation. The population size is two
    % times the initial population.
    
    [main_pop,~] = size(chromosome);
    [offspring_pop,~] = size(offspring_chromosome);
    % intermediate_chromosome is a concatenation of current population and
    % the offspring population.
    intermediate_chromosome(1:main_pop,:) = chromosome; % original parents (50*6)
    intermediate_chromosome(main_pop + 1 : main_pop + offspring_pop,1 : M+V) = ...
        offspring_chromosome; % childen (50*6)

    % Non-domination-sort of intermediate population
    % The intermediate population is sorted again based on non-domination sort
    % before the replacement operator is performed on the intermediate
    % population.
    intermediate_chromosome = ...
        non_domination_sort_mod(intermediate_chromosome, M, V);
    % Perform Selection
    % Once the intermediate population is sorted only the best solution is
    % selected based on it rank and crowding distance. Each front is filled in
    % ascending order until the addition of population size is reached. The
    % last front is included in the population based on the individuals with
    % least crowding distance
    chromosome = replace_chromosome(intermediate_chromosome, M, V, pop, ...
        cfg.nsga.min_unique_chromosomes);
    generation_seconds = toc(generation_timer);
    unique_chromosomes = size(unique(chromosome(:, 1:V), 'rows'), 1);
    required_unique = min(cfg.nsga.min_unique_chromosomes, pop);
    if unique_chromosomes < required_unique
        error('NSGA2:InsufficientPopulationDiversity', ...
            'Generation %d contains %d unique chromosomes; %d are required.', ...
            i, unique_chromosomes, required_unique);
    end
    fprintf(['Generation %d/%d completed: %.1f s, %d offspring, ' ...
        '%d unique chromosomes.\n'], i, gen, generation_seconds, ...
        offspring_pop, unique_chromosomes);
    drawnow;
    
    % save and plot the iterations
    Iteration((i-1)*pop+1:i*pop,:) = chromosome;
    if cfg.nsga.save_interval > 0 && (mod(i, cfg.nsga.save_interval) == 0 || i == gen)
        save(cfg.files.chromosome_iter, 'Iteration', '-ASCII');
    end

    if cfg.nsga.plot_interval > 0 && ~mod(i, cfg.nsga.plot_interval)
        disp('==========================================================');
        fprintf('%d generations completed\n',i);
        disp('==========================================================');
        if M == 2
            plot(chromosome(:,V + 1),chromosome(:,V + 2),'o');
            grid on
            title('NSGA-II internal penalized objectives')
        elseif M == 3
            plot3(chromosome(:,V + 1),chromosome(:,V + 2),chromosome(:,V + 3),'o');
            grid on
            title('NSGA-II internal penalized objectives')
        end
        pause(0.5)
    end
end

chromosome = abs(chromosome); % 100*8 matrix
%% Result
% Save the result in ASCII text format.
save(cfg.files.chromosome_final, 'chromosome', '-ASCII', '-double');
% Build the raw-objective cache from exactly the serialized population that
% main_function/export_selected_solution will load. Otherwise ASCII rounding
% changes the exact cache signature and triggers an unnecessary second pass.
chromosome = load(cfg.files.chromosome_final);

%% Visualize the final front using raw, unpenalized physical objectives.
raw_candidates = list_pareto_candidates(chromosome, cfg);
if M == 2 && ~isempty(raw_candidates)
    plot(raw_candidates(:,2), raw_candidates(:,3), 'o');
    grid on
    xlabel('Time'); ylabel('Error');
    title('Pareto front: raw objectives without penalties')
elseif M == 3 && ~isempty(raw_candidates)
    plot3(raw_candidates(:,2), raw_candidates(:,3), raw_candidates(:,4), 'o');
    grid on
    xlabel('Time'); ylabel('Error'); zlabel('Vibration');
    title('Pareto front: raw objectives without penalties')
end

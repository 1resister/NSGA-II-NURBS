function f = initialize_variables( ...
        pop, M, V, min_range, max_range, num_weights, reference_decision)

%% function f = initialize_variables(N, M, V, min_tange, max_range) 
% This function initializes the chromosomes. Each chromosome has the
% following at this stage
%       * set of decision variables
%       * objective function values
% 
% where,
% N - Population size
% M - Number of objective functions
% V - Number of decision variables
% min_range - A vector of decimal values which indicate the minimum value
% for each decision variable.
% max_range - Vector of maximum possible values for decision variables.

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

lower_bound = min_range;
upper_bound = max_range;
if nargin < 6 || isempty(num_weights)
    num_weights = V;
end
if nargin < 7 || isempty(reference_decision)
    reference_decision = ones(1, V);
end
if numel(reference_decision) ~= V || any(~isfinite(reference_decision))
    error('reference_decision must contain one finite value per variable.');
end

% K is the total number of array elements. For ease of computation decision
% variables and objective functions are concatenated to form a single
% array. For crossover and mutation only the decision variables are used
% while for selection, only the objective variable are utilized.

K = M + V;
f = zeros(pop, K);

%% Initialize each chromosome
% For each chromosome perform the following (N is the population size)
for i = 1 : pop
    % Initialize the decision variables based on the configured bounds.
    % The first group is NURBS weights. The remaining variables contain
    % time/velocity-node scale factors followed by endpoint-jerk utilization.
    if i == 1
        % Keep one unit-weight/unit-speed reference. The dynamically fitted
        % baseline with zero endpoint jerk is represented in every fresh
        % optimization.
        f(i, 1:V) = min(max(reference_decision, lower_bound), upper_bound);
    else
        for j = 1 : V
            precision = 3;
            if j <= num_weights
                precision = 2;
            end
            f(i, j) = round(lower_bound(j) + ...
                (upper_bound(j) - lower_bound(j))*rand(1), precision);
        end
    end
    % For ease of computation and handling data the chromosome also has the
    % vlaue of the objective function concatenated at the end. The elements
    % V + 1 to K has the objective function valued. 
    % The function evaluate_objective takes one chromosome at a time,
    % infact only the decision variables are passed to the function along
    % with information about the number of objective functions which are
    % processed and returns the value for the objective functions. These
    % values are now stored at the end of the chromosome itself.
    f(i, V + 1 : K) = evaluate_objective(f(i, :), M, V);
end

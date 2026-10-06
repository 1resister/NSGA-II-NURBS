function [S, V, A, J, snap] = evaluate_septic_transition( ...
        time, s0, v0, v1, T, j0, j1)
% Evaluate a seventh-order velocity transition.
%
% The legacy five-argument call keeps zero acceleration, jerk and snap at
% both ends. Optional j0/j1 replace only the endpoint jerk conditions;
% acceleration and snap remain zero. Adjacent SPUs can therefore share one
% physical jerk value at an internal LMSC without changing the zero-jerk
% baseline used by NSGA-II.

if nargin < 6 || isempty(j0)
    j0 = 0;
end
if nargin < 7 || isempty(j1)
    j1 = 0;
end
if ~isscalar(j0) || ~isscalar(j1) || ...
        ~isfinite(j0) || ~isfinite(j1)
    error('Endpoint jerk values must be finite scalars.');
end

if T <= 0
    S = s0 + v1 * time;
    V = v1 * ones(size(time));
    A = zeros(size(time));
    J = zeros(size(time));
    snap = zeros(size(time));
    return
end

tau = min(max(time / T, 0), 1);
delta_v = v1 - v0;
jerk_tolerance = 10 * eps(max([1, abs(j0), abs(j1)]));
if abs(j0) <= jerk_tolerance && abs(j1) <= jerk_tolerance
    % Preserve the original expression exactly for regression compatibility.
    g = 35*tau.^4 - 84*tau.^5 + 70*tau.^6 - 20*tau.^7;
    integral_g = 7*tau.^5 - 14*tau.^6 + 10*tau.^7 - 2.5*tau.^8;
    V = v0 + delta_v .* g;
    A = delta_v / T .* (140*tau.^3 .* (1 - tau).^3);
    J = delta_v / T^2 .* ...
        (420*tau.^2 .* (1 - tau).^2 .* (1 - 2*tau));
    S = s0 + v0 .* time + delta_v * T .* integral_g;
    if nargout >= 5
        snap = delta_v / T^3 .* (840*tau - 5040*tau.^2 + ...
            8400*tau.^3 - 4200*tau.^4);
    end
    return
end

% V(tau) = sum(k=0:7) b(k+1)*tau^k.  In normalized time, the
% acceleration/jerk/snap boundary values are multiplied by T, T^2 and T^3.
% b0..b3 are fixed by the four start conditions.  The remaining four
% coefficients are obtained from the end conditions.
b = zeros(8, 1);
b(1) = v0;
b(3) = 0.5 * T^2 * j0;
persistent endpoint_inverse
if isempty(endpoint_inverse)
    endpoint_matrix = [ ...
        1,   1,   1,   1; ...
        4,   5,   6,   7; ...
        12, 20,  30,  42; ...
        24, 60, 120, 210];
    endpoint_inverse = inv(endpoint_matrix);
end
rhs = [ ...
    v1 - b(1) - b(3); ...
    -2 * b(3); ...
    T^2 * j1 - 2 * b(3); ...
    0];
b(5:8) = endpoint_inverse * rhs;

V = zeros(size(tau));
A = zeros(size(tau));
J = zeros(size(tau));
snap = zeros(size(tau));
S = s0 * ones(size(tau));
for k = 0:7
    V = V + b(k + 1) .* tau.^k;
    S = S + T * b(k + 1) / (k + 1) .* tau.^(k + 1);
    if k >= 1
        A = A + k * b(k + 1) / T .* tau.^(k - 1);
    end
    if k >= 2
        J = J + k * (k - 1) * b(k + 1) / T^2 .* tau.^(k - 2);
    end
    if k >= 3
        snap = snap + k * (k - 1) * (k - 2) * b(k + 1) / T^3 .* ...
            tau.^(k - 3);
    end
end
end

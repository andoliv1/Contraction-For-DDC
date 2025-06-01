function u = compute_feedback_using_geodesic(X_opt, N, rho_1, Q, G,vars)
% computeController - Numerically integrates
% u = ∫₀¹ -1/2 * rho(γ(s)) * B' * W(γ(s)) * ∂γ/∂s ds
%
% Inputs:
%   X_opt  - n x (N+1) matrix of geodesic points γ(s_i)
%   N      - number of intervals in [0,1]
%   rho_fun - function handle @(x) -> scalar rho(x)
%   W_fun  - function handle @(x) -> n x n matrix W(x)
%   B      - n x m control input matrix
%
% Output:
%   u      - m x 1 control vector

    ds = 1 / N;
    [n, ~] = size(X_opt);
    [~, m] = size(G);
    u_sum = zeros(m,1);

    for i = 1:N
        x_curr = X_opt(:, i);
        x_next = X_opt(:, i+1);

        dx_ds = (x_next - x_curr) / ds;          % Approximate ∂γ/∂s at s_i        
        Q_curr_inv = inv(replace(Q,vars,x_curr));     % Evaluate W at γ(s_i)
        rho_curr = replace(rho_1,vars,x_curr);    % Evaluate ρ at γ(s_i)

        integrand = -0.5 * rho_curr * (G') * Q_curr_inv * dx_ds;  % m x 1 vector

        u_sum = u_sum + integrand * ds;
    end

    u = u_sum;
end

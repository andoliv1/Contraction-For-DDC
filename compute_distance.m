function [d,energy] = compute_distance(Q,x,p,j,W_Val,N)
[X_opt,f_opt] = compute_geodesic_casadi(p,[0;0;0],W_Val,N);
ds = 1 / N;
d = 0;
energy_int = 0;
for i = 1:N
    x_curr = X_opt(:, i);
    x_next = X_opt(:, i+1);

    dx_ds = (x_next - x_curr) / ds;          % Approximate ∂γ/∂s at s_i        
    Q_curr_inv = inv(replace(Q,x,x_curr));     % Evaluate W at γ(s_i)

    integrand = sqrt(dx_ds' * Q_curr_inv * dx_ds);  % m x 1 vector
    energy_int = energy_int + (dx_ds' * Q_curr_inv * dx_ds - f_opt)^2;
    d = d + integrand * ds;
end

energy = 1/(f_opt)*sqrt(energy_int);

end
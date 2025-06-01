function [X_opt,f_opt] = compute_geodesic_casadi(x0, xT, W_val, N)
    % Compute geodesic under nonlinear Riemannian metric G(x)
    % Inputs:
    %   x0     : 3x1 initial point
    %   xT     : 3x1 final point
    %   W_val  : 12x12 weight matrix (positive definite)
    %   N      : number of discretization intervals
    % W_val = zeros(12,12);
    % W_val(1:3, 1:3) = eye(3);
    
    import casadi.*
    % tic;
    if nargin < 4
        N = 40;  % Default discretization
    end

    n = 3;       % Dimension of state space
    % d = 4;       % Feature length: [1; x1; x2; x3]
    T = 1.0;     % Final time
    dt = T / N;

    % Define symbolic variables
    x = SX.sym('x', n, 1);             % x = [x1; x2; x3]
    v = [1; x(1); x(2); x(3)];         % Feature map
    % Precompute fixed W_val metric
    phi_sym = kron(v, eye(n));               % Symbolic φ(x)
    G_expr = inv(phi_sym' * W_val * phi_sym);     % Now W_val is directly used
    G_fun = Function('G_fun', {x}, {G_expr}); % Only takes x now

    % Optimization variables
    X = SX.sym('X', n, N+1);          % States
    w = X(:);
    J = 0;

    % Build cost: geodesic length
    for k = 1:N
        xk = X(:,k);
        uk = (X(:,k+1) - X(:,k)) / dt;
        Gk = G_fun(xk);
        J = J + (uk' * Gk * uk) * dt;
    end

    % Boundary constraints
    g = [X(:,1) - x0; X(:,end) - xT];
    % Set up NLP
    nlp = struct('x', w, 'f', J, 'g', g);
    
    % Solver options to suppress logs
    opts = struct;
    opts.print_time = false;                        % Disable timing info
    opts.ipopt.print_level = 0;                     % IPOPT log level 0 = no output
    opts.ipopt.sb = 'yes';                          % Suppress IPOPT banner

% Create solver
    solver = nlpsol('solver', 'ipopt', nlp, opts);
    % Initial guess: straight line
    w0 = zeros(n, N+1);
    for i = 1:N+1
        w0(:,i) = x0 + (i-1)/N*(xT - x0);
    end

    % Solve
    sol = solver('x0', w0(:), 'lbg', zeros(size(g)), 'ubg', zeros(size(g)));
    w_opt = full(sol.x);
    X_opt = reshape(w_opt, n, N+1);
    f_opt = full(sol.f);

    % elapsed = toc;
    % fprintf('it took this amount of time %.6f s\n', elapsed);


    % Optional: 3D plot. Uncomment this to see what approximate geodesics
    % look like
    % figure;
    % plot3(X_opt(1,:), X_opt(2,:), X_opt(3,:), 'b-', 'LineWidth', 2); hold on;
    % plot3(x0(1), x0(2), x0(3), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    % plot3(xT(1), xT(2), xT(3), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
    % xlabel('x_1'); ylabel('x_2'); zlabel('x_3');
    % title('Geodesic under nonlinear 3D metric');
    % axis equal; grid on;
end

%% find contraction metric in model based way
% yalmip('clear')
% K = sdpvar(2,4);
% partial_u = -(1/2)*G'*inv(Q)*rho_1;
% Phi = [x(1);x(2);x(1)^2;x(1)^3];
% c = coefficients(partial_u - K*jacobian(Phi,x),x);
options = sdpsettings('solver','mosek','verbose',1,'MSK_DPAR_DATA_TOL_X',1e-4, 'MSK_DPAR_INTPNT_CO_TOL_PFEAS', 1e-4,'MSK_DPAR_DATA_TOL_QIJ', 1e-4); % Below is nonlinear metric
% sol = optimize(c == 0,[],options,K(:));
% value(K)


F = sdpvar(2,4);
optimize(F*(phi_data') - (x_dot') == 0,[],options,F(:)) 



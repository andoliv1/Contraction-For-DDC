function [L] = check_results(vars,power,Q)
    F = [vars(2) - (3/2)*vars(1)^2 - (1/2)*vars(1)^3; 3*vars(1) - vars(2)];
    v = monolist(vars',power,0);
    L = sdpvar(length(v)*2,length(v)*2);
    P = kron(v,eye(2))'*L*kron(v,eye(2));
    Contraction = Q*(jacobian(F,vars)') + (jacobian(F,vars)')*Q + Q + jacobian(Q,vars(1))*F(1) + jacobian(Q,vars(2))*F(2);
    % coefficients(P,vars)
    sdisplay(coefficients(P,vars))
    sdisplay(coefficients(Contraction,vars))
    % coefficients(Contraction,vars)
    constraints = [L <= 1e-2 ,...
        sos(-Contraction,vars) == 0];
    options = sdpsettings('solver','mosek','verbose',1);
    optimize(constraints, [], options, L(:))
    % coefficients(P,vars);
end
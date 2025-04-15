function [] = verify_contraction_metric(x,y,rho_1,Q,G,autonomous)
    % Define the differential equations
    global name_of_system
    dxdt = select_dynamical_system(name_of_system,x);
    if autonomous == 1
        Z = jacobian(dxdt,x)*Q + Q*(jacobian(dxdt,x))' - jacobian(Q,x(1))*dxdt(1) ...
            - jacobian(Q,x(2))*dxdt(2) + 0.05*Q;
    else
        Z = jacobian(dxdt,x)*Q + Q*(jacobian(dxdt,x))' - jacobian(Q,x(1))*dxdt(1) ...
            - jacobian(Q,x(2))*dxdt(2) + 0.05*Q - rho_1*G*G';
    end
    % + y'*y*sigma*(4 - x'*x)
    options = sdpsettings('solver','mosek','verbose',1); % Below is nonlinear metric
    F = sos(-y'*Z*y);
    solvesos(F,[],options);
end
function [] = verify_contraction_metric(x,y,rho_1,Q,G,sigma,autonomous)
    % Define the differential equations
    % dxdt = [-x(2) - 3/2*(x(1)^2) - 1/2*(x(1)^3);3*x(1) - x(2)];
    % dxdt = [0.4285   -0.4298; 0.4018 1.3036]*[x(1);x(2)];
    % dxdt = [4 - 0.5*x(1) + x(2) - (1-x(2))*x(1); -x(1) + (1-x(2))*x(1)];
    % dxdt = [-x(1) + x(1)*x(2); -x(2)];
    dxdt = [-x(1) + x(3); x(1)^2 - x(2) - 2*x(1)*x(3) + x(3); -x(2)];
    if autonomous == 1
        Z = jacobian(dxdt,x)*Q + Q*(jacobian(dxdt,x))' - jacobian(Q,x(1))*dxdt(1) ...
            - jacobian(Q,x(2))*dxdt(2) - jacobian(Q,x(3))*dxdt(3) + 0.05*Q;
    else
        Z = jacobian(dxdt,x)*Q + Q*(jacobian(dxdt,x))' - jacobian(Q,x(1))*dxdt(1) ...
            - jacobian(Q,x(2))*dxdt(2) - jacobian(Q,x(3))*dxdt(3) + 0.05*Q - rho_1*G*G';
    end
    % + y'*y*sigma*(4 - x'*x)
    options = sdpsettings('solver','mosek','verbose',1); % Below is nonlinear metric
    sdisplay(Z)
    F = sos(-y'*Z*y);
    solvesos(F,[],options);
end
function dxdt = sys_poly(t,x,eps,feedback,data_points,dimensions,vars,rho_1,G,Q,W)
global x_glo
global index
global update
global t_prev
global interval
global trajectories
if feedback == 1
    t
    N = 100;
    [X_opt,~] = compute_geodesic_casadi([0;0;0],x,W,N);
    u = compute_feedback_using_geodesic(X_opt,N,rho_1,Q,G,vars);
    dxdt = [-x(1) + x(3); x(1)^2 - x(2) - 2*x(1)*x(3) + x(3); -x(2)] +  G*u;
else
    dxdt = [-x(1) + x(3); x(1)^2 - x(2) - 2*x(1)*x(3) + x(3); -x(2)];
    if(((mod(index,data_points/trajectories) ~= 1) || update) && t - t_prev > interval)
        x_glo(index,1) = t;
        x_glo(index,2:2+dimensions -1) = x;
        x_glo(index, 2+dimensions:end) = dxdt;
        index = index + 1;
        t_prev = t;
        update = false;
    end
end
end
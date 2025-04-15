function dxdt = sys_poly(t,x,eps,feedback,data_points,dimensions,vars,rho_1,G,Q)
global x_glo
global index
global update
global t_prev
global interval
global trajectories
if feedback == 1
    t
    u = compute_feedback(rho_1,Q,vars,[x(1);x(2);x(3)],G);
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
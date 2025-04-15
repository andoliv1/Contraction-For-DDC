function dxdt = sys_poly(t,x,eps,feedback,data_points,vars,rho_1,G,Q)
global x_glo
global index
global update
global name_of_system
global t_prev
global interval
global trajectories
dxdt = select_dynamical_system(name_of_system,x);
t
if feedback == 1
    u = compute_feedback(rho_1,inv(Q),vars,[x(1);x(2)],G);
    dxdt = dxdt +  G*u;
else
    if(((mod(index,data_points/trajectories) ~= 1) || update) && t - t_prev > interval)
        x_glo(index,1) = t;
        x_glo(index,2) = x(1);
        x_glo(index,3) = x(2);
        x_glo(index, 4) = dxdt(1);
        x_glo(index, 5) = dxdt(2);
        index = index + 1;
        t_prev = t;
        update = false;
    end
end
end
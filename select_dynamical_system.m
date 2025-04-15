function [ode] = select_dynamical_system(name,x)
    if (strcmp(name,'unstable-linear'))
        ode = [0.4285   -0.4298; 0.4018 1.3036]*[x(1);x(2)];
    elseif strcmp(name,'unstable-nonlinear')
        ode = [-x(2) - 3/2*(x(1)^2) - 1/2*(x(1)^3);3*x(1) + x(2)];
    elseif strcmp(name,'sontag')
        ode = [4 - 0.5*x(1) + x(2) - (1-x(2))*x(1); -x(1) + (1-x(2))*x(1)];
    else
        ode = [-x(1) + x(1)*x(2); -x(2)];
    end
end
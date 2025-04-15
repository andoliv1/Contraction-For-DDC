%% Steps that I need to do to implement my problem
% 1 - Create data from an ODE
% 2 - Form the feasibility matrices that is 2\partial \phiWyy^T - \Sigma
% u_i,k\Phi_{i,k}^T = 0...
% 3 - Create the inequality in (x,y) and solve it as a system of polynomial
% inequalities
% The above problem cannot be formulated as a 
clearvars -global
clear all; close all
clc
rng('default')
%% Step 1 - generate data from an ODE which is defined in sys_poly function
global x_glo
global t_prev
global interval
global update
global index
global trajectories
data_points = 60;
dimensions = 3;
trajectories = 6;
x_glo = zeros(data_points,2*dimensions+1);
interval = 0.1;
index = 1;
update = true;
G = [0;0;1]; 
tspan = 0:0.01:1;          % excitation time
for i=1:trajectories
    initial = rand(dimensions,1).*cos(pi*rand(dimensions,1));    % initial state   for some systems, where to start matters
    t_prev = 0;
    update = true;
    [t, x] = ode15s(@(t,x) sys_poly(t,x,0,0,data_points,dimensions), tspan, initial);
end
update = false;
%% Corrupt the state derivative with noise
eps_1 = max(max(abs(x_glo(1:end,2+dimensions))));
eps_2 = max(max(abs(x_glo(1:end,2+dimensions+1))));
eps_3 = max(max(abs(x_glo(1:end,2+dimensions+2))));
eps = max([eps_1;eps_2;eps_3])/15;
for i = 1:data_points
    x_glo(i,2+dimensions:end) = x_glo(i,2+dimensions:end) + eps*(2*rand(1,dimensions)-1);
end
%% build A,B,xi
T = data_points;      % of sample for design, increase if noise large
[phi_data,x_dot,data_fit] = data_matrix(T,eps,dimensions);
%% clear yalmip
yalmip('clear')
%% find contraction metric in data driven way
tic
degree = 1;    %% this is the variable used to indicate the maximum degree of the feedback controller
degree_rho = 3;
[x,y,mu,rho_1,Q,W,c,rho_Y] = ddc_poly(phi_data,T,data_fit,degree,degree_rho,G);
toc
%% verify contraction metric in oringal model
verify_contraction_metric(x,y,rho_2,M,G,sigma,0)
%% Create trajectories under the contraction feedback controller found
Z = [];
colors = lines(10);
for i = 1:5
    tspan = 0:0.1:5;
    % initial = 10*rand(3,1);
    % initial = [cos(pi*rand(1,1))*10*rand(1,1);cos(pi*rand(1,1))*10*rand(1,1);cos(pi*rand(1,1))*10*rand(1,1)];
    initial = [-9 + 4*(i-1);-9 + 4*(i-1);-9 + 4*(i-1)];
    [t,Z(i,:,:,:)] = ode15s(@(t,z) sys_poly(t,z,eps,1,data_points,dimensions,x,rho_1,G,Q), tspan, initial);
    'here'
end
%% evaluate the feedback control action to show all trajectories are contracting
set(groot, ['Default', 'Line', 'LineWidth'], 2);
set(groot, 'DefaultAxesFontSize', 14);
set(groot, 'DefaultAxesFontWeight', 'bold');
figure();
for i = 1:5
    plot(t,Z(i,:,1),'-*', 'Color',colors(i,:));
    if i == 1
        xlabel('T');
        ylabel('X_1');
        xlim([0,5])
    end
    hold on;
end
% saveas(gcf,'XCoordinatePlot-NonlinearNonAutonomous-3D-W-Rho','jpg')
hold off
figure();
for i = 1:5
    plot(t,Z(i,:,2),'-*', 'Color',colors(i,:));
    if i == 1
        xlabel('T');
        ylabel('X_2');
        xlim([0,5])
    end
    hold on;
end
% saveas(gcf,'YCoordinatePlot-NonlinearNonAutonomous-3D-W-Rho','jpg')
hold off
figure();
for i = 1:5
    plot(t,Z(i,:,3),'-*', 'Color',colors(i,:));
    if i == 1
        xlabel('T');
        ylabel('X_3');
        xlim([0,5])
    end
    hold on;
end
% saveas(gcf,'ZCoordinatePlot-NonLinearNonAutonomous-3D-W-Rho','jpg')
hold off;
figure();
for i = 1:5
    plot3(Z(i,:,1),Z(i,:,2),Z(i,:,3),'-*', 'Color',colors(i,:));
    if i == 1
        xlabel('X_1');
        ylabel('X_2');
        zlabel('X_3');
        grid on
    end
    hold on;
end
% saveas(gcf,'3DTrajectoriesPlot-NonlinearNonAutonomous-W-Rho','jpg')
hold off;
%% Compute distance of trajectories
 l = [];
 for i = 1:5
     for j=1:51
         l(i,j) = compute_distance(Q,x,squeeze(Z(i,j,:)),j);
     end
     'here'
 end
 %% Plot distance
 t_span = 0:0.1:5;
 for i = 1:5
     semilogy(t_span,l(i,:),'-*', 'Color',colors(i,:));
     if i == 1
         xlabel('Time');
         ylabel('~log(d(x(t),0)');
         xlim([0,5]);
     end
     hold on;
 end
 hold off;
 saveas(gcf,'ApproxDistEven-W-Rho','jpg')
%% Plot Many Different Solutions to our ODE
Z = [];
colors = lines(10);
figure;
hold on;
for i=1:10
    i
    tspan = 0:1:10;
    initial = 4*rand(2,1);
    eps = 0;             % no noise
    [t,Z(i,:,:)] = ode15s(@(t,x) sys_poly(t,x,eps,0,100), tspan, initial);
    Z(i,:,:)
    plot(Z(i,:,1), Z(i,:,2), '-o', 'Color',colors(i,:));
    % plotpp(odefun,'tspan',10,'quivercolor', [0.6,0.6,0.6],'xlim',[-50,50],'ylim',[-50,50]);
end
xlabel('X_1');
ylabel('X_2');
saveas(gcf,'TrajectoriesPlot-NonLinearNonAutonomous-Unstable','jpg')
grid on;
hold off;
figure();
for i = 1:10
    plot(t,Z(i,:,2),'-*', 'Color',colors(i,:));
    if i == 1
        xlabel('Time');
        ylabel('X_2');
        xlim([0,10]);
    end
    hold on;
end
saveas(gcf,'YCoordinatePlot-NonLinearNonAutonomous-Unstable','jpg')
hold off;
figure();
for i = 1:10
    plot(t,Z(i,:,1),'-*', 'Color',colors(i,:));
    if i == 1
        xlabel('Time');
        ylabel('X_1');
        xlim([0,10]);
    end
    hold on;
end
saveas(gcf,'XCoordinatePlot-NonLinearNonAutonomous-Unstable','jpg')
hold off;
%% find contraction metric in model based way
degree = 1;
[M,C,x,y,rho_2] =  find_contraction_metric(degree,G);
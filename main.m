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
global name_of_system
global t_prev
global interval
global trajectories
interval = 0.01;
name_of_system = 'unstable-nonlinear'; %% only use this name or 'unstable-linear'
trajectories = 6;
data_points = 60;
x_glo = zeros(data_points,5);
global index
global update
index = 1;
G = eye(2,2);
% G = [-0.7826 0.7731; -0.5110 0.033];
% G = [0;1];
tspan = 0:interval:1;          % excitation time
for i=1:trajectories
    % initial = [cos(pi*rand(1,1))*100*rand(1,1);cos(pi*rand(1,1))*100*rand(1,1)];    % initial state   for some systems, where to start matters
    initial = rand(2,1);    % initial state   for some systems, where to start matters
    t_prev = 0;
    update = true;
    [t, x] = ode15s(@(t,x) sys_poly(t,x,0,0,data_points), tspan, initial);
end
%% Corrupt the state derivative with noise
sup_cord_1 = max(abs(x_glo(1:end,4)));
sup_cord_2 = max(abs(x_glo(1:end,5)));
eps = max(sup_cord_1,sup_cord_2)/15;
for i = 1:data_points
    x_glo(i,4) = x_glo(i,4) + (eps)*(2*rand(1) - 1);
    x_glo(i,5) = x_glo(i,5) + (eps)*(2*rand(1) - 1);
end
%% build A,B,xi
T = data_points;      % of sample for design, increase if noise large
[phi_data,x_dot,data_fit] = data_matrix(T,eps);
%% clear yalmip
yalmip('clear')
%% find contraction metric in data driven way
tic
degree = 0;    %% this is the variable used to indicate the maximum degree of the feedback controller
degree_rho = 1;
[x,y,mu,rho_1,Q,W] = ddc_poly(phi_data,T,data_fit,degree,degree_rho,G);
toc
%% verify contraction metric in oringal model
verify_contraction_metric(x,y,rho_1,Q,G,0)
%% Create trajectories under the contraction feedback controller found
Z = [];
colors = lines(10);
for i = 1:10
    'hello'
    tspan = 0:0.01:1;
    initial = [cos(pi*rand(1,1))*10*rand(1,1);cos(pi*rand(1,1))*10*rand(1,1)];
    [t,Z(i,:,:)] = ode15s(@(t,z) sys_poly(t,z,eps,1,data_points,x,rho_1,G,Q), tspan, initial);
end
%% evaluate the feedback control action to show all trajectories are contracting
set(groot, ['Default', 'Line', 'LineWidth'], 2);
set(groot, 'DefaultAxesFontSize', 14);
set(groot, 'DefaultAxesFontWeight', 'bold');
figure();
for i = 1:10
    plot(t,Z(i,:,1),'-*', 'Color',colors(i,:));
    if i == 1
        xlabel('T');
        ylabel('X_1');
        xlim([0,1])
    end
    hold on;
end
saveas(gcf,'Plots/XCoordinatePlot-NonlinearNonAutonomous','jpg')
hold off
figure();
for i = 1:10
    plot(t,Z(i,:,2),'-*', 'Color',colors(i,:));
    if i == 1
        xlabel('T');
        ylabel('X_2');
        xlim([0,1])
    end
    hold on;
end
saveas(gcf,'Plots/YCoordinatePlot-NonLinearNonAutonomous','jpg')
hold off;
figure();
for i = 1:10
    plot(Z(i,:,1),Z(i,:,2),'-*', 'Color',colors(i,:));
    if i == 1
        xlabel('X_1');
        ylabel('X_2');
    end
    hold on;
end
saveas(gcf,'Plots/TrajectoriesPlot-NonlinearNonAutonomous','jpg')
hold off;

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
%% find contraction metric in model based way - THIS IS OPTIONAL JUST FOR DOUBLE CHECKING
degree = 2;
[M,C,x,y] =  find_contraction_metric(degree,G);
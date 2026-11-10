clear;
clc;
clf;

%% params
l0 = 0.08;
k = 200;
a = 0.08;
b = 0.08;
m = 1;
Ic = (1/6)*m*2*a*2*b;
g = -9.8;
theta = 23*pi/180;

k_list = 2*[100;100;100;100;100;100;100;100];

% %Global Positions of fixed spring points, P_world
P_world = [.0038354, .175; .0809284,.175;.175,-.0038354;.175,-.0809284; -.0038354,-.175;-.0809284,-.175;-.175,.0038354;-.175,.0809284]';
% 
% %Positions of spring points connected to inner box, P_box
P_box = [.0423819,.10489888;.0423819,.10489888;.10489888,-.0423819; .10489888,-.0423819; -.0423819,-.10489888;-.0423819,-.10489888; -.10489888,.0423819;-.10489888,.0423819]';

box_params = struct();
box_params.m = m;
box_params.I = Ic;
box_params.g = g;
box_params.k_list = k_list;
box_params.l0_list = l0*ones(length(k_list),1);
box_params.P_world = P_world;
box_params.P_box = P_box;


%% Find equilibrium state
my_rate_func = @(V) box_rate_func(0,V,box_params);

% Initial guess for equilibrium
V_guess = [0;0;0;0;0;0];
% thresholds
ftol = 1e-14; % ftol: termination threshold (stop when abs(f(x_{i}))<ftol
dxtol = 1e-14; % dxtol: termination threshold (stop when interval abs(x_{i+1}-x_i) < dxtol)
dxmax = 1e14;
max_iter = 100; % number of iterations per trial
num_iter = 100; % number of trials we would like to perform
[Veq, exit_flag] = multidim_newton_solver(my_rate_func,V_guess,dxtol,ftol, max_iter,dxmax, 1);
% disp(exit_flag); disp(Veq'); disp(norm(my_rate_func(Veq)));
    % equil x = 0.1231, equil y = 0.7501, equil theta = -0.1539

%% need to find jacobian

f = @(V) box_rate_func(0, V, box_params);
J_approx = approximate_jacobian(f,Veq);
A = J_approx;

dx0 = 2e-3; %2mm
dy0 = 2e-3; %2mm
dtheta0 = pi/180; % 1 degrees
vx0 = 2e-3; %2mm/s
vy0 = 2e-3; %2mm/s
vtheta0 = 0.01; % rad/s

%small number used to scale initial perturbation
% epsilon_list = [1];
% epsilon_list = [1, 0.1, 0.01, 0.001];
epsilon_list = [1, 10];
tspan = [0 1.5]; % 0 to 1 seconds
% initial perturbation
dV0 = [dx0; dy0; dtheta0; vx0; vy0; vtheta0];
my_rate_func = @(t,V) box_rate_func(t,V,box_params);
my_linear_rate = @(t,V) J_approx*(V-Veq);

% Plot formatting
LW_plot = 2;  % Line width
FS = 16;      % Font size

% for loop
for i = 1:length(epsilon_list)

    epsilon = epsilon_list(i);

    % Initial condition
    V0 = Veq + epsilon*dV0;

    % Integrate nonlinear system
    [t_nl,V_nl] = ode45(my_rate_func,tspan,V0);

    % Integrate linearized system
    [t_lin,V_lin] = ode45(my_linear_rate,tspan,V0);

    fprintf('Epsilon = %g\n',epsilon);
    fprintf('Initial Y = %.6f\n',V0(2));
    fprintf('Final nonlinear Y = %.6f\n',V_nl(end,2));
    fprintf('Final linearized Y = %.6f\n',V_lin(end,2));

    % Create figure
    figure('Position',[100 100 1000 500]);

    % X position
    subplot(3,1,1);
    plot(t_nl,V_nl(:,1),'b','LineWidth',LW_plot);
    hold on;
    plot(t_lin,V_lin(:,1),'r--','LineWidth',LW_plot);
    ylabel('$x$','Interpreter','latex','FontSize',FS);
    grid on;
    set(gca,'FontSize',FS,'TickLabelInterpreter','latex');

    % Y position
    subplot(3,1,2);
    plot(t_nl,V_nl(:,2),'b','LineWidth',LW_plot);
    hold on;
    plot(t_lin,V_lin(:,2),'r--','LineWidth',LW_plot);
    ylabel('$y$','Interpreter','latex','FontSize',FS);
    grid on;
    set(gca,'FontSize',FS,'TickLabelInterpreter','latex');

    % Orientation
    subplot(3,1,3);
    plot(t_nl,V_nl(:,3),'b','LineWidth',LW_plot);
    hold on;
    plot(t_lin,V_lin(:,3),'r--','LineWidth',LW_plot);
    ylabel('$\theta$','Interpreter','latex','FontSize',FS);
    xlabel('$t$ (s)','Interpreter','latex','FontSize',FS);
    grid on;
    set(gca,'FontSize',FS,'TickLabelInterpreter','latex');

    % Shared legend
    legend({'Nonlinear','Linearized'}, ...
        'Interpreter','latex', ...
        'FontSize',FS-1, ...
        'Location','best');

    % Overall title
    sgtitle(sprintf('Linearized vs. Nonlinear ODE: $\\epsilon = %g$%', ...
        epsilon), ...
    'Interpreter','latex','FontSize',FS+2);


end

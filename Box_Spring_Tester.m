

clc;
clf;
clear;

%Equilibrium length of all springs, lo = 0.08m
l0=.08;
%Spring Constant, k = 400N/m
k=100;
%horizontal distance between box's center and corner, a = 0.08m
a=.08;
%vertical distance between box's center and corner, b = 0.08m
b=.08;
%Mass of box, m = 1kg
m=1;
%Mass moment of inertia of box
I=(1/6)*m*2*a*2*b;
%Gravity = -9.8 m/s^2
g=-9.8;
%Initial rotation of box, θ = 23°
theta=23*pi/180;
%List of spring constants, k = 100N/m
k_list=[100;100;100;100;100;100;100;100];

%Global Positions of fixed spring points, P_world
P_world = [.0038354, .175; .0809284,.175;.175,-.0038354;.175,-.0809284; -.0038354,-.175;-.0809284,-.175;-.175,.0038354;-.175,.0809284];

%Positions of spring points connected to inner box, P_box
P_box = [.0423819,.10489888;.0423819,.10489888;.10489888,-.0423819; .10489888,-.0423819; -.0423819,-.10489888;-.0423819,-.10489888; -.10489888,.0423819;-.10489888,.0423819];

box_params = struct();
box_params.m = m;
box_params.I = I;
box_params.g = g;
box_params.k_list = k_list;
box_params.P_world = P_world;
box_params.l0_list = l0.*ones(length(k_list),1);
box_params.P_box = P_box;

t=0;
t_max=30;
tspan=[0 t_max];
t_stepcount=3000;
t_range=linspace(0,t_max,t_stepcount);
h_ref=(t_max+1)/t_stepcount;

%Initial C.O.M. position of box
x=0;
y=0;
%Initial box velocity
dxdt=0;
dydt=0;
%Initial box angular velocity
dthetadt=0;
%State vector
V0=[x;y;theta;dxdt;dydt;dthetadt];
F_list=zeros(length(k_list),2);

%% VAN DER HOUWEN'S/WRAY'S METHOD!!!
BT_struct_V = struct();
BT_struct_V.A = [0 0 0; 8/15 0 0; 1/4 5/12 0];
BT_struct_V.B = [1/4; 0; 3/4];
BT_struct_V.C = [0; 8/15; 2/3];

% %% Initialize Simulation
% num_zigs = 5;
% w = .1;
% hold on;
% spring_plot_struct = initialize_spring_plot(num_zigs,w);
% axis equal; axis square;
% axis([-.5,.5,-.5,.5]);
% 
% drawnow;


[t_V,X_V,h_avg_V, num_evals_V] = Box_RK_fixed_step_integration(@box_rate_func,tspan,V0,h_ref,BT_struct_V, box_params);
figure();
plot(t_V,X_V(:,1))
figure();
plot(t_V,X_V(:,2))



% box_corners=P_box(1:2:end,:);
% spring_links=box_params.P_box+box_params.P_world;

% for n = 1:length(k_list)
%     PA_n=P_box(n,:);
%     PB_n=P_world(n,:);
%     F_list(n,:)=compute_spring_force(k,l0,PA_n,PB_n);
% end
% 
% dVdt=box_rate_func(t,V0,box_params,F_list);








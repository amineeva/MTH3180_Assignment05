%Computes the linear and angular acceleration of the box
%given its current position and orientation
%INPUTS:
%x: current x position of the box
%y: current y position of the box
%theta: current orientation of the box
%box_params: a struct containing the parameters that describe the system
%Fields:
%box_params.m: mass of the box
%box_params.I: moment of inertia w/respect to centroid
%box_params.g: acceleration due to gravity
%box_params.k_list: list of spring stiffnesses
%box_params.l0_list: list of spring natural lengths
%box_params.P_world: 2 x n list of static mounting
% points for the spring (in the world frame)
%box_params.P_box: 2 x n list of mounting points
% for the spring (in the box frame)
%
%OUTPUTS
%ax: x acceleration of the box
%ay: y acceleration of the box
%atheta: angular acceleration of the box
function [ax,ay,atheta] = compute_accel(x,y,theta,box_params)
    P_now = compute_rbt(x,y,theta,box_params.P_box);   % 2 x n
    F_net = [0; -box_params.m*box_params.g];
    torque = 0;
    for i = 1:length(box_params.k_list)
        PA = box_params.P_world(:,i);   % fixed end
        PB = P_now(:,i);                % box end
        F = compute_spring_force(box_params.k_list(i), box_params.l0_list(i), PA, PB);
        r = PB - [x;y];
        F_net  = F_net + F;
        torque = torque + (r(1)*F(2) - r(2)*F(1));
    end
    ax = F_net(1)/box_params.m;
    ay = F_net(2)/box_params.m;
    atheta = torque/box_params.I;
end

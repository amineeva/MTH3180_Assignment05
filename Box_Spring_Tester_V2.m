function Box_Spring_Tester_V2()
    
    %Equilibrium length of all springs, lo = 0.08m
    l0=.08;
    %Spring Constant, k = 400N/m
    k=200;
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
    k_list=2*[100;100;100;100;100;100;100;100];
    
    %Global Positions of fixed spring points, P_world
    P_world = [.0038354, .175; .0809284,.175;.175,-.0038354;.175,-.0809284; -.0038354,-.175;-.0809284,-.175;-.175,.0038354;-.175,.0809284];
    
    %Positions of spring points connected to inner box, P_box
    P_box = [.0423819,.10489888;.0423819,.10489888;.10489888,-.0423819; .10489888,-.0423819; -.0423819,-.10489888;-.0423819,-.10489888; -.10489888,.0423819;-.10489888,.0423819];
    boundary_pts = [.0423819,.10489888; .10489888,-.0423819; -.0423819,-.10489888; -.10489888,.0423819;.0423819,.10489888];
    
    box_params = struct();
    box_params.m = m;
    box_params.I = I;
    box_params.g = g;
    box_params.k_list = k_list;
    box_params.P_world = P_world;
    box_params.l0_list = l0.*ones(length(k_list),1);
    box_params.P_box = P_box;
    
    t=0;
    t_max=10;
    tspan=[0 t_max];
    t_stepcount=1000;
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

    rate_func = @(tin,Vin) box_rate_func(tin,Vin,box_params);

    tlist = linspace(tspan(1),tspan(end),1001);

    [~,Vlist] = ode45(rate_func,tlist,V0);

    %% VIDEO STUFF %%
fname='Box_Spring_Vid.avi';
mypath1 = 'C:\Users\msmith2\';
mypath2 = 'Downloads';

input_fname = [mypath1, mypath2,fname];
writerObj = VideoWriter(fname,'Motion JPEG AVI');
open(writerObj);
    %% Initialize Simulation (also setup for video)
    fig1 = figure(1);
    set(fig1,'units','pixels','position',[0 0 1440 1080]);

    num_zigs = 5;
    w = .01;
    hold on;
    title('Vibrating Box Simulation','Interpreter','latex',FontSize=28)
    xlabel('X-position (meters)','Interpreter','latex',FontSize=20)
    ylabel('Y-position (meters)','Interpreter','latex',FontSize=20)
    
    spring_plots = {};
    box_plot = plot(0,0,'k','linewidth',2);

    num_springs = size(P_box,1);

    for n = 1:num_springs
        spring_plots{n} = initialize_spring_plot(num_zigs,w);
    end

    axis equal; axis square;


    % Vlist = 0*Vlist;

    for m = 1:length(tlist)
        
        P_box_world = compute_rbt(Vlist(m,1),Vlist(m,2),Vlist(m,3),P_box);
        boundary_pts_world = compute_rbt(Vlist(m,1),Vlist(m,2),Vlist(m,3),boundary_pts);

        set(box_plot,'xdata',boundary_pts_world(:,1),'ydata',boundary_pts_world(:,2));
        for n=1:num_springs
            update_spring_plot(spring_plots{n},P_world(n,:)',P_box_world(n,:)');
        end

        drawnow

        %% VIDEO STUFF 2, ELECTRIC BOOGALOO
        current_frame = getframe(fig1);
        writeVideo(writerObj,current_frame);
        pause(.01)
    end
%% Close Video
close(writerObj)

    % axis([-.5,.5,-.5,.5]);

    % drawnow;
    % % 
    % % 
    % plot(t_V,X_V(:,1))
    % figure();
    % plot(t_V,X_V(:,2))

end








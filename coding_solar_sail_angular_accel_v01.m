clc; clear; close all;

%% coding solar sail angular accel as a funct of payload mass cog movement using inchworm motors
% main parameters to adjust are dx_mp (x movement of payload) and
% theta_sol_max (angle of incoming solar radiation, also angle at 
% which grav force acts)
%%
F_sol_orig = 9.072e-6; % [N], solar radiation force, here approx 9uN
Y = 2; % [m], total y distance between sail and payload cg, here held constant
X = 0; % [m], total initial x distance between sail and payload cg

% solar sail (ss) constants
A_ss = 1; % [m^2], Area of solar sail, here 1m^2
r_ss = sqrt(A_ss/pi()); % [m], radius of solar sail
t_ss = 4.5e-6; % [m], thickness of mylar solar sail
rho_ss = 1390; % [kg/m^3], density of mylar solar sail 
m_ss = A_ss*t_ss*rho_ss; % [kg], mass of solar sail

% carbon fiber (cf) constants [todo: integrate these values into torque]
rho_cf = 2250; % [kg/m^3], density of cf composite 
r_cf = 150e-6; % [m], radius of carbon fiber pultruded composite, here 300um
A_cf = pi()*r_cf^2; % [m^2], area of carbon fiber pultruded composite wire
C_cf = 2*pi()*(r_ss)^2; % [m], circumference of solar sail where cf is surrounding
L_cf = 2; % [m], length of carbon fiber lines for control
V_cf_ring = A_cf*C_cf; % [m^3], volume of cf around ss 
V_cf_length = A_cf*4*L_cf; % [m^3], volume of cf lengths/tethers
% m_cf_ring = rho_cf*V_cf_ring; 
% m_cf_length = rho_cf*V_cf_length;
% m_cf = m_cf_ring + m_cf_length;

% currently just setting cf masses to zero until I fix the eqtns here
m_cf_ring = 0;
m_cf_length = 0;
m_cf = 0;

% mass of payload (mp) constants and payload mass fraction
m_mp = 0.010; % [kg], projected payload mass, here 10 grams
m_tot = m_mp + m_ss + m_cf; % [kg], total mass
nu = m_mp/m_tot; % [kg/kg], mass efficiency (aka mass fraction) of payload

% ss efficiency and forces
n = 0.9; % solar sail efficiency
P_rad = 5.04e-6; % [N/m^2], solar radiation @ 1 AU
F_sol_0 = 2*n*P_rad*A_ss; % [N], force on solar sail due to rad pressure

% characteristic acceleration
a_c = F_sol_0/((1-nu)*m_tot); % [], characteristic acceleration

% solar constants for grav forces
G = 6.674e-11; % [m3 kg-1 s-2] gravitational constant;
m_sol = 1.989e30; % [kg] mass of sun;
r_sol = 150e9; % [m] distance from earth to sun;
mu = G*m_sol/r_sol^2;

% grav forces
F_g_ss = m_ss*mu; % force due to grav from sun on ss
F_g_mp = m_mp*mu; % force due to grav from sun on mp

num = 101; % counter used to set max of arrays

% movement by inchworm motors, but keeping x_cg at zero and moving ss and
% mp instead
dx_mp = linspace(0,10,num)/1000; % [m], change in x position of payload mass
x_mp =  (1 - nu)*dx_mp; % [m], x seperation of mp from cg
x_ss = -nu*dx_mp; % [m], x Seperation of ss from cg
x_cg = x_ss/nu + x_mp/(1-nu); % [m], x position of cg, here trying to set to zero

% y position of ss and mp using mass fraction
dY = Y*ones(1,num);
y_mp = -(1-nu)*dY; % [m], y position of payload
y_ss = nu*dY; % [m], y position of solar sail
y_cg = zeros(1,num); % [m], y position of cg

% initializing values
x_mp_0 = ones(num,num);
y_mp_0 = ones(num,num);

Torque = zeros(num,num);
alpha = zeros(num,num);
omega = zeros(num,num);
theta = zeros(num,num);

t_max = 3000; % [s], maximum time span

for i = 1:num % loop counting through x_mp
    
    t = linspace(0,t_max,num); % [s], time step to see progress as updated from forces on bliss
    dt = t(2) - t(1);
    
    % resetting initial anglular values
    alpha_0 = 0; % [deg/s/s], original angular accel of bliss, here set to zero
    omega_0 = 0; % [deg/s], original angular speed of bliss, here set to zero
    theta_sol_max = 0; % [deg], maximum angle that theta_sol will sweep through to
    theta_sol = linspace(0,theta_sol_max,num); % [deg], angle of sun line
    theta_sol_orig = theta_sol(num);
    
    % mass moment of inertia (resistance to rotation)
    I_ss_z(i) = m_ss*(y_ss(i)^2 + x_ss(i)^2 + 0.25*pi()*r_ss^2);
    I_mp_z(i) = m_mp*(x_mp(i)^2 + y_mp(i)^2);
    I_cf_z(i) = m_cf_ring/2*r_ss^2 + 4*m_cf_length/12*L_cf^2;
    I_tot(i) = I_ss_z(i) + I_mp_z(i) + I_cf_z(i); 

    for j = 1:num % loop for stepping through time
        
        theta_0(j) = theta_sol(num); % need theta_0 to be current solar angle
        
        % here we assume that the back half of the ss is not reflective but
        % still has an absorbance efficinency equal to the reflective
        % efficiency of n = 0.90
        if mod(abs(theta_0(j)),270) > 90
            F_sol(i,j) = n*P_rad*A_ss;
        elseif mod(abs(theta_0(j)),270) <= 90
            F_sol(i,j) = 2*n*P_rad*A_ss;
        end

        % Force due to sol on ss in cartesian
        F_x_sol = F_sol(i,j)*sind(theta_sol(num))*cosd(theta_sol(num));
        F_y_sol = F_sol(i,j)*cosd(theta_sol(num));%*cosd(theta_sol(num));
        
        % Force due to grav on ss in cart
        F_x_g_ss = -F_g_ss*sind(theta_sol(num));
        F_y_g_ss = -F_g_ss*cosd(theta_sol(num));
        
        % Force due to grav on mp in cart
        F_x_g_mp = -F_g_mp*sind(theta_sol(num));
        F_y_g_mp = -F_g_mp*cosd(theta_sol(num));
        
        % Torsion due to forces on ss, mp, and cf[todo]
        T_sol(i,j) = F_x_sol*y_ss(i) - F_y_sol*x_ss(i);
        T_g_ss(i,j) = F_x_g_ss*y_ss(i) - F_y_g_ss*x_ss(i);
        T_g_mp(i,j) = F_x_g_mp*y_mp(i) - F_y_g_mp*x_mp(i);
        % torsion due to grav on cf [todo]
%         T_g_cf(i,j) = F_g_cf*sind(theta_sol(num))*abs(y_cf(i)) ...
%                 + F_g_cf*cosd(theta_sol(num))*abs(x_cf(i));
        T_g_cf(i,j) = 0; % setting to zero until i fix above
        T_tot(i,j) = T_sol(i,j) + T_g_ss(i,j) + T_g_mp(i,j) + T_g_cf(i,j);
            
        % rotational dynamics
        alpha(i,j) = alpha_0 + rad2deg(T_tot(i,j)/I_tot(i));
        omega(i,j) = omega_0 + alpha(i,j)*dt;
        theta(i,j) = theta_0(j) + omega(i,j)*dt + 0.5*(alpha(i,j)*dt^2);
        
        % updating initial rotational dynamics
        omega_0 = omega(i,j);
        theta_sol(num) = theta(i,j);

%         % linear forces
%         F_x_tot = F_x_sol + F_x_g_ss + F_x_g_mp;
%         F_y_tot = F_y_sol + F_y_g_ss + F_y_g_mp;
%             
%         % linear accelerations
%         a_x_lin = a_0x + F_x_tot/m_tot;
%         a_y_lin = a_0y + F_y_tot/m_tot;
%         
%         % linear speeds
%         v_x_lin(j) = v_0x + a_x_lin*t(j);
%         v_y_lin(j) = v_0y + a_y_lin*t(j);
%         
%         % linear displacements
%         x_tot(j) = x_0 + v_x_lin(j)*t(j) + a_x_lin*t(j)^2;
%         y_tot(j) = y_0 + v_y_lin(j)*t(j) + a_y_lin*t(j)^2;
    end
end

figure(201)
subplot(2,3,1)
ps1 = surf(x_mp*10^3,t,alpha');
title(strcat('Angular Acceleration, dx_m_p = ', sprintf('%.0f',dx_mp(num)*1e3),'mm,   \theta_s_o_l_,_i_n_i_t = ',sprintf('%.0f',theta_sol_orig),' deg'))
xlabel('x_m_p (mm)')
ylabel('t (s)')
zlabel('\alpha_a_c_c_e_l (deg/s/s)')
c = colorbar;
% c.Label.String = 'm_p (kg)';
% set(gca, 'xdir', 'reverse');
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
% set(gca, 'ZScale', 'log')
% set(gca,'ColorScale','log')
ps1.EdgeColor = 'none';
ax = gca;
ax.View = [45 45];

% figure(202)
subplot(2,3,2)
ps1 = surf(x_mp*10^3,t,omega');
title(strcat('Angular Speed, dx_m_p = ', sprintf('%.0f',dx_mp(num)*1e3),'mm,   \theta_s_o_l_,_i_n_i_t = ',sprintf('%.0f',theta_sol_orig),' deg'))
xlabel('x_m_p (mm)')
ylabel('t (s)')
zlabel('\omega (deg/s)')
c = colorbar;
% c.Label.String = 'm_p (kg)';
% set(gca, 'xdir', 'reverse');
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
% set(gca, 'ZScale', 'log')
% set(gca,'ColorScale','log')
ps1.EdgeColor = 'none';
ax = gca;
ax.View = [45 45];

% figure(203)
subplot(2,3,3)
ps1 = surf(x_mp*10^3,t,theta');
title(strcat('Angle, dx_m_p = ', sprintf('%.0f',dx_mp(num)*1e3),'mm,   \theta_s_o_l_,_i_n_i_t = ',sprintf('%.0f',theta_sol_orig),' deg'))
xlabel('x_m_p (mm)')
ylabel('t (s)')
zlabel('\theta (deg)')
c = colorbar;
% c.Label.String = 'm_p (kg)';
% set(gca, 'xdir', 'reverse');
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
% set(gca, 'ZScale', 'log')
% set(gca,'ColorScale','log')
ps1.EdgeColor = 'none';
ax = gca;
ax.View = [45 45];

x_low = 2;
% figure(301)
subplot(2,3,4)
plot(t,alpha(num,:),t,alpha(x_low,:), '--r')
title(strcat('Angular Accel, \theta_s_o_l_,_i_n_i_t = ',sprintf('%.0f',theta_sol_orig),' deg'))
xlabel('t (s)')
ylabel('\alpha_a_c_c_e_l (deg/s/s)')
legend(strcat('dx_m_p = ', sprintf('%.3f',dx_mp(num)*1e3),'mm'), strcat('dx_m_p = ', sprintf('%.3f',dx_mp(x_low)*1e3),'mm')) 

% figure(302)
subplot(2,3,5)
plot(t,omega(num,:),t,omega(x_low,:), '--r')
title(strcat('Angular Speed, \theta_s_o_l_,_i_n_i_t = ',sprintf('%.0f',theta_sol_orig),' deg'))
xlabel('t (s)')
ylabel('\omega (deg/s)')
legend(strcat('dx_m_p = ', sprintf('%.3f',dx_mp(num)*1e3),'mm'), strcat('dx_m_p = ', sprintf('%.3f',dx_mp(x_low)*1e3),'mm')) 

% figure(303)
subplot(2,3,6)
plot(t,theta(num,:),t,theta(x_low,:), '--r')
title(strcat('Clock Angle, \theta_s_o_l_,_i_n_i_t = ',sprintf('%.0f',theta_sol_orig),' deg'))
xlabel('t (s)')
ylabel('\theta (deg)')
legend(strcat('dx_m_p = ', sprintf('%.3f',dx_mp(num)*1e3),'mm'), strcat('dx_m_p = ', sprintf('%.3f',dx_mp(x_low)*1e3),'mm')) 


% figure(305)
% plot(t, T_tot(num,:), t, T_sol(num,:),'--r', t, T_g_ss(num,:), '--k', t, T_g_mp(num,:),'--b')
% title(strcat('Torque, dx_m_p = ', sprintf('%.0f',dx_mp(num)*1e3),'mm,   \theta_s_o_l_,_i_n_i_t = ',sprintf('%.0f',theta_sol_orig),' deg'))
% xlabel('t (s)')
% ylabel('\tau (Nm)')
% legend('\tau_t_o_t', '\tau_s_o_l', '\tau_g_s_s', '\tau_g_m_p')
% 
% figure(304)
% plot(x_mp,I_tot)
% xlabel('t (s)')
% ylabel('I_t_o_t (kgm^2)')
% 
% figure(306)
% surf(x_mp,t, F_sol)
% xlabel('x')
% ylabel('t')
% zlabel('F_s_o_l')

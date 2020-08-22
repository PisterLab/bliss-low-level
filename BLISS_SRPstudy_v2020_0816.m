%% Berkeley Low-cost Interplanetary Solar Sail (BLISS) project
% This file is used to study the effect of solar radiation pressure on a
% spacecraft initially starting in a geosynchronous (GSO) orbit.
% There are four bodies with the following forces acting upon them:
% Sun   : assumed to be fixed in heliocentric inertial frame
% Earth : Sun's gravity
% Sail  : Sun's gravity, Earth's gravity, force from Solar Radiation Pressure (SRP)
% NEO   : Sun's gravity, Earth's gravity
%% Clear variables, workspace, figures
clear all; clc; clf;
%% Physical constants
m_AU     = (149597870.7)*10^3;                         % [m/AU]     Length of an Astronomical Unit in meters
mu_Sun   = (1.3271244004193938e11)*10^9;               % [m^3/s^2]  Standard gravitational parameter of Sun
mu_Earth = (398600.4418)*10^9;                         % [m^3/s^2]  Standard gravitational parameter of Earth
rE       = (6378.137)*10^3;                            % [m]        Earth equatorial radius
rGEO     = rE + (35786*10^3);                          % [m]        Distance of GSO orbit from Earth's center
%% Satellite parameters (10 gram, 1 square meter area sail)
m        = 0.01;                                       % [kg]       Mass of satellite
A_SRP    = 1.00;                                       % [m^2]      SRP effective area
C_SRP    = 2.00;                                       % [ ]        SRP coefficent, C_SRP = 1 for perfectly absorbing, C_SRP = 2 for perfectly reflective  solar sail
%% Solar Radiation Pressure (SRP) as a function of distance from Sun
% values collected from Wikipedia page on SRP
Rvec = [0.20 0.39 0.72 1.00 1.52 3.00 5.20];           % [AU]       Distance from Sun
Pvec = [227.0 59.7 17.5 9.08 3.93 1.01 0.34]/2*10^-6;  % [N/m^2]    Radiation pressure
% looks exponential, create an exponential fit
fit = @(x) x(1).*exp(x(2).*Rvec) + x(3).*exp(x(4).*Rvec);
[xsol,fval] = fminsearch(@(x) norm(Pvec-fit(x)), [5e-4; 0; -1.5e-4; -1.0] );
Pvecfit = @(R) xsol(1)*exp(xsol(2).*R) + xsol(3)*exp(xsol(4).*R);
figure(1)
lnwidth = 2;
fsize   = 14;
clf;
hold on;
grid on;
plot(Rvec,Pvec,'o','linewidth',2)
Rvectry = 0.1:0.1:6;
plot(Rvectry,Pvecfit(Rvectry),'linewidth',2)
xlabel('Distance from Sun [AU]'); ylabel('Solar Radiation Pressure [N/m^2]');
title('SRP as a function of Distance from Sun')
legend('Data','Fit')
set(gca, 'FontSize', fsize,'FontWeight','bold')
P_SRP = @(r) Pvecfit(r/m_AU);   % function of distance in meters
%% Initial conditions and ODE solver settings
r_Sun_Earth = m_AU;
x0_Earth    = [r_Sun_Earth;0;0;0;sqrt(mu_Sun/r_Sun_Earth);0];
Rot         = @(a) [cos(a) sin(a) 0;-sin(a) cos(a) 0;0 0 1];   % x0_NEO = [m_AU+100*rGEO;0;0;0;-1.3*sqrt(mu_Sun/(m_AU+100*rGEO));0];
x0_Sail     = kron(eye(2),Rot(0.0*pi/180))*[r_Sun_Earth+100*rGEO;0;0;0;sqrt(mu_Sun/(r_Sun_Earth+100*rGEO));0];
% x0_NEO      = kron(eye(2),Rot(42*pi/180))*[-7.470427658352538e+11;2.078004922068783e+11;0;0;0;0];
x0_NEO      = kron(eye(2),Rot(42*pi/180))*[-7.470427658352538e+11;2.078004922068783e+11;0;0;0;0];
x0          = [x0_Earth;x0_Sail;x0_NEO];
dt          = 200;
tEnd        = 2*365*86164;%690*86164;
options     = odeset('AbsTol',1e-12,'RelTol',1e-12);
%% Solve for trajectories
tStart = tic;
[tSim,xSim] = ode23s(@(t,x) fsim(t,x,mu_Sun,mu_Earth,m,P_SRP,C_SRP,A_SRP),0:dt:tEnd,x0,options);
toc(tStart)
%% Extract Inputs

Plot_inputs = 1;   % set to 1 to plot input figure. Takes some seconds.
if Plot_inputs
    tStart = tic;
    numt = numel(tSim);
    [xdotSim,avec_Sun_EarthSim,avec_Sun_SailSim,avec_Sun_NEOSim,avec_Earth_SailSim,avec_Earth_NEOSim,avec_SRP_SailSim] = cellfun(@(t,x) fsim(t,x,mu_Sun,mu_Earth,m,P_SRP,C_SRP,A_SRP), num2cell(tSim), num2cell(xSim',1)','uni',0);
    avec_Sun_EarthSim_mat  = reshape(cell2mat(avec_Sun_EarthSim),3,numt);
    avec_Sun_SailSim_mat   = reshape(cell2mat(avec_Sun_SailSim),3,numt);
    avec_Sun_NEOSim_mat    = reshape(cell2mat(avec_Sun_NEOSim),3,numt);
    avec_Earth_SailSim_mat = reshape(cell2mat(avec_Earth_SailSim),3,numt);
    avec_Earth_NEOSim_mat  = reshape(cell2mat(avec_Earth_NEOSim),3,numt);
    avec_SRP_SailSim_mat   = reshape(cell2mat(avec_SRP_SailSim),3,numt);
    
    xdotSim_mat = cell2mat(xdotSim);
    
    dx = reshape(xdotSim_mat',18,numt);
    dx_Earth    = dx(1:6,:);
    dx_Sail     = dx(7:12,:);
    dx_NEO      = dx(13:18,:);
    atot_Earth  = dx_Earth(4:6,:);
    atot_Sail   = dx_Sail(4:6,:);
    atot_NEO    = dx_NEO(4:6,:);
    atot_Earth_mag = vecnorm(atot_Earth,2,1);
    atot_Sail_mag  = vecnorm(atot_Sail,2,1);
    atot_NEO_mag = vecnorm(atot_NEO,2,1);
    
    a_grav_Sun_Earth_mag  = vecnorm(avec_Sun_EarthSim_mat,2,1);
    a_grav_Sun_Sail_mag   = vecnorm(avec_Sun_SailSim_mat,2,1);
    a_grav_Sun_NEO_mag    = vecnorm(avec_Sun_NEOSim_mat,2,1);
    a_grav_Earth_Sail_mag = vecnorm(avec_Earth_SailSim_mat,2,1);
    a_grav_Earth_NEO_mag  = vecnorm(avec_Earth_NEOSim_mat,2,1);
    a_SRP_Sail_mag        = vecnorm(avec_SRP_SailSim_mat,2,1);
    toc(tStart)
    
    figure(2);
    clf;
    lnwidth = 2;
    fsize   = 14;
    plot(tSim/86164,a_grav_Sun_Sail_mag,'color',[243 130 53]/255,'linewidth',lnwidth);
    hold on
    plot(tSim/86164,a_grav_Earth_Sail_mag,'color',[0 0.447 0.741],'linewidth',lnwidth);
    plot(tSim/86164,a_SRP_Sail_mag,'color',[0.6350 0.0780 0.1840],'linewidth',lnwidth);
    plot(tSim/86164,atot_Sail_mag,'color','k','linestyle','--','linewidth',lnwidth);
    ytickformat('%.1f')
    xlabel('time [Days]'); ylabel('Specific Forces [m/s^2]');
    title('Specific Force Magnitudes acting on Sail')
    legend('Sun gravity','Earth gravity','SRP','Total','location','best')
    grid on;
    set(gca, 'FontSize', fsize,'FontWeight','bold')
    
    figure(3);
    clf;
    lnwidth = 2;
    fsize   = 14;
    plot(tSim/86164,a_grav_Sun_NEO_mag,'color',[243 130 53]/255,'linewidth',lnwidth);
    hold on
    plot(tSim/86164,a_grav_Earth_NEO_mag,'color',[0 0.447 0.741],'linewidth',lnwidth);
    % plot(tSim/86164,atot_NEO_mag,'color','k','linewidth',lnwidth);
    ytickformat('%.1f')
    xlabel('time [Days]'); ylabel('Specific Forces [m/s^2]');
    title('Specific Force Magnitudes acting on NEO')
    legend('Sun gravity','Earth gravity','location','best')
    grid on;
    set(gca, 'FontSize', fsize,'FontWeight','bold')
    
    figure(4);
    clf;
    lnwidth = 2;
    fsize   = 14;
    plot(tSim/86164,a_grav_Sun_Earth_mag,'color',[243 130 53]/255,'linewidth',lnwidth);
    % hold on
    % plot(tSim/86164,atot_Earth_mag,'color','k','linewidth',lnwidth);
    % ytickformat('%.1f')
    xlabel('time [Days]'); ylabel('Specific Forces [m/s^2]');
    title('Specific Forces acting on Earth')
    legend('Sun gravity','location','best')
    grid on;
    set(gca, 'FontSize', fsize,'FontWeight','bold')
end
%% Closest Point of Approach
rvec_Sun_Sail     = xSim(:,7:9);
rvec_Sun_NEO      = xSim(:,13:15);
rvec_Sail_NEO     = rvec_Sun_NEO - rvec_Sun_Sail;
r_Sail_NEO        = vecnorm(rvec_Sail_NEO,2,2);
[CPA_Sail_NEO,t_CPA_Sail_NEO] = min(r_Sail_NEO);
disp(['Closest point of approach to NEO is ',num2str(CPA_Sail_NEO/1000,'%0.1f'),' [km] at Day ',num2str(tSim(t_CPA_Sail_NEO)/86164,'%0.1f')])
%% Animation
xtrajSim = xSim(:,1:6:end)';
ytrajSim = xSim(:,2:6:end)';
ztrajSim = xSim(:,3:6:end)';

recordAnimation2File = 0;   % set to 1 to save a movie file to the animation subfolder
figure(99)
clf;
animateSolarSystem(xtrajSim,ytrajSim,ztrajSim,tSim,tSim(end),recordAnimation2File)
%% Model
function [xdot,avec_Sun_Earth,avec_Sun_Sail,avec_Sun_NEO,avec_Earth_Sail,avec_Earth_NEO,avec_SRP_Sail] = fsim(t,x,mu_Sun,mu_Earth,m,P_SRP,C_SRP,A_SRP)
% States
rvec_Sun_Earth    = x(1:3);
vvec_Sun_Earth    = x(4:6);
rvec_Sun_Sail     = x(7:9);
vvec_Sun_Sail     = x(10:12);
rvec_Sun_NEO      = x(13:15);
vvec_Sun_NEO      = x(16:18);

% Absolute distances w.r.t Sun
r_Sun_Earth       = norm(rvec_Sun_Earth);
r_Sun_Sail        = norm(rvec_Sun_Sail);
r_Sun_NEO         = norm(rvec_Sun_NEO);

% Relative positions and distances
rvec_Earth_Sail   = rvec_Sun_Sail - rvec_Sun_Earth;
rvec_Earth_NEO    = rvec_Sun_NEO - rvec_Sun_Earth;
r_Earth_Sail      = norm(rvec_Earth_Sail);
r_Earth_NEO       = norm(rvec_Earth_NEO);

% Sun's gravity on Earth, Sail, and NEO
avec_Sun_Earth    = -mu_Sun/r_Sun_Earth^3 * rvec_Sun_Earth;
avec_Sun_Sail     = -mu_Sun/r_Sun_Sail^3  * rvec_Sun_Sail;
avec_Sun_NEO      = zeros(3,1);%-mu_Sun/r_Sun_NEO^3   * rvec_Sun_NEO;

% Earth's gravity on Sail and NEO
avec_Earth_Sail   = -mu_Earth/r_Earth_Sail^3 * rvec_Earth_Sail;
avec_Earth_NEO    = -mu_Earth/r_Earth_NEO^3 * rvec_Earth_NEO;

% SRP on Sail (A_SRP is a constant 1.0 m^2 right now but through control of
% sail orientation, we can have this vary between 0 <= A_SRP <= 1.
avec_SRP_Sail     = (1/m) * P_SRP(r_Sun_Sail) * C_SRP * A_SRP * rvec_Sun_Sail/r_Sun_Sail;

% State derivatives
rdotvec_Sun_Earth = vvec_Sun_Earth;
rdotvec_Sun_Sail  = vvec_Sun_Sail;
rdotvec_Sun_NEO   = vvec_Sun_NEO;
vdotvec_Sun_Earth = avec_Sun_Earth;
vdotvec_Sun_Sail  = avec_Sun_Sail + avec_Earth_Sail + avec_SRP_Sail;
vdotvec_Sun_NEO   = avec_Sun_NEO + avec_Earth_NEO;

xdot = [rdotvec_Sun_Earth;vdotvec_Sun_Earth;rdotvec_Sun_Sail;vdotvec_Sun_Sail;rdotvec_Sun_NEO;vdotvec_Sun_NEO];
end
%% E2E sim

%% clear everything
clear; close all; clc;

%% DEFINE INTIAL CONDITIONS AND RELEVANT CONSTANTS
tol = 1e-10;

addpath('./mean_osc/');
addpath('../Documents/Github/aa279d-project/helper_functions/');
addpath('../Documents/Github/aa279d-project/pset1');

% Earth physical constants from Vallado
muearth = 3.986 * 10^5;              % [km^3/sec^2]
reearth =   6378.1363;               % [km] mean equatorial radius
J2earth =      0.0010826266;         % J2

% determine orbital elements
SV4_init_QNSOE = [6944, -4e-5, 1.6e-3, 99.4, -151.1, -47.9];
a_init_SV4 = SV4_init_QNSOE(1);
ex_init_SV4 = SV4_init_QNSOE(2);
ey_init_SV4 = SV4_init_QNSOE(3);
inclination_init_SV4 = SV4_init_QNSOE(4);
RAAN_init_SV4 = SV4_init_QNSOE(5);
u_init_SV4 = SV4_init_QNSOE(6);

% Determine classical OEs
[a_init_SV4, e_init_SV4, inclination_init_SV4, RAAN_init_SV4, w_init_SV4, f_init_SV4, M_init_SV4] = ...
    QNSOE2OE(a_init_SV4, ex_init_SV4, ey_init_SV4, inclination_init_SV4, RAAN_init_SV4, u_init_SV4);

% FOR THIS PSET FORCE e=0 CHANGE WHEN NECESSARY
% e_init_SV4 = 0;

n_init_Sv4 = sqrt(muearth / a_init_SV4^3); % rad/s
T_SV4 = 2*pi*sqrt(a_init_SV4^3/muearth);
nOrbits = 20;
t0 = 0;
tf = T_SV4*nOrbits;
step_size = 100;

%% DETERMINE CHIEF INITIAL CONDITIONS AND SETUP NUMERICAL INTEGRATION

% determine initial velocity and position 
[r_init_SV4, v_init_SV4] = OE2ECI(a_init_SV4, e_init_SV4, inclination_init_SV4, RAAN_init_SV4, w_init_SV4, f_init_SV4);
initial_state_SV4 = [r_init_SV4; v_init_SV4];

tspan = 0:step_size:tf; % [sec] 
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-12); % higher accuracy to see through numerical errors

%% DEFINE DEPUTY INITIAL CONDITIONS
SV3_init_ROE = [-1, -79328, 42, 452, 36, 827];

% should the a component be multiplied again?

% set exact a component to 0 for this problem set
SV3_init_ROE(1) = 0;

% assuming init in degrees here, can check e_x, e_y defs too
SV3_init_QNSOE = ROE2QNSOE(a_init_SV4, ex_init_SV4, ey_init_SV4, inclination_init_SV4, RAAN_init_SV4, ...
    M_init_SV4 + w_init_SV4, SV3_init_ROE(1), SV3_init_ROE(2), SV3_init_ROE(3), ...
    SV3_init_ROE(4), SV3_init_ROE(5), SV3_init_ROE(6));

[a_init_SV3, e_init_SV3, inclination_init_SV3, RAAN_init_SV3, w_init_SV3, f_init_SV3, M_init_SV3] = ...
    QNSOE2OE(SV3_init_QNSOE(1), SV3_init_QNSOE(2), SV3_init_QNSOE(3), ...
             SV3_init_QNSOE(4), SV3_init_QNSOE(5), SV3_init_QNSOE(6));

[r_init_SV3, v_init_SV3] = OE2ECI(a_init_SV3, e_init_SV3, inclination_init_SV3, RAAN_init_SV3, w_init_SV3, f_init_SV3);

%% TRANSFORM DEPUTY TO CHIEF RTN

% define chief RTN frame
[rot_RTN, w_RTN] = ECI2RTN(r_init_SV4, v_init_SV4, a_init_SV4, e_init_SV4, f_init_SV4);

% take SV3 to chief RTN frame
[rho_init_SV3_RTN, rho_dot_init_SV3_RTN] = rvECI2rvRTN(r_init_SV3 - r_init_SV4, v_init_SV3 - v_init_SV4, rot_RTN, w_RTN);
[r_init_Sv4_RTN, v_init_Sv4_RTN]        = rvECI2rvRTN(r_init_SV4, v_init_SV4, rot_RTN, w_RTN);

%% RUN SIMULATION FOR SV3 AND SV4 IN THE ABSOLUTE FRAME

% determine initial velocity and position 
initial_state_SV4 = [r_init_SV4; v_init_SV4];
initial_state_SV3 = [r_init_SV3; v_init_SV3];

%% RUN NUMERICAL INTEGRATION
% simulate
[t_SV4, output_SV4_ECI] = ode45(@orbital_EOM_ECI_J2, tspan, initial_state_SV4, options);
[t_SV3, output_SV3_ECI] = ode45(@orbital_EOM_ECI_J2, tspan, initial_state_SV3, options);

% Convert to the chiefs RTN frame for each time step
r_SV4_ECI = [output_SV4_ECI(:,1)'; output_SV4_ECI(:,2)'; output_SV4_ECI(:,3)'];
v_SV4_ECI = [output_SV4_ECI(:,4)'; output_SV4_ECI(:,5)'; output_SV4_ECI(:,6)'];
r_SV3_ECI = [output_SV3_ECI(:,1)'; output_SV3_ECI(:,2)'; output_SV3_ECI(:,3)'];
v_SV3_ECI = [output_SV3_ECI(:,4)'; output_SV3_ECI(:,5)'; output_SV3_ECI(:,6)'];

% preallocate arrays for transformed positions and velocities
r_SV4_RTN = zeros(3, length(t_SV4));
v_SV4_RTN = zeros(3, length(t_SV4));
r_SV3_RTN = zeros(3, length(t_SV4));
v_SV3_RTN = zeros(3, length(t_SV4));

% Transform to the chiefs RTN frame
for i = 1:length(t_SV4)
    % find rotation from ECI to chiefs RTN frame
    [rot_RTN, w_RTN] = ECI2RTN(r_SV4_ECI(:, i), v_SV4_ECI(:, i), a_init_SV4, e_init_SV4, f_init_SV4);
    [r_SV4_RTN(:, i), v_SV4_RTN(:, i)] = rvECI2rvRTN(r_SV4_ECI(:, i), v_SV4_ECI(:, i), rot_RTN, w_RTN);
    [r_SV3_RTN(:, i), v_SV3_RTN(:, i)] = rvECI2rvRTN(r_SV3_ECI(:, i), v_SV3_ECI(:, i), rot_RTN, w_RTN);
end

% Find deputy with respect to chief orbit
rho_SV3_SV4_RTN = r_SV3_RTN - r_SV4_RTN;

%% PLOT RESULTS OF ABSOLUTE FRAME SIMULATION IN RTN FRAME
figure(1); clf;
sgtitle('Relative Motion in RTN Frame using absolute orbits');

subplot(2, 2, 1);
plot(rho_SV3_SV4_RTN(2, :), rho_SV3_SV4_RTN(1, :));
xlabel('Tangential (km)');
ylabel('Radial (km)');
axis equal;
grid on;

subplot(2, 2, 2);
plot(rho_SV3_SV4_RTN(3, :), rho_SV3_SV4_RTN(1, :));
xlabel('Normal (km)');
ylabel('Radial (km)');
axis equal;
grid on;

subplot(2, 2, 3);
plot(rho_SV3_SV4_RTN(3, :), rho_SV3_SV4_RTN(2, :));
xlabel('Normal (km)');
ylabel('Tangential (km)');
axis equal;
grid on;

subplot(2, 2, 4);
plot3(rho_SV3_SV4_RTN(2, :), rho_SV3_SV4_RTN(3, :), rho_SV3_SV4_RTN(1, :));
xlabel('Tangential (km)');
ylabel('Normal (km)');
zlabel('Radial (km)');
title('Relative Motion in 3D');
grid on;

%% FODE SETUP COMPLETE %%

%% Begin PS8 Material

% check initialization for roe vs aroe
% Desired ROE state 
aroe_des = [0, 7920, -297, -315, 16, 1144]; % [m]
roe_des  = aroe_des/(a_init_SV4*1e3);

%% Initializations and Setup
% inits of x0, P0, Q, R2
% setup values histories for all time steps

n = 6; % state dim

% allocate logs
Nsteps = numel(tspan);
x_plus = zeros(n, Nsteps);
x_minus = zeros(n, Nsteps);
P_plus = zeros(n, n, Nsteps);
P_minus = zeros(n, n, Nsteps);
K = zeros(n, n, Nsteps);
h = zeros(n, Nsteps);
y = zeros(n, Nsteps);
roeHist = zeros(n, Nsteps);
preResiduals = zeros(n, Nsteps);
postResiduals = zeros(n, Nsteps);
chiefOEs = zeros(n, Nsteps);

% control inits
trackingErrors = zeros(n, Nsteps);
controlInputs = zeros(2, Nsteps);
controlStateHist = zeros(Nsteps-1, 18);

% seed for now
seed = 42;
rng(seed);

% QNSOE2ROE for deputy ROEs - BIG BLOCK
SV3_OE = ECI2OE(r_SV3_ECI(:, 1), v_SV3_ECI(:, 1), muearth);
SV4_OE = ECI2OE(r_SV4_ECI(:, 1), v_SV4_ECI(:, 1), muearth);
ac = SV4_OE.a;
wc = wrapTo360(SV4_OE.w);
ex = SV4_OE.e * cosd(wc);
ey = SV4_OE.e * sind(wc);
inc = SV4_OE.i;
SV4_f = wrapTo360(SV4_OE.anom);
SV3_f = wrapTo360(SV3_OE.anom);
SV4_M = rad2deg(true2mean(deg2rad(SV4_f), SV4_OE.e));
SV3_M = rad2deg(true2mean(deg2rad(SV3_f), SV3_OE.e));
SV4_OE_vec = [SV4_OE.a*1e3, SV4_OE.e, deg2rad(SV4_OE.i), deg2rad(SV4_OE.Om), ...
              deg2rad(wc), true2mean(deg2rad(SV4_f), SV4_OE.e)];
SV3_OE_vec = [SV3_OE.a*1e3, SV3_OE.e, deg2rad(SV3_OE.i), deg2rad(SV3_OE.Om), ...
              deg2rad(wrapTo360(SV3_OE.w)), true2mean(deg2rad(SV3_f), SV3_OE.e)];
SV4_OE_mean_noUnits = osc2mean(SV4_OE_vec); % meters and rad
SV3_OE_mean_noUnits = osc2mean(SV3_OE_vec);
SV4_OE_mean = [SV4_OE_mean_noUnits(1)/1e3, SV4_OE_mean_noUnits(2), ...
               rad2deg(SV4_OE_mean_noUnits(3)), ...
               rad2deg(wrapTo2Pi(SV4_OE_mean_noUnits(4))), ...
               rad2deg(wrapTo2Pi(SV4_OE_mean_noUnits(5))), ...
               rad2deg(wrapTo2Pi(mean2true(SV4_OE_mean_noUnits(6), SV4_OE_mean_noUnits(2))))]';
SV3_OE_mean = [SV3_OE_mean_noUnits(1)/1e3, SV3_OE_mean_noUnits(2), ...
               rad2deg(SV3_OE_mean_noUnits(3)), ...
               rad2deg(wrapTo2Pi(SV3_OE_mean_noUnits(4))), ...
               rad2deg(wrapTo2Pi(SV3_OE_mean_noUnits(5))), ...
               rad2deg(wrapTo2Pi(mean2true(SV3_OE_mean_noUnits(6), SV3_OE_mean_noUnits(2))))]';
chiefOEs(:, 1) = SV4_OE_mean_noUnits';
% get qnsoes
SV3_QNSOE_struct = OE2QNSOE(SV3_OE_mean(1), SV3_OE_mean(2), SV3_OE_mean(3), ...
                            SV3_OE_mean(4), SV3_OE_mean(5), SV3_OE_mean(6));
SV4_QNSOE_struct = OE2QNSOE(SV4_OE_mean(1), SV4_OE_mean(2), SV4_OE_mean(3), ...
                            SV4_OE_mean(4), SV4_OE_mean(5), SV4_OE_mean(6));
SV3_QNSOE = [SV3_QNSOE_struct.a, SV3_QNSOE_struct.ex, SV3_QNSOE_struct.ey, ...
             SV3_QNSOE_struct.inc, SV3_QNSOE_struct.RAAN, wrapTo360(SV3_QNSOE_struct.u)];
SV4_QNSOE = [SV4_QNSOE_struct.a, SV4_QNSOE_struct.ex, SV4_QNSOE_struct.ey, ...
             SV4_QNSOE_struct.inc, SV4_QNSOE_struct.RAAN, wrapTo360(SV4_QNSOE_struct.u)];
x_hat = QNSOE2ROE(SV4_QNSOE, SV3_QNSOE);
roeHist(:, 1) = x_hat';
aroe_init1 = x_hat';
% END x_hat CALCULATION

% cov and noise
P_hat = eye(n)*10; % 10m 
Q = P_hat/1000; 
R = P_hat;

x_plus(:, 1) = x_hat' + sqrtm(P_hat)*randn(n,1);
P_plus(:, :, 1) = P_hat;

% init tracking error
trackingError0 = roe_des' - x_hat';
trackingErrors(:, 1) = trackingError0;

% init ECI
state_SV4 = initial_state_SV4;
state_SV3 = initial_state_SV3;
r_SV4_ECI = [output_SV4_ECI(1,1)'; output_SV4_ECI(1,2)'; output_SV4_ECI(1,3)'];
v_SV4_ECI = [output_SV4_ECI(1,4)'; output_SV4_ECI(1,5)'; output_SV4_ECI(1,6)'];
r_SV3_ECI = [output_SV3_ECI(1,1)'; output_SV3_ECI(1,2)'; output_SV3_ECI(1,3)'];
v_SV3_ECI = [output_SV3_ECI(1,4)'; output_SV3_ECI(1,5)'; output_SV3_ECI(1,6)'];

% TEST mean vs osc
oscSV4 = zeros(n, Nsteps);
oscSV4(:, 1) = [ac, hypot(ex, ey), inc, wc, SV4_OE.Om, SV4_M]';

a_gnss_km = 26560;
e_gnss = 0.01;
i_gnss = 55;
RAAN0_gnss = 0;
w_gnss = 0;
Mref_gnss = 0;
N_gnss = 12;
T_gnss = 3;
F_gnss = 1;

gnss_sats = generate_walker_constellation( ...
    a_gnss_km, e_gnss, i_gnss, RAAN0_gnss, w_gnss, Mref_gnss, ...
    N_gnss, T_gnss, F_gnss);

[r_gnss_eci, v_gnss_eci] = propagate_walker_J2( ...
    gnss_sats, tspan, muearth, reearth, J2earth);

sigma_rho = 3.0; % pseudorange std

% LS init
x0_c_m = r_SV4_ECI * 1e3;
b0_c   = 0;
x0_d_m = r_SV3_ECI * 1e3;
b0_d   = 0;

%% loop and filter
for k = 1:Nsteps-1
    dt = tspan(k+1) - tspan(k);

    %%% calculate the control action %%%
    % [chief_oe; roe; state_ref]
    state = [chiefOEs(:, k); x_plus(:, k)/(a_init_SV4*1e3); roe_des']; % x_plus is dimensionalized

    [statedot, uInp, trackingErr, da_des] = ps6_279d_meanLongControl(0, state); % t is unused so just say 0
    %[statedot, uInp, trackingErr] = ps9_279d_unconstrained(0, state); % UNCONSTRAINED CASE TO TEST
    %uInp = zeros(2,1);

    % set history
    controlStateHist(k, :) = state;
    trackingErrors(:, k) = trackingErr;
    controlInputs(:, k) = uInp;

    % add uInp in RTN to ECI
    r_SV4_ECI = state_SV4(1:3);
    v_SV4_ECI = state_SV4(4:6);
    r_SV3_ECI = state_SV3(1:3);
    v_SV3_ECI = state_SV3(4:6);
    [rot_RTN, w_RTN] = ECI2RTN(r_SV4_ECI, v_SV4_ECI, chiefOEs(1, k), chiefOEs(2, k), ...
        rad2deg(mean2true(deg2rad(chiefOEs(6, k)), chiefOEs(2, k))));
    [r_SV3_RTN, v_SV3_RTN] = rvECI2rvRTN(r_SV3_ECI, v_SV3_ECI, rot_RTN, w_RTN);
    v_SV3_RTN = v_SV3_RTN + [0, uInp(1), uInp(2)]'*dt; % add impulse

    % back to ECI
    [r_SV3_ECI_new, v_SV3_ECI_new] = rvRTN2rvECI(r_SV3_RTN, v_SV3_RTN, rot_RTN, w_RTN);

    %%% simulate numerically one step %%%
    % adjust state with velocity from control action 
    state_SV3 = [r_SV3_ECI_new; v_SV3_ECI_new];
    newTimeSpan = [tspan(k), tspan(k+1)];
    [t_SV4, output_SV4_ECI] = ode45(@orbital_EOM_ECI_J2, newTimeSpan, state_SV4, options);
    [t_SV3, output_SV3_ECI] = ode45(@orbital_EOM_ECI_J2, newTimeSpan, state_SV3, options);
    
    % Convert kth r, v to be used in state estimation
    r_SV4_ECI = [output_SV4_ECI(1,1)'; output_SV4_ECI(1,2)'; output_SV4_ECI(1,3)'];
    v_SV4_ECI = [output_SV4_ECI(1,4)'; output_SV4_ECI(1,5)'; output_SV4_ECI(1,6)'];
    r_SV3_ECI = [output_SV3_ECI(1,1)'; output_SV3_ECI(1,2)'; output_SV3_ECI(1,3)'];
    v_SV3_ECI = [output_SV3_ECI(1,4)'; output_SV3_ECI(1,5)'; output_SV3_ECI(1,6)'];

    % for ground truth take full prop
    state_SV4 = output_SV4_ECI(end, :)';
    state_SV3 = output_SV3_ECI(end, :)';

    %%% get current OEs for truth ROE, STM %%%
    SV3_OE = ECI2OE(r_SV3_ECI, v_SV3_ECI, muearth); % deputy
    SV4_OE = ECI2OE(r_SV4_ECI, v_SV4_ECI, muearth); % chief

    ac = SV4_OE.a;
    wc = SV4_OE.w;
    ex = SV4_OE.e * cosd(wc);
    ey = SV4_OE.e * sind(wc);
    inc = SV4_OE.i;
    SV4_f = wrapTo360(SV4_OE.anom);
    SV3_f = wrapTo360(SV3_OE.anom);
    SV4_M = rad2deg(true2mean(deg2rad(SV4_f), SV4_OE.e));
    SV3_M = rad2deg(true2mean(deg2rad(SV3_f), SV3_OE.e));

    SV4_OE_vec = [SV4_OE.a*1e3, SV4_OE.e, deg2rad(SV4_OE.i), deg2rad(SV4_OE.Om), ...
                  deg2rad(wc), true2mean(deg2rad(SV4_f), SV4_OE.e)];
    SV3_OE_vec = [SV3_OE.a*1e3, SV3_OE.e, deg2rad(SV3_OE.i), deg2rad(SV3_OE.Om), ...
                  deg2rad(SV3_OE.w), true2mean(deg2rad(SV3_f), SV3_OE.e)];

    % get mean from osc
    SV4_OE_mean_noUnits = osc2mean(SV4_OE_vec); % m, rad
    SV3_OE_mean_noUnits = osc2mean(SV3_OE_vec);
    SV4_mean_f = rad2deg(wrapTo2Pi(mean2true(SV4_OE_mean_noUnits(6), SV4_OE_mean_noUnits(2))));
    SV4_OE_mean = [SV4_OE_mean_noUnits(1)/1e3, SV4_OE_mean_noUnits(2), ...
                   rad2deg(SV4_OE_mean_noUnits(3)), ...
                   wrapTo360(rad2deg(SV4_OE_mean_noUnits(4))), ...
                   rad2deg(SV4_OE_mean_noUnits(5)), ...
                   SV4_mean_f]';
    SV3_OE_mean = [SV3_OE_mean_noUnits(1)/1e3, SV3_OE_mean_noUnits(2), ...
                   rad2deg(SV3_OE_mean_noUnits(3)), ...
                   wrapTo360(rad2deg(SV3_OE_mean_noUnits(4))), ...
                   rad2deg(SV3_OE_mean_noUnits(5)), ...
                   rad2deg(wrapTo2Pi(mean2true(SV3_OE_mean_noUnits(6), SV3_OE_mean_noUnits(2))))]';

    % set histories
    chiefOEs(:, k+1) = SV4_OE_mean; 
    oscSV4(:, k+1)   = [ac, hypot(ex, ey), inc, wc, SV4_OE.Om, SV4_M]';

    % get qnsoes for truth ROE
    SV3_QNSOE_struct = OE2QNSOE(SV3_OE_mean(1), SV3_OE_mean(2), SV3_OE_mean(3), ...
                                SV3_OE_mean(4), SV3_OE_mean(5), SV3_OE_mean(6));
    SV4_QNSOE_struct = OE2QNSOE(SV4_OE_mean(1), SV4_OE_mean(2), SV4_OE_mean(3), ...
                                SV4_OE_mean(4), SV4_OE_mean(5), SV4_OE_mean(6));

    SV3_QNSOE = [SV3_QNSOE_struct.a, SV3_QNSOE_struct.ex, SV3_QNSOE_struct.ey, ...
                 SV3_QNSOE_struct.inc, SV3_QNSOE_struct.RAAN, wrapTo360(SV3_QNSOE_struct.u)];
    SV4_QNSOE = [SV4_QNSOE_struct.a, SV4_QNSOE_struct.ex, SV4_QNSOE_struct.ey, ...
                 SV4_QNSOE_struct.inc, SV4_QNSOE_struct.RAAN, wrapTo360(SV4_QNSOE_struct.u)];
    curROE    = (QNSOE2ROE(SV4_QNSOE, SV3_QNSOE))'; % OUTPUT IN METERS!!
    roeHist(:, k+1) = curROE; % set true ROEs

    r_sats_km = squeeze(r_gnss_eci(:, :, k+1)); % 3 x N_gnss
    x_sv_m    = r_sats_km(1,:)' * 1e3;
    y_sv_m    = r_sats_km(2,:)' * 1e3;
    z_sv_m    = r_sats_km(3,:)' * 1e3;
    Nsat_gnss = size(r_sats_km, 2);

    % Chief pseudorange geometry + noise
    ranges_c_km = vecnorm(r_sats_km - r_SV4_ECI, 2, 1);
    ranges_c_m  = ranges_c_km * 1e3;
    corr_pr_c_m = ranges_c_m' + sigma_rho * randn(Nsat_gnss, 1);

    % Deputy pseudoranges
    ranges_d_km = vecnorm(r_sats_km - r_SV3_ECI, 2, 1);
    ranges_d_m  = ranges_d_km * 1e3;
    corr_pr_d_m = ranges_d_m' + sigma_rho * randn(Nsat_gnss, 1);

    max_iter = 10;
    tol_ls   = 1.0;

    % LS for chief
    [x_c_est_m, b_c_est_m, ~] = gnss_localization_ls( ...
        x0_c_m, b0_c, ...
        x_sv_m, y_sv_m, z_sv_m, corr_pr_c_m, ...
        max_iter, tol_ls, false);
    x0_c_m = x_c_est_m;
    b0_c   = b_c_est_m;

    % LS for deputy
    [x_d_est_m, b_d_est_m, ~] = gnss_localization_ls( ...
        x0_d_m, b0_d, ...
        x_sv_m, y_sv_m, z_sv_m, corr_pr_d_m, ...
        max_iter, tol_ls, false);
    x0_d_m = x_d_est_m;
    b0_d   = b_d_est_m;

    % Convert LS results to km
    r_c_est_km = x_c_est_m / 1e3;
    r_d_est_km = x_d_est_m / 1e3;

    v_c_for_ROE_kmps = v_SV4_ECI;
    v_d_for_ROE_kmps = v_SV3_ECI;

    [roe_meas_k, ~, ~] = gnss_epoch_to_roe( ...
        r_c_est_km, v_c_for_ROE_kmps, ...
        r_d_est_km, v_d_for_ROE_kmps, ...
        muearth);

    %%% State estimation %%%
    % STM
    STM = ps8_279d_getJ2STM(dt, SV4_QNSOE_struct.a, SV4_QNSOE_struct.ex, SV4_QNSOE_struct.ey, SV4_QNSOE_struct.inc);
    Bc  = ps6_279d_getReducedControlInput(SV4_QNSOE_struct.a, SV4_QNSOE_struct.ex, SV4_QNSOE_struct.ey, ...
                                          SV4_QNSOE_struct.inc, SV4_OE_mean(5), SV4_mean_f);

    % time update
    x_minus(:, k+1) = STM * x_plus(:, k) + SV4_OE_mean_noUnits(1)*[Bc; zeros(1, 2)]*(uInp.*[1; 1])*dt; 
    P_minus(:, :, k+1) = STM * P_plus(:, :, k) * STM' + Q; 

    %%% measurement update
    y(:, k+1) = roe_meas_k;
    h(:, k+1) = x_minus(:, k+1);

    H = eye(n);
    S_k = H * P_minus(:, :, k+1) * H' + R;
    K(:, :, k+1) = P_minus(:, :, k+1) * H' / S_k;

    x_plus(:, k+1) = x_minus(:, k+1) + K(:, :, k+1) * (y(:, k+1) - h(:, k+1));
    P_plus(:, :, k+1) = (eye(n) - K(:, :, k+1)*H) * P_minus(:, :, k+1);

    % residuals
    residual_pre  = y(:, k+1) - h(:, k+1);
    residual_post = y(:, k+1) - x_plus(:, k+1);
    preResiduals(:, k+1)  = residual_pre;
    postResiduals(:, k+1) = residual_post;
end

%% Plot Results and Eval Performance
stateLabels = {'a \delta a', 'a \delta \lambda', 'a \delta e_x', 'a \delta e_y', 'a \delta i_x', 'a \delta i_y'};
oeLabels    = {'a [km]',    'e [-]',         'i [deg]',       '\Omega [deg]',  '\omega [deg]',  'M [deg]'};

%% Estimation Error
figure; clf;
for i = 1:6
    subplot(3,2,i);
    error_i = x_plus(i,2:end) - roeHist(i,2:end);
    sigma_i = squeeze( sqrt(P_plus(i,i,2:end)) )';  % standard deviation
    
    plot(tspan(2:end)/T_SV4, error_i, 'b', 'LineWidth', 1.5); hold on;
    plot(tspan(2:end)/T_SV4, 3*sigma_i,  'r--');
    plot(tspan(2:end)/T_SV4, -3*sigma_i, 'r--');
    hold off;
    
    xlabel('Time [orbital periods]');
    ylabel(['Error in ', stateLabels{i}]);
    title(['Error & \sigma Bounds: ', stateLabels{i}]);
    grid on;
end

hL = legend({'Estimation Error','+3\sigma','-3\sigma'}, ...
            'FontSize', 8, ...
            'Location','best');
hL.Box = 'off';

%% Residual Scatter Plots
figure('Units','normalized','Position',[0.1 0.1 0.8 0.8]); clf;
for i = 1:6
    subplot(3,2,i);
    scatter(tspan(2:end)/T_SV4, preResiduals(i,2:end), 10, 'r'); hold on;
    scatter(tspan(2:end)/T_SV4, postResiduals(i,2:end), 10, 'g');
    hold off;
    
    xlabel('Time [orbital periods]');
    ylabel(['Residual in ', stateLabels{i}]);
    title(['Measurement Residuals: ', stateLabels{i}]);
    grid on;
end
legend('Pre-fit','Post-fit','Location','best');

%% True vs Estimated States
figure('Units','normalized','Position',[0.1 0.1 0.8 0.8]); clf;
for i = 1:6
    subplot(3,2,i);
    plot(tspan(2:end)/T_SV4, roeHist(i,2:end), 'b', 'LineWidth', 1.5); hold on;
    plot(tspan(2:end)/T_SV4, x_plus(i,2:end),  'r', 'LineWidth', 1.5); hold on;
    plot(tspan(2:end)/T_SV4, aroe_des(i)*ones(size(tspan(2:end))), 'g', 'LineWidth', 1.5);
    hold off;
    
    xlabel('Time [orbital periods]');
    ylabel(stateLabels{i});
    title(['True vs Estimated: ', stateLabels{i}]);
    grid on;
end
legend('True','Estimated','Desired','Location','best');

%% Chief Orbital Elements
figure; clf;
for i = 1:6
    subplot(3,2,i);
    plot(tspan(2:end)/T_SV4, chiefOEs(i,2:end), 'b', 'LineWidth', 1.5); hold on;
    plot(tspan(2:end)/T_SV4, oscSV4(i,2:end),   'r', 'LineWidth', 1.5);
    hold off;
    
    xlabel('Time [orbital periods]');
    ylabel(oeLabels{i});
    title(['Chief OE: ', oeLabels{i}]);
    grid on;
end
legend('Mean','Osc','Location','best');

%% ROEs
figure; clf;
for i = 1:6
    subplot(3,2,i);
    plot(tspan(2:end)/T_SV4, roeHist(i,2:end), 'b', 'LineWidth', 1.5); hold on;
    plot(tspan(2:end)/T_SV4, aroe_des(i)*ones(size(tspan(2:end))), 'g', 'LineWidth', 1.5); hold off;
    
    xlabel('Time [orbital periods]');
    ylabel(stateLabels{i});
    title(['ROE: ', stateLabels{i}]);
    grid on;
end

%% Stats
P_ss = P_plus(:, :, end);
disp('The steady-state covariance from final orbit:');
disp('P =');
disp(P_ss);

sigma_est = sqrt(diag(P_ss))';
fprintf('σ = [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f] m\n', sigma_est);

t_final  = tspan(end);
t_start  = t_final - T_SV4;
idx_final = find(tspan >= t_start);

error_final = x_plus(:, idx_final) - roeHist(:, idx_final);
mu_ss    = mean(error_final, 2)';
sigma_ss = std(error_final, 0, 2)';

fprintf('μss = [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f] m\n', mu_ss);
fprintf('σss = [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f] m\n', sigma_ss);

pre_mean  = mean(preResiduals(:, idx_final), 2)';
pre_std   = std(preResiduals(:, idx_final), 0, 2)';
post_mean = mean(postResiduals(:, idx_final), 2)';
post_std  = std(postResiduals(:, idx_final), 0, 2)';

fprintf('\nTable: Mean and standard deviation statistics of the pre-fit and post-fit residuals over final orbit.\n');
fprintf('%12s %10s %10s %10s %10s %10s %10s\n', '', 'aδa', 'aδλ', 'aδex', 'aδey', 'aδix', 'aδiy');
fprintf('%12s %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n', 'Pre-fit mean', pre_mean);
fprintf('%12s %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n', 'Pre-fit std.', pre_std);
fprintf('%12s %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n', 'Post-fit mean', post_mean);
fprintf('%12s %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n', 'Post-fit std.', post_std);

%% Visualize de, di vectors
figure;
plot(controlStateHist(:, 10)*a_init_SV4*1e3, controlStateHist(:, 9)*a_init_SV4*1e3);
hold on;
plot(aroe_init1(4), aroe_init1(3), Color="r", Marker="square", DisplayName='Initial');
plot(aroe_des(4), aroe_des(3), Color="g", Marker="o", DisplayName='Target');
hold off;
xlabel('dex (m)');
ylabel('dey (m)');
title('Relative Eccentricity Vector');

figure;
plot(controlStateHist(:, 12)*a_init_SV4*1e3, controlStateHist(:, 11)*a_init_SV4*1e3);
hold on;
plot(aroe_init1(6), aroe_init1(5), Color="r", Marker="square", DisplayName='Initial');
plot(aroe_des(6), aroe_des(5), Color="g", Marker="o", DisplayName='Target');
hold off;
xlabel('dix (m)');
ylabel('diy (m)');
title('Relative Inclination Vector');

figure;
plot(controlStateHist(:, 8)*a_init_SV4*1e3, controlStateHist(:, 7)*a_init_SV4*1e3);
hold on;
plot(aroe_init1(2), aroe_init1(1), Color="r", Marker="square", DisplayName='Initial');
plot(aroe_des(2), aroe_des(1), Color="g", Marker="o", DisplayName='Target');
hold off;
xlabel('dlambda (m)');
ylabel('da (m)');
title('Relative sma and lambda');

%% Control Actions and Control Tracking Error
T = T_SV4;
figure;
plot(tspan/T, trackingErrors(1, :)*a_init_SV4*1e3); hold on;
plot(tspan/T, trackingErrors(2, :)*a_init_SV4*1e3);
plot(tspan/T, trackingErrors(3, :)*a_init_SV4*1e3);
plot(tspan/T, trackingErrors(4, :)*a_init_SV4*1e3);
plot(tspan/T, trackingErrors(5, :)*a_init_SV4*1e3);
plot(tspan/T, trackingErrors(6, :)*a_init_SV4*1e3);
hold off;
xlabel('Number of Orbital Periods');
ylabel('Control Tracking Error (m)');
legend('da', 'dlambda', 'dex', 'dey', 'dix', 'diy', Location='southeast');
title('Control Tracking Error vs. Time');
grid on;

figure('Name','Maneuver schedule'), clf, hold on
plot(tspan/T, controlInputs(1,:)*1e3*step_size)
plot(tspan/T, controlInputs(2,:)*1e3*step_size)
xlabel('time  [reconfiguration period]'), ylabel('impulse  [m/s^2]')
legend('Δv_t', 'Δv_n'), grid on, title('Impulse schedule (T, N)')

dVstep = vecnorm(controlInputs,2,1);
cumDV  = cumsum(dVstep, 2);

figure('Name','Cumulative DeltaV'), clf
stairs([0 tspan]/T, [0 cumDV*1e3*step_size], 'LineWidth',1.2)
xlabel('time  [reconfiguration period]'), ylabel('Cumulative DeltaV [m/s]')
grid on, title('Cumulative DeltaV')

figure; clf;
for i = 1:6
    plot(tspan(2:end)/T_SV4, aroe_des(i)*ones(size(tspan(2:end)))-roeHist(i,2:end)); hold on;
end
xlabel('Time [orbital periods]');
ylabel('Error [m]')
grid on;
legend('da', 'dlambda', 'dex', 'dey', 'dix', 'diy', Location='southeast');
title('Ground Truth Control Tracking Error vs. Time');

dV_lb_ip  = 0.5 * n_init_Sv4 * hypot(aroe_des(3)-aroe_init1(3), aroe_des(4)-aroe_init1(4)) / sqrt(1 - e_init_SV4^2);
dV_lb_oop = n_init_Sv4 * (1 - e_init_SV4) * hypot(aroe_des(5)-aroe_init1(5), aroe_des(6)-aroe_init1(6)) / sqrt(1 - e_init_SV4^2);



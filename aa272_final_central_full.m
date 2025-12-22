%% E2E sim with GNSS front-end (chief + deputy)

%% clear everything
clear; close all; clc;

%% DEFINE INTIAL CONDITIONS AND RELEVANT CONSTANTS
tol = 1e-10;

%addpath(genpath(pwd))

addpath('./mean_osc/');
addpath('./helper_functions/');
addpath('./numerical_int/');
addpath('./aa272-final-util/');
addpath('./lyapunov_controller/')

% Earth physical constants from Vallado
muearth = 3.986 * 10^5;
reearth = 6378.1363;
J2earth = 0.0010826266;

%% GNSS FRONT-END CONFIG
c_light = 299792458; % [m/s]
chip_rate = 1.023e6; % [chips/s]
fc = 4.1304e6; % [Hz]
fs = 16e6; % [Hz]
duration_obs = 0.05; % [s]
noise_std = 0.5; % RF noise std

% PLL / DLL loop gains
pll_params.Kp = 50;
pll_params.Ki = 200;
dll_params.Kp_code = 5.0;
dll_params.Ki_code = 1.0;

% LS params
max_iter_ls = 10;
tol_ls = 1.0; 

%% determine orbital elements for chief (SV4)
SV4_init_QNSOE = [6944, -4e-5, 1.6e-3, 99.4, -151.1, -47.9];
a_init_SV4 = SV4_init_QNSOE(1);
ex_init_SV4 = SV4_init_QNSOE(2);
ey_init_SV4 = SV4_init_QNSOE(3);
inclination_init_SV4 = SV4_init_QNSOE(4);
RAAN_init_SV4 = SV4_init_QNSOE(5);
u_init_SV4 = SV4_init_QNSOE(6);

% Determine classical OEs
[a_init_SV4, e_init_SV4, inclination_init_SV4, RAAN_init_SV4, ...
    w_init_SV4, f_init_SV4, M_init_SV4] = QNSOE2OE( ...
        a_init_SV4, ex_init_SV4, ey_init_SV4, ...
        inclination_init_SV4, RAAN_init_SV4, u_init_SV4);

n_init_Sv4 = sqrt(muearth / a_init_SV4^3); % rad/s
T_SV4      = 2*pi*sqrt(a_init_SV4^3/muearth);
nOrbits    = 1/5;
t0         = 0;
tf         = T_SV4*nOrbits;
step_size  = 100;

%% DETERMINE CHIEF INITIAL CONDITIONS AND SETUP NUMERICAL INTEGRATION

% determine initial velocity and position 
[r_init_SV4, v_init_SV4] = OE2ECI( ...
    a_init_SV4, e_init_SV4, inclination_init_SV4, ...
    RAAN_init_SV4, w_init_SV4, f_init_SV4);
initial_state_SV4 = [r_init_SV4; v_init_SV4];

tspan   = 0:step_size:tf;
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-12); % higher accuracy

%% GNSS Walker constellation propagation
a_gnss_km = 26560;
e_gnss = 0.01;
i_gnss_deg = 55;
RAAN0_gnss = 0;
w_gnss_deg = 0;
Mref_gnss = 0;
N_gnss = 32;
T_gnss = 4;
F_gnss = 1;
tic;
sats_gnss = generate_walker_constellation( ...
    a_gnss_km, e_gnss, i_gnss_deg, ...
    RAAN0_gnss, w_gnss_deg, Mref_gnss, ...
    N_gnss, T_gnss, F_gnss);

% Propagate GNSS sats
t_vec_gnss = tspan - 0.125 .* ones(size(tspan)); % rough account for propagation time
[r_gnss_eci, v_gnss_eci] = propagate_walker_J2( ...
    sats_gnss, t_vec_gnss, muearth, reearth, J2earth);

Nsat_gnss = size(r_gnss_eci, 2); 

g2_taps_table = [ ...
    2,6;   3,7;   4,8;   5,9;   1,9;   2,10;  1,8;   2,9;   ... % PRN 1-8
    3,10;  2,3;   3,4;   5,6;   6,7;   7,8;   8,9;   9,10;  ... % PRN 9-16
    1,4;   2,5;   3,6;   4,7;   5,8;   6,9;   1,3;   4,6;   ... % PRN 17-24
    5,7;   6,8;   7,9;   8,10;  1,6;   2,7;   3,8;   4,9    ... % PRN 25-32
];

% PRN codes for each GNSS sat
prn_pm_list = cell(1, Nsat_gnss);
for i = 1:Nsat_gnss
    prn_id = mod(i-1, 32) + 1;
    
    idx1 = g2_taps_table(prn_id, 1);
    idx2 = g2_taps_table(prn_id, 2);
    
    % Generate the C/A code
    prn_bits       = create_prn(idx1, idx2, 1023);
    prn_pm_list{i} = map_binary(prn_bits);
end

%% DEFINE DEPUTY INITIAL CONDITIONS (SV3)
SV3_init_ROE = [-1, -79328, 42, 452, 36, 827];

% set exact a component to 0 for this problem set
SV3_init_ROE(1) = 0;

% assuming init in degrees here
SV3_init_QNSOE = ROE2QNSOE( ...
    a_init_SV4, ex_init_SV4, ey_init_SV4, ...
    inclination_init_SV4, RAAN_init_SV4, M_init_SV4 + w_init_SV4, ...
    SV3_init_ROE(1), SV3_init_ROE(2), SV3_init_ROE(3), ...
    SV3_init_ROE(4), SV3_init_ROE(5), SV3_init_ROE(6));

[a_init_SV3, e_init_SV3, inclination_init_SV3, RAAN_init_SV3, ...
    w_init_SV3, f_init_SV3, M_init_SV3] = QNSOE2OE( ...
        SV3_init_QNSOE(1), SV3_init_QNSOE(2), SV3_init_QNSOE(3), ...
        SV3_init_QNSOE(4), SV3_init_QNSOE(5), SV3_init_QNSOE(6));

[r_init_SV3, v_init_SV3] = OE2ECI( ...
    a_init_SV3, e_init_SV3, inclination_init_SV3, ...
    RAAN_init_SV3, w_init_SV3, f_init_SV3);

%% TRANSFORM DEPUTY TO CHIEF RTN

% define chief RTN frame
[rot_RTN, w_RTN] = ECI2RTN(r_init_SV4, v_init_SV4, a_init_SV4, e_init_SV4, f_init_SV4);

% take SV3 to chief RTN frame
[rho_init_SV3_RTN, rho_dot_init_SV3_RTN] = rvECI2rvRTN( ...
    r_init_SV3 - r_init_SV4, v_init_SV3 - v_init_SV4, rot_RTN, w_RTN);
[r_init_Sv4_RTN, v_init_Sv4_RTN] = rvECI2rvRTN( ...
    r_init_SV4, v_init_SV4, rot_RTN, w_RTN);

%% RUN SIMULATION FOR SV3 AND SV4 IN THE ABSOLUTE FRAME

initial_state_SV4 = [r_init_SV4; v_init_SV4];
initial_state_SV3 = [r_init_SV3; v_init_SV3];

%% RUN NUMERICAL INTEGRATION (truth ECI orbits)
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
    [rot_RTN, w_RTN] = ECI2RTN(r_SV4_ECI(:, i), v_SV4_ECI(:, i), a_init_SV4, e_init_SV4, f_init_SV4);
    [r_SV4_RTN(:, i), v_SV4_RTN(:, i)] = rvECI2rvRTN(r_SV4_ECI(:, i), v_SV4_ECI(:, i), rot_RTN, w_RTN);
    [r_SV3_RTN(:, i), v_SV3_RTN(:, i)] = rvECI2rvRTN(r_SV3_ECI(:, i), v_SV3_ECI(:, i), rot_RTN, w_RTN);
end

% Find deputy with respect to chief orbit
rho_SV3_SV4_RTN = r_SV3_RTN - r_SV4_RTN;

%% PLOT RESULTS OF ABSOLUTE FRAME SIMULATION IN RTN FRAME (sanity check)
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

%% FODE HEALTH CHECK COMPLETE %%

%% Begin Simulation Setup

% Desired ROE state 
aroe_des = [0, 7920, -297, -315, 16, 1144];
roe_des  = aroe_des/(a_init_SV4*1e3);

%% Initializations and Setup

n = 6; % state dim

Nsteps   = numel(tspan);
x_plus   = zeros(n, Nsteps);
x_minus  = zeros(n, Nsteps);
P_plus   = zeros(n, n, Nsteps);
P_minus  = zeros(n, n, Nsteps);
K        = zeros(n, n, Nsteps);
h        = zeros(n, Nsteps);
y        = zeros(n, Nsteps);
roeHist  = zeros(n, Nsteps);
preResiduals  = zeros(n, Nsteps);
postResiduals = zeros(n, Nsteps);
chiefOEs = zeros(n, Nsteps);

% control inits
trackingErrors   = zeros(n, Nsteps);
controlInputs    = zeros(2, Nsteps);
controlStateHist = zeros(Nsteps-1, 18);

% seed for now
seed = 42;
rng(seed);

%% Initial ROE estimate x_hat

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

SV4_OE_vec = [SV4_OE.a*1e3, SV4_OE.e, deg2rad(SV4_OE.i), ...
              deg2rad(SV4_OE.Om), deg2rad(wc), ...
              true2mean(deg2rad(SV4_f), SV4_OE.e)];
SV3_OE_vec = [SV3_OE.a*1e3, SV3_OE.e, deg2rad(SV3_OE.i), ...
              deg2rad(SV3_OE.Om), deg2rad(wrapTo360(SV3_OE.w)), ...
              true2mean(deg2rad(SV3_f), SV3_OE.e)];

SV4_OE_mean_noUnits = osc2mean(SV4_OE_vec); % meters and rad
SV3_OE_mean_noUnits = osc2mean(SV3_OE_vec);

SV4_OE_mean = [SV4_OE_mean_noUnits(1)/1e3, ...
               SV4_OE_mean_noUnits(2), ...
               rad2deg(SV4_OE_mean_noUnits(3)), ...
               rad2deg(wrapTo2Pi(SV4_OE_mean_noUnits(4))), ...
               rad2deg(wrapTo2Pi(SV4_OE_mean_noUnits(5))), ...
               rad2deg(wrapTo2Pi(mean2true(SV4_OE_mean_noUnits(6), SV4_OE_mean_noUnits(2))))]';
SV3_OE_mean = [SV3_OE_mean_noUnits(1)/1e3, ...
               SV3_OE_mean_noUnits(2), ...
               rad2deg(SV3_OE_mean_noUnits(3)), ...
               rad2deg(wrapTo2Pi(SV3_OE_mean_noUnits(4))), ...
               rad2deg(wrapTo2Pi(SV3_OE_mean_noUnits(5))), ...
               rad2deg(wrapTo2Pi(mean2true(SV3_OE_mean_noUnits(6), SV3_OE_mean_noUnits(2))))]';

chiefOEs(:, 1) = SV4_OE_mean_noUnits';

SV3_QNSOE_struct = OE2QNSOE(SV3_OE_mean(1), SV3_OE_mean(2), SV3_OE_mean(3), ...
                            SV3_OE_mean(4), SV3_OE_mean(5), SV3_OE_mean(6));
SV4_QNSOE_struct = OE2QNSOE(SV4_OE_mean(1), SV4_OE_mean(2), SV4_OE_mean(3), ...
                            SV4_OE_mean(4), SV4_OE_mean(5), SV4_OE_mean(6));
SV3_QNSOE = [SV3_QNSOE_struct.a, SV3_QNSOE_struct.ex, SV3_QNSOE_struct.ey, ...
             SV3_QNSOE_struct.inc, SV3_QNSOE_struct.RAAN, wrapTo360(SV3_QNSOE_struct.u)];
SV4_QNSOE = [SV4_QNSOE_struct.a, SV4_QNSOE_struct.ex, SV4_QNSOE_struct.ey, ...
             SV4_QNSOE_struct.inc, SV4_QNSOE_struct.RAAN, wrapTo360(SV4_QNSOE_struct.u)];

x_hat        = QNSOE2ROE(SV4_QNSOE, SV3_QNSOE);  % a*delta form, [m]
roeHist(:,1) = x_hat';
aroe_init1   = x_hat';

% cov and noise
P_hat = eye(n)*10; % 10 m
Q     = P_hat/100;
R     = P_hat*0.1;

% initialize EKF state and covariance
x_plus(:, 1)     = x_hat' + sqrtm(P_hat)*randn(n,1);
P_plus(:, :, 1)  = P_hat;

% init control tracking error
trackingError0     = roe_des' - x_hat';
trackingErrors(:,1)= trackingError0;

% init ECI states
state_SV4 = initial_state_SV4;
state_SV3 = initial_state_SV3;
r_SV4_ECI = [output_SV4_ECI(1,1)'; output_SV4_ECI(1,2)'; output_SV4_ECI(1,3)'];
v_SV4_ECI = [output_SV4_ECI(1,4)'; output_SV4_ECI(1,5)'; output_SV4_ECI(1,6)'];
r_SV3_ECI = [output_SV3_ECI(1,1)'; output_SV3_ECI(1,2)'; output_SV3_ECI(1,3)'];
v_SV3_ECI = [output_SV3_ECI(1,4)'; output_SV3_ECI(1,5)'; output_SV3_ECI(1,6)'];

% LS init
x0_c_m = r_SV4_ECI * 1e3;
x0_d_m = r_SV3_ECI * 1e3;
b0_c   = 0;
b0_d   = 0;
pdop_c_hist = zeros(Nsteps, 1);
pdop_d_hist = zeros(Nsteps, 1);

% fallback out of GNSS counter
fallback_counter = 0;

% TEST mean vs osc
oscSV4 = zeros(n, Nsteps);
oscSV4(:, 1) = [ac, hypot(ex, ey), inc, wc, SV4_OE.Om, SV4_M]';

%% loop and filter
for k = 1:Nsteps-1
    dt = tspan(k+1) - tspan(k);

    %%% calculate the control action %%%
    % [chief_oe; roe; state_ref]
    state = [chiefOEs(:, k); x_plus(:, k)/(a_init_SV4*1e3); roe_des']; % x_plus dimensionalized

    [statedot, uInp, trackingErr, da_des] = ps6_279d_meanLongControl(0, state); % t unused
    % [statedot, uInp, trackingErr] = ps9_279d_unconstrained(0, state);
    % uInp = zeros(2,1);

    controlStateHist(k, :) = state;
    trackingErrors(:, k)   = trackingErr;
    controlInputs(:, k)    = uInp;

    % add uInp in RTN to ECI
    r_SV4_ECI = state_SV4(1:3);
    v_SV4_ECI = state_SV4(4:6);
    r_SV3_ECI = state_SV3(1:3);
    v_SV3_ECI = state_SV3(4:6);

    [rot_RTN, w_RTN] = ECI2RTN( ...
        r_SV4_ECI, v_SV4_ECI, ...
        chiefOEs(1, k), chiefOEs(2, k), ...
        rad2deg(mean2true(deg2rad(chiefOEs(6, k)), chiefOEs(2, k))));

    [r_SV3_RTN, v_SV3_RTN] = rvECI2rvRTN(r_SV3_ECI, v_SV3_ECI, rot_RTN, w_RTN);
    v_SV3_RTN = v_SV3_RTN + [0, uInp(1), uInp(2)]'*dt; % add impulse

    % back to ECI for SV3
    [r_SV3_ECI_new, v_SV3_ECI_new] = rvRTN2rvECI(r_SV3_RTN, v_SV3_RTN, rot_RTN, w_RTN);
    state_SV3 = [r_SV3_ECI_new; v_SV3_ECI_new];

    %%% simulate numerically one step truth orbits %%%
    newTimeSpan = [tspan(k), tspan(k+1)];
    [~, output_SV4_ECI] = ode45(@orbital_EOM_ECI_J2, newTimeSpan, state_SV4, options);
    [~, output_SV3_ECI] = ode45(@orbital_EOM_ECI_J2, newTimeSpan, state_SV3, options);

    r_SV4_ECI = [output_SV4_ECI(1,1)'; output_SV4_ECI(1,2)'; output_SV4_ECI(1,3)'];
    v_SV4_ECI = [output_SV4_ECI(1,4)'; output_SV4_ECI(1,5)'; output_SV4_ECI(1,6)'];
    r_SV3_ECI = [output_SV3_ECI(1,1)'; output_SV3_ECI(1,2)'; output_SV3_ECI(1,3)'];
    v_SV3_ECI = [output_SV3_ECI(1,4)'; output_SV3_ECI(1,5)'; output_SV3_ECI(1,6)'];

    state_SV4 = output_SV4_ECI(end, :)';
    state_SV3 = output_SV3_ECI(end, :)';

   %% GNSS FRONT-END
    % GNSS satellite states
    r_sat_eci_km   = squeeze(r_gnss_eci(:, :, k+1)); % 3 x Nsat_gnss
    v_sat_eci_kmps = squeeze(v_gnss_eci(:, :, k+1));
    
    % Ground truth for "perfect acquisition"
    x0_c_m = r_SV4_ECI * 1e3;
    x0_d_m = r_SV3_ECI * 1e3;
    
    % ================= CHIEF (SV4) =================
    [true_delays_c, true_dopplers_c, ranges_c_m] = compute_gnss_delays_dopplers( ...
        r_sat_eci_km, v_sat_eci_kmps, ...
        r_SV4_ECI, v_SV4_ECI*0, ...  % no dop rn
        fc, chip_rate);
    
    gnss_c = track_multisatellite( ...
        prn_pm_list, chip_rate, fc, fs, duration_obs, ...
        true_delays_c, true_dopplers_c, noise_std, ...
        pll_params, dll_params, true); 
    
    corr_pr_c_m = gnss_c.pseudorange_m(end, :).';
    
    % Code solution
    [x_c_est_code, b_c_est_m, ~, pdop_c] = gnss_localization_ls( ...
        x0_c_m, b0_c, ...
        r_sat_eci_km(1,:)'*1e3, ...
        r_sat_eci_km(2,:)'*1e3, ...
        r_sat_eci_km(3,:)'*1e3, ...
        corr_pr_c_m, ...
        max_iter_ls, tol_ls, false);
    
    gnss_c.wavelength_m = 0.19029; % set before but dbl check
    % wavelength = gnss_c.wavelength_m;
    
    x_c_est_m = x_c_est_code;
    b0_c = b_c_est_m;
    
    % ================= DEPUTY (SV3) =================
    [true_delays_d, true_dopplers_d, ranges_d_m] = compute_gnss_delays_dopplers( ...
        r_sat_eci_km, v_sat_eci_kmps, ...
        r_SV3_ECI, v_SV3_ECI*0, ...
        fc, chip_rate);
    
    gnss_d = track_multisatellite( ...
        prn_pm_list, chip_rate, fc, fs, duration_obs, ...
        true_delays_d, true_dopplers_d, noise_std, ...
        pll_params, dll_params, true);
    
    corr_pr_d_m = gnss_d.pseudorange_m(end, :).';
    gnss_c.wavelength_m = 0.19029;
    
    % Code solution
    [x_d_est_code, b_d_est_m, ~, pdop_d] = gnss_localization_ls( ...
        x0_d_m, b0_d, ...
        r_sat_eci_km(1,:)'*1e3, ...
        r_sat_eci_km(2,:)'*1e3, ...
        r_sat_eci_km(3,:)'*1e3, ...
        corr_pr_d_m, ...
        max_iter_ls, tol_ls, false);
    
    % Differential solution
    delta_r_m = solve_baseline_dd_integer(r_SV4_ECI*1e3, r_sat_eci_km*1e3, gnss_c, gnss_d);
    
    % not proper but ease of conversion
    x_d_est_m = (r_SV4_ECI * 1e3) + delta_r_m; % x_d_est_code; %
    % x0_d_m = x_d_est_m;u
    b0_d   = b_d_est_m;
    
    % Convert LS results to km for ROE calculation
    r_c_est_km = x_c_est_m / 1e3;
    r_d_est_km = x_d_est_m / 1e3;
    
    v_c_for_ROE_kmps = v_SV4_ECI;
    v_d_for_ROE_kmps = v_SV3_ECI;
    
    SV3_OE = ECI2OE(r_SV3_ECI, v_SV3_ECI, muearth);
    SV4_OE = ECI2OE(r_SV4_ECI, v_SV4_ECI, muearth);

    ac  = SV4_OE.a;
    wc  = SV4_OE.w;
    ex  = SV4_OE.e * cosd(wc);
    ey  = SV4_OE.e * sind(wc);
    inc = SV4_OE.i;
    SV4_f = wrapTo360(SV4_OE.anom);
    SV3_f = wrapTo360(SV3_OE.anom);
    SV4_M = rad2deg(true2mean(deg2rad(SV4_f), SV4_OE.e));
    SV3_M = rad2deg(true2mean(deg2rad(SV3_f), SV3_OE.e));

    SV4_OE_vec = [SV4_OE.a*1e3, SV4_OE.e, deg2rad(SV4_OE.i), ...
                  deg2rad(SV4_OE.Om), deg2rad(wc), ...
                  true2mean(deg2rad(SV4_f), SV4_OE.e)];
    SV3_OE_vec = [SV3_OE.a*1e3, SV3_OE.e, deg2rad(SV3_OE.i), ...
                  deg2rad(SV3_OE.Om), deg2rad(SV3_OE.w), ...
                  true2mean(deg2rad(SV3_f), SV3_OE.e)];

    SV4_OE_mean_noUnits = osc2mean(SV4_OE_vec); % meters and rad
    SV3_OE_mean_noUnits = osc2mean(SV3_OE_vec);

    SV4_mean_f = rad2deg(wrapTo2Pi(mean2true(SV4_OE_mean_noUnits(6), SV4_OE_mean_noUnits(2))));
    SV4_OE_mean = [SV4_OE_mean_noUnits(1)/1e3, ...
                   SV4_OE_mean_noUnits(2), ...
                   rad2deg(SV4_OE_mean_noUnits(3)), ...
                   wrapTo360(rad2deg(SV4_OE_mean_noUnits(4))), ...
                   rad2deg(SV4_OE_mean_noUnits(5)), ...
                   SV4_mean_f]';
    SV3_OE_mean = [SV3_OE_mean_noUnits(1)/1e3, ...
                   SV3_OE_mean_noUnits(2), ...
                   rad2deg(SV3_OE_mean_noUnits(3)), ...
                   wrapTo360(rad2deg(SV3_OE_mean_noUnits(4))), ...
                   rad2deg(SV3_OE_mean_noUnits(5)), ...
                   rad2deg(wrapTo2Pi(mean2true(SV3_OE_mean_noUnits(6), SV3_OE_mean_noUnits(2))))]';

    chiefOEs(:, k+1) = SV4_OE_mean;
    oscSV4(:, k+1)   = [ac, hypot(ex, ey), inc, wc, SV4_OE.Om, SV4_M]';

    SV3_QNSOE_struct = OE2QNSOE(SV3_OE_mean(1), SV3_OE_mean(2), SV3_OE_mean(3), ...
                                SV3_OE_mean(4), SV3_OE_mean(5), SV3_OE_mean(6));
    SV4_QNSOE_struct = OE2QNSOE(SV4_OE_mean(1), SV4_OE_mean(2), SV4_OE_mean(3), ...
                                SV4_OE_mean(4), SV4_OE_mean(5), SV4_OE_mean(6));

    SV3_QNSOE = [SV3_QNSOE_struct.a, SV3_QNSOE_struct.ex, SV3_QNSOE_struct.ey, ...
                 SV3_QNSOE_struct.inc, SV3_QNSOE_struct.RAAN, wrapTo360(SV3_QNSOE_struct.u)];
    SV4_QNSOE = [SV4_QNSOE_struct.a, SV4_QNSOE_struct.ex, SV4_QNSOE_struct.ey, ...
                 SV4_QNSOE_struct.inc, SV4_QNSOE_struct.RAAN, wrapTo360(SV4_QNSOE_struct.u)];

    curROE          = (QNSOE2ROE(SV4_QNSOE, SV3_QNSOE))'; % [m]
    roeHist(:,k+1)  = curROE;

    %%% State estimation EKF %%%
    STM = ps8_279d_getJ2STM(dt, SV4_QNSOE_struct.a, SV4_QNSOE_struct.ex, ...
                            SV4_QNSOE_struct.ey, SV4_QNSOE_struct.inc);
    Bc  = ps6_279d_getReducedControlInput( ...
               SV4_QNSOE_struct.a, SV4_QNSOE_struct.ex, SV4_QNSOE_struct.ey, ...
               SV4_QNSOE_struct.inc, SV4_OE_mean(5), SV4_mean_f);

    % time update
    x_minus(:, k+1)    = STM * x_plus(:, k) + ...
                         SV4_OE_mean_noUnits(1)*[Bc; zeros(1, 2)]*(uInp.*[1; 1])*dt;
    P_minus(:, :, k+1) = STM * P_plus(:, :, k) * STM' + Q;

    % CHECK MEASUREMENT OK
    if norm(r_c_est_km) > 1e6
        roe_meas_k = curROE + sqrtm(R)*randn(n,1);
        fallback_counter = fallback_counter + 1;
        if fallback_counter == 1
            disp(k);
        end
    else
        % GNSS-based ROE measurement (meters)
        [roe_meas_k, ~, ~] = gnss_epoch_to_roe( ...
            r_c_est_km, v_c_for_ROE_kmps, ...
            r_d_est_km, v_d_for_ROE_kmps, ...
            muearth);

        % log PDOP
        pdop_d_hist(k) = pdop_d;
        pdop_c_hist(k) = pdop_c;
    end

    % measurement update
    y(:, k+1) = roe_meas_k;
    h(:, k+1) = x_minus(:, k+1);
    H         = eye(n);

    S           = H * P_minus(:, :, k+1)*H' + R;
    K(:, :,k+1) = P_minus(:, :, k+1) * H' / S;
    x_plus(:,k+1) = x_minus(:,k+1) + K(:, :,k+1) * (y(:,k+1) - h(:,k+1));
    P_plus(:, :,k+1) = (eye(n)-K(:, :,k+1)*H)*P_minus(:, :,k+1)* ...
                       (eye(n)-K(:, :,k+1)*H)' + K(:, :,k+1)*R*K(:, :,k+1)';

    residual_pre  = y(:, k+1) - h(:, k+1);
    residual_post = y(:, k+1) - x_plus(:, k+1);

    preResiduals(:, k+1)  = residual_pre;
    postResiduals(:,k+1)  = residual_post;
end

%% Plot Results and Eval Performance
stateLabels = {'a \delta a', 'a \delta \lambda', 'a \delta e_x', ...
               'a \delta e_y', 'a \delta i_x', 'a \delta i_y'};
oeLabels    = {'a [km]', 'e [-]', 'i [deg]', '\Omega [deg]', '\omega [deg]', 'M [deg]'};

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
            'FontSize', 8, 'Location','best');
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
legend('True','Estimated','Location','best');

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
legend('True','Osc SV4','Location','best');

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
disp('P ='); disp(P_ss);

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

fprintf('\nTable: Mean and standard deviation of pre/post-fit residuals over final orbit.\n');
fprintf('%12s %10s %10s %10s %10s %10s %10s\n', '', 'aδa', 'aδλ', 'aδex', 'aδey', 'aδix', 'aδiy');
fprintf('%12s %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n', 'Pre-fit mean', pre_mean);
fprintf('%12s %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n', 'Pre-fit std.', pre_std);
fprintf('%12s %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n', 'Post-fit mean', post_mean);
fprintf('%12s %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n', 'Post-fit std.', post_std);

%% Visualize de, di vectors
figure;
plot(controlStateHist(:, 10)*a_init_SV4*1e3, controlStateHist(:, 9)*a_init_SV4*1e3);
hold on;
plot(aroe_init1(4), aroe_init1(3), 'rs', 'DisplayName','Initial');
plot(aroe_des(4), aroe_des(3), 'go', 'DisplayName','Target');
hold off;
xlabel('dex (m)'); ylabel('dey (m)');
title('Relative Eccentricity Vector'); grid on;

figure;
plot(controlStateHist(:, 12)*a_init_SV4*1e3, controlStateHist(:, 11)*a_init_SV4*1e3);
hold on;
plot(aroe_init1(6), aroe_init1(5), 'rs', 'DisplayName','Initial');
plot(aroe_des(6), aroe_des(5), 'go', 'DisplayName','Target');
hold off;
xlabel('dix (m)'); ylabel('diy (m)');
title('Relative Inclination Vector'); grid on;

figure;
plot(controlStateHist(:, 8)*a_init_SV4*1e3, controlStateHist(:, 7)*a_init_SV4*1e3);
hold on;
plot(aroe_init1(2), aroe_init1(1), 'rs', 'DisplayName','Initial');
plot(aroe_des(2), aroe_des(1), 'go', 'DisplayName','Target');
hold off;
xlabel('dlambda (m)'); ylabel('da (m)');
title('Relative sma and lambda'); grid on;

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
legend('da', 'dlambda', 'dex', 'dey', 'dix', 'diy', 'Location','southeast');
title('Control Tracking Error vs. Time');
grid on;

figure('Name','Maneuver schedule'); clf; hold on
plot(tspan/T, controlInputs(1,:)*1e3*step_size);
plot(tspan/T, controlInputs(2,:)*1e3*step_size);
xlabel('time  [reconfiguration period]');
ylabel('impulse  [m/s^2]');
legend('\Delta v_t', '\Delta v_n');
grid on; title('Impulse schedule (T, N)');

dVstep = vecnorm(controlInputs,2,1);
cumDV  = cumsum(dVstep, 2);
figure('Name','Cumulative DeltaV'); clf;
stairs([0 tspan]/T, [0 cumDV*1e3*step_size], 'LineWidth',1.2);
xlabel('time  [reconfiguration period]');
ylabel('Cumulative DeltaV [m/s]');
grid on; title('Cumulative DeltaV');

figure; clf;
for i = 1:6
    plot(tspan(2:end)/T_SV4, aroe_des(i)*ones(size(tspan(2:end)))-roeHist(i,2:end)); hold on;
end
xlabel('Time [orbital periods]');
ylabel('Error [m]');
grid on;
legend('da', 'dlambda', 'dex', 'dey', 'dix', 'diy', 'Location','southeast');
title('Ground Truth Control Tracking Error vs. Time');

%% PDOP Plots
figure;
plot(tspan/T, pdop_c_hist);
xlabel('time  [reconfiguration period]');
ylabel('Chief PDOP');
grid on; title('Chief PDOP Over Time');

figure;
plot(tspan/T, pdop_d_hist);
xlabel('time  [reconfiguration period]');
ylabel('Deputy PDOP');
grid on; title('Deputy PDOP Over Time');

%% Final Printouts
% Lower bounds on delta-V (for context)
dV_lb_ip  = 0.5 * n_init_Sv4 * hypot(aroe_des(3)-aroe_init1(3), aroe_des(4)-aroe_init1(4)) / sqrt(1 - e_init_SV4^2);
dV_lb_oop = n_init_Sv4 * (1 - e_init_SV4) * hypot(aroe_des(5)-aroe_init1(5), aroe_des(6)-aroe_init1(6)) / sqrt(1 - e_init_SV4^2);
fprintf('\nDelta-V lower bounds (ip / oop) [km/s]: %.6f / %.6f\n', dV_lb_ip, dV_lb_oop);

fprintf('Fallback Percentage of: %.2f\n', fallback_counter/(Nsteps-1));
end_t = toc; disp(end_t);

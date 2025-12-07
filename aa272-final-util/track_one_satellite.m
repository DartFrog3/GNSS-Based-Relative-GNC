%% FULL E2E Function for Single Sat -- GNSS chain from RF signal → baseband → correlators / DLL / PLL

function results = track_one_satellite( ...
    prn_pm, chip_rate, fc, fs, duration, ...
    true_delay_chips, true_doppler_hz, noise_std, ...
    pll_params, dll_params, use_sdcp)

    if nargin < 11
        use_sdcp = false;
    end

    % Complete single sat tracking loop

    % get transmitted signal
    phi0 = 0;
    [t, s_rf, ~] = build_gnss_signal( ...
        prn_pm, chip_rate, fc, fs, ...
        true_delay_chips, true_doppler_hz, ...
        phi0, duration, noise_std);

    Nsamp = numel(s_rf);
    dt = 1/fs; % by sample rate

    % baseband mixing
    fc_hat0 = fc + true_doppler_hz;
    [I, Q] = mix_to_baseband(s_rf, fc_hat0, fs, 0);

    % Setup for tracking
    Tcoh = 1e-3;
    Ns_coh = round(Tcoh * fs);
    Nblock = floor(Nsamp / Ns_coh);
    t_blocks = ((0:Nblock-1) * Tcoh);
    c_light = 299792458;
    wavelength_m = c_light / fc;

    Nchips = numel(prn_pm);
    code_phase_hat = true_delay_chips;
    code_int = 0.0;
    code_freq_hat = chip_rate;
    pll = pll_init(fc_hat0);

    % prellocate
    code_phase_hist = zeros(1, Nblock);
    pseudorange_hist = zeros(1, Nblock);
    dll_disc_hist = zeros(1, Nblock);
    carrier_phase_hist = zeros(1, Nblock);
    doppler_hist = zeros(1, Nblock);
    carrier_phase_cycles_hist = zeros(1, Nblock);
    phase_unwrapped = 0;
    phase_prev_wrapped = 0;

    % Actually track
    for k = 1:Nblock
        idx_start = (k-1)*Ns_coh + 1;
        idx_end = k*Ns_coh;

        I_blk = I(idx_start:idx_end);
        Q_blk = Q(idx_start:idx_end);
        t_blk = t(idx_start:idx_end);

        % Create local replicas
        t_code = t_blk * chip_rate - code_phase_hat;
        t_code_P = t_code;
        t_code_E = t_code + 0.5;
        t_code_L = t_code - 0.5;

        idxP = mod(floor(t_code_P), Nchips) + 1;
        idxE = mod(floor(t_code_E), Nchips) + 1;
        idxL = mod(floor(t_code_L), Nchips) + 1;

        code_P = prn_pm(idxP);
        code_E = prn_pm(idxE);
        code_L = prn_pm(idxL);

        % Correlators
        IP = sum(I_blk .* code_P);
        QP = sum(Q_blk .* code_P);
        IE = sum(I_blk .* code_E);
        QE = sum(Q_blk .* code_E);
        IL = sum(I_blk .* code_L);
        QL = sum(Q_blk .* code_L);

        % from lec11 slide9:
        Power_E = IE^2 + QE^2; 
        Power_L = IL^2 + QL^2;

        % early minus late Power Discriminator
        D_dll = (Power_E - Power_L) / (Power_E + Power_L + eps);
        
        dll_disc_hist(k) = D_dll;

        code_int = code_int + dll_params.Ki_code * D_dll * Tcoh;
        code_freq_hat = chip_rate + dll_params.Kp_code * D_dll + code_int;
        code_phase_hat = code_phase_hat + code_freq_hat * Tcoh;

        % Store estimated code phase and pseudorange
        code_phase_hist(k) = code_phase_hat;
        tau_hat = code_phase_hat / chip_rate;
        pseudorange_hist(k) = c_light * tau_hat;

        % PLL update
        [pll, phase_err, ~, ~] = pll_update(pll, IP, QP, Tcoh, pll_params);

        carrier_phase_hist(k) = pll.phase;
        doppler_hist(k) = pll.freq - pll.nominal_freq; % dopp

        % track cycles for relative solving for cdgnss
        % still slip vulnerable since not absolute?
        if use_sdcp
            phase_wrapped = pll.phase;
            phase_diff = phase_wrapped - phase_prev_wrapped;
            
            if phase_diff > pi
                phase_diff = phase_diff - 2*pi;
            elseif phase_diff < -pi
                phase_diff = phase_diff + 2*pi;
            end
            
            phase_unwrapped = phase_unwrapped + phase_diff;
            phase_prev_wrapped = phase_wrapped;
            carrier_phase_cycles_hist(k) = phase_unwrapped / (2*pi);
        else
            carrier_phase_cycles_hist(k) = NaN;
        end
    end

    % Pack out
    results.t_blocks_s        = t_blocks;
    results.code_phase_chips  = code_phase_hist;
    results.pseudorange_m     = pseudorange_hist;
    results.dll_disc          = dll_disc_hist;
    results.carrier_phase_rad = carrier_phase_hist;
    results.doppler_est_hz    = doppler_hist;
    results.true_delay_chips  = true_delay_chips;
    results.true_doppler_hz   = true_doppler_hz;
    results.carrier_phase_cycles = carrier_phase_cycles_hist;
    results.wavelength_m         = wavelength_m;
end


%% PLL DOUBLE CHECK THIS W BOOK
function pll = pll_init(nominal_freq_hz)
    % PLL struct setup

    pll.phase = 0;
    pll.freq = nominal_freq_hz;
    pll.integrator = 0;
    pll.nominal_freq = nominal_freq_hz;
end

function [pll, phase_error, nco_cos, nco_sin] = pll_update(pll, I, Q, dt, params)
    % Costas discriminator (sect. 12.3.2 of book)
    phase_error = I*Q;

    %PI update
    pll.integrator = pll.integrator + params.Ki * phase_error * dt;
    freq_correction = params.Kp * phase_error + pll.integrator;

    % Updated NCO frequency [Hz]
    pll.freq = pll.nominal_freq + freq_correction;
    pll.phase = pll.phase + 2*pi * pll.freq * dt;   % [rad]

    pll.phase = atan2(sin(pll.phase), cos(pll.phase));
    nco_cos = cos(pll.phase);
    nco_sin = sin(pll.phase);
end



%% SEND SIGNAL helper
function [t, s_rf, c_seq] = build_gnss_signal(prn_pm, chip_rate, fc, fs, ...
                                              delay_chips, doppler_hz, ...
                                              phi0, duration, noise_std)
    % Build simulated RF signal

    % noised signal option
    if nargin < 9
        noise_std = 0;
    end

    Nsamp = round(duration * fs);
    t = (0:Nsamp-1) / fs;  % [s]

    Nchips = numel(prn_pm);
    t_code = t * chip_rate - delay_chips;
    t_code_wrapped = mod(t_code, Nchips);

    chip_idx = floor(t_code_wrapped);
    chip_idx = chip_idx + 1;
    c_seq = prn_pm(chip_idx);

    % include doppler
    f_eff = fc + doppler_hz;
    phase = 2*pi*f_eff.*t + phi0;
    carrier = cos(phase);

    % element wise
    s_rf_clean = c_seq .* carrier;

    if noise_std > 0
        s_rf = s_rf_clean + noise_std * randn(size(s_rf_clean));
    else
        s_rf = s_rf_clean;
    end
end


%% To I-Q domain as already digital ofc
function [I, Q] = mix_to_baseband(s_rf, fc_hat, fs, phi0)
    N = numel(s_rf);
    t = (0:N-1) / fs;

    lo_cos = cos(2*pi*fc_hat .* t + phi0);
    lo_sin = -sin(2*pi*fc_hat .* t + phi0);  % conj?

    I = s_rf .* lo_cos;
    Q = s_rf .* lo_sin;
end



%% Multi Sat E2E function
function gnss = track_multisatellite( ...
    prn_pm_list, chip_rate, fc, fs, duration, ...
    true_delays_chips, true_dopplers_hz, noise_std, ...
    pll_params, dll_params, use_sdcp) 

    if nargin < 11
        use_sdcp = false;
    end

    % multisat executor for track_one_satellite
    
    Nsat = numel(prn_pm_list);
    
    % Run first to setup prealloc
    prn_pm1 = prn_pm_list{1};
    delay1 = true_delays_chips(1);
    doppler1 = true_dopplers_hz(1);
    
    res1 = track_one_satellite( ...
        prn_pm1, chip_rate, fc, fs, duration, ...
        delay1, doppler1, noise_std, ...
        pll_params, dll_params, use_sdcp);
        
    % Preallocate
    t_blocks = res1.t_blocks_s;
    Nblock = numel(t_blocks);
    
    pseudorange_m = zeros(Nblock, Nsat);
    carrier_phase_rad = zeros(Nblock, Nsat);
    carrier_phase_cycles = zeros(Nblock, Nsat);
    doppler_est_hz = zeros(Nblock, Nsat);
    code_phase_chips = zeros(Nblock, Nsat);
    dll_disc = zeros(Nblock, Nsat);
    
    % Store first result
    pseudorange_m(:,1) = res1.pseudorange_m(:);
    carrier_phase_rad(:,1) = res1.carrier_phase_rad(:);
    carrier_phase_cycles(:,1) = res1.carrier_phase_cycles(:);
    doppler_est_hz(:,1) = res1.doppler_est_hz(:);
    code_phase_chips(:,1) = res1.code_phase_chips(:);
    dll_disc(:,1) = res1.dll_disc(:);
    
    % Track remaining satellites
    for s = 2:Nsat
        prn_pm  = prn_pm_list{s};
        delay_s = true_delays_chips(s);
        dop_s   = true_dopplers_hz(s);
        
        res = track_one_satellite( ...
            prn_pm, chip_rate, fc, fs, duration, ...
            delay_s, dop_s, noise_std, ...
            pll_params, dll_params, use_sdcp); 
            
        % timing array guard
        if numel(res.t_blocks_s) ~= Nblock
            error('All satellites must use same int duration');
        end
        
        pseudorange_m(:,s) = res.pseudorange_m(:);
        carrier_phase_rad(:,s) = res.carrier_phase_rad(:);
        carrier_phase_cycles(:,s) = res.carrier_phase_cycles(:);
        doppler_est_hz(:,s) = res.doppler_est_hz(:);
        code_phase_chips(:,s) = res.code_phase_chips(:);
        dll_disc(:,s) = res.dll_disc(:);
    end
    
    % Pack struct
    gnss.t_blocks_s        = t_blocks;
    gnss.pseudorange_m     = pseudorange_m;
    gnss.carrier_phase_rad = carrier_phase_rad;
    gnss.carrier_phase_cycles = carrier_phase_cycles;
    gnss.doppler_est_hz    = doppler_est_hz;
    gnss.code_phase_chips  = code_phase_chips;
    gnss.dll_disc          = dll_disc;
    gnss.true_delays_chips = true_delays_chips;
    gnss.true_dopplers_hz  = true_dopplers_hz;
    gnss.wavelength_m      = res1.wavelength_m; 
end
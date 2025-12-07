function [true_delays_chips, true_dopplers_hz, ranges_m] = ...
    compute_gnss_delays_dopplers(r_sat_eci_km, v_sat_eci_kmps, ...
                                 r_user_eci_km, v_user_eci_kmps, ...
                                 fc_hz, chip_rate)

    % get true delay and Doppler (if used) for constellation

    % consts
    c_mps = 299792458;
    c_kmps = c_mps / 1000;

    Nsat = size(r_sat_eci_km, 2);

    true_delays_chips = zeros(1, Nsat);
    true_dopplers_hz = zeros(1, Nsat);
    ranges_m = zeros(1, Nsat);

    for s = 1:Nsat
        rs = r_sat_eci_km(:, s);
        vs = v_sat_eci_kmps(:, s);

        rho_vec_km = rs - r_user_eci_km;
        range_km = norm(rho_vec_km);
        los_hat = rho_vec_km / range_km;

        v_rel_kmps = vs - v_user_eci_kmps;
        range_rate_kmps = dot(v_rel_kmps, los_hat);

        ranges_m(s) = range_km * 1000;
        tau_s = ranges_m(s) / c_mps;
        true_delays_chips(s) = tau_s * chip_rate;

        % prob wont use but doppler: f_D = -(fc/c) * range_rate
        true_dopplers_hz(s) = -(fc_hz / c_kmps) * range_rate_kmps;
    end
end

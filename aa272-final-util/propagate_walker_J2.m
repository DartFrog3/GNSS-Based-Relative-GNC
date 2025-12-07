function [r_eci, v_eci] = propagate_walker_J2(sats, t_vec, mu, Re, J2)
    % Analytic j2 prop for all constellation members, returns ECI

    N = numel(sats);
    Nt = numel(t_vec);

    r_eci = zeros(3, N, Nt);
    v_eci = zeros(3, N, Nt);

    for k = 1:N
        a_km  = sats(k).a_km;
        e = sats(k).e;
        i_rad = deg2rad(sats(k).i_deg);
        Om0 = deg2rad(sats(k).RAAN_deg);
        w0 = deg2rad(sats(k).w_deg);
        M0 = deg2rad(sats(k).M_deg);

        a = a_km;
        n = sqrt(mu / a^3);
        p = a * (1 - e^2);

        % from 279d
        prefactor = J2 * (Re^2 / p^2);

        Om_dot = -1.5 * n * prefactor * cos(i_rad);
        w_dot = 0.75 * n * prefactor * (5*cos(i_rad)^2 - 1);
        M_dot = n + 0.75 * n * prefactor * sqrt(1 - e^2) * (3*cos(i_rad)^2 - 1);

        for it = 1:Nt
            t = t_vec(it);

            Om_t = Om0 + Om_dot * t;
            w_t = w0 + w_dot * t;
            M_t = M0 + M_dot * t;

            M_wrap = wrapTo2Pi(M_t);
            E = kepler_E(M_wrap, e);
            f = 2*atan2( sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2) );

            % OE2ECI accepts deg
            i_deg_t = rad2deg(i_rad);
            Om_deg_t = rad2deg(wrapTo2Pi(Om_t));
            w_deg_t = rad2deg(wrapTo2Pi(w_t));
            f_deg_t = rad2deg(wrapTo2Pi(f));

            [r, v] = OE2ECI(a_km, e, i_deg_t, Om_deg_t, w_deg_t, f_deg_t);

            r_eci(:, k, it) = r;
            v_eci(:, k, it) = v;
        end
    end
end

function E = kepler_E(M, e)
    % Newton-raphson for E anom as in hw

    E = M;
    for iter = 1:20 % prob dont need a max_iter for now
        f = E - e*sin(E) - M;
        fp = 1 - e*cos(E);
        dE = -f / fp;
        E = E + dE;
        if abs(dE) < 1e-12
            break;
        end
    end
end

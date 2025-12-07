function sats = generate_walker_constellation( ...
    a_km, e, i_deg, RAAN0_deg, w_deg, Mref_deg, ...
    N, T, F)
    % create walker-delta constellation

    S = N / T;
    dRAAN = 360 / T;

    sats(1:N) = struct('a_km', [], 'e', [], 'i_deg', [], ...
                       'RAAN_deg', [], 'w_deg', [], 'M_deg', []);

    for k = 0:N-1
        p = floor(k / S);
        j = mod(k, S);

        RAAN_k = RAAN0_deg + p * dRAAN;
        M_k    = Mref_deg + (360/N) * k + (360*F/T) * p;

        % set OEs
        sats(k+1).a_km = a_km;
        sats(k+1).e = e;
        sats(k+1).i_deg = i_deg;
        sats(k+1).RAAN_deg = wrapTo360(RAAN_k);
        sats(k+1).w_deg = wrapTo360(w_deg);
        sats(k+1).M_deg = wrapTo360(M_k);
    end
end

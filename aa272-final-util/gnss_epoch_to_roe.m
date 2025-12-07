%% Simple conversion from ECI+chief to mean ROE

% inclusions TODO

function [roe_meas, chief_OE_mean, dep_OE_mean] = gnss_epoch_to_roe( ...
    r_chief_ECI, v_chief_ECI, ...
    r_dep_ECI,   v_dep_ECI, ...
    muearth)

    %% 1. Osculating OEs from ECI
    SVc = ECI2OE(r_chief_ECI, v_chief_ECI, muearth);
    SVd = ECI2OE(r_dep_ECI,   v_dep_ECI,   muearth);

    % Chief
    ac = SVc.a;
    ec = SVc.e;
    ic = SVc.i;
    Omc = SVc.Om;
    wc = wrapTo360(SVc.w);
    fc = wrapTo360(SVc.anom);
    Mc = rad2deg(true2mean(deg2rad(fc), ec));

    % Deputy
    ad = SVd.a;
    ed = SVd.e;
    id = SVd.i;
    Omd = SVd.Om;
    wd = wrapTo360(SVd.w);
    fd = wrapTo360(SVd.anom);
    Md = rad2deg(true2mean(deg2rad(fd), ed));

    %% 2. osc oe states
    SVc_osc = [ac*1e3; ec; deg2rad(ic); deg2rad(Omc); deg2rad(wc); deg2rad(Mc)];
    SVd_osc = [ad*1e3; ed; deg2rad(id); deg2rad(Omd); deg2rad(wd); deg2rad(Md)];

    % Convert to mean elements
    SVc_mean_noUnits = osc2mean(SVc_osc);
    SVd_mean_noUnits = osc2mean(SVd_osc);

    % Chief
    a_c_mean_km = SVc_mean_noUnits(1)/1e3;
    e_c_mean = SVc_mean_noUnits(2);
    i_c_mean = rad2deg(SVc_mean_noUnits(3));
    Om_c_mean = wrapTo360(rad2deg(SVc_mean_noUnits(4)));
    w_c_mean = rad2deg(SVc_mean_noUnits(5));
    f_c_mean = rad2deg(wrapTo2Pi(mean2true(SVc_mean_noUnits(6), e_c_mean)));

    chief_OE_mean = [a_c_mean_km; e_c_mean; i_c_mean; Om_c_mean; w_c_mean; f_c_mean];

    % Deputy
    a_d_mean_km = SVd_mean_noUnits(1)/1e3;
    e_d_mean = SVd_mean_noUnits(2);
    i_d_mean = rad2deg(SVd_mean_noUnits(3));
    Om_d_mean = wrapTo360(rad2deg(SVd_mean_noUnits(4)));
    w_d_mean = rad2deg(SVd_mean_noUnits(5));
    f_d_mean = rad2deg(wrapTo2Pi(mean2true(SVd_mean_noUnits(6), e_d_mean)));

    dep_OE_mean = [a_d_mean_km; e_d_mean; i_d_mean; Om_d_mean; w_d_mean; f_d_mean];

    %% 3. OE2QNSOE

    % Chief
    SVc_QNSOE_struct = OE2QNSOE( ...
        chief_OE_mean(1), chief_OE_mean(2), chief_OE_mean(3), ...
        chief_OE_mean(4), chief_OE_mean(5), chief_OE_mean(6));
    SVc_QNSOE = [SVc_QNSOE_struct.a, ...
                 SVc_QNSOE_struct.ex, ...
                 SVc_QNSOE_struct.ey, ...
                 SVc_QNSOE_struct.inc, ...
                 SVc_QNSOE_struct.RAAN, ...
                 wrapTo360(SVc_QNSOE_struct.u)];

    % Deputy
    SVd_QNSOE_struct = OE2QNSOE( ...
        dep_OE_mean(1), dep_OE_mean(2), dep_OE_mean(3), ...
        dep_OE_mean(4), dep_OE_mean(5), dep_OE_mean(6));
    SVd_QNSOE = [SVd_QNSOE_struct.a, ...
                 SVd_QNSOE_struct.ex, ...
                 SVd_QNSOE_struct.ey, ...
                 SVd_QNSOE_struct.inc, ...
                 SVd_QNSOE_struct.RAAN, ...
                 wrapTo360(SVd_QNSOE_struct.u)];

    %% 4. QNSOE2ROE
    roe_meas = QNSOE2ROE(SVc_QNSOE, SVd_QNSOE);
    roe_meas = roe_meas(:);
end

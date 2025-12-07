function NOE = ROE2QNSOE(a_o, e_x_o, e_y_o, i_o, RAAN_o, u_o, d_a, d_lambda, d_e_x, d_e_y, d_i_x, d_i_y)
    % Convert Relative Orbital Elements (ROE) to Quasi Nonsingular Orbital Elements (QNSOE). 
    % Inputs:
    % QNSOE of the chief orbit (orbit o) (km, unitless, unitless, deg, deg, deg)
    % ROE of the deputy orbit (meters)
    a_o = a_o*1e3; % convert to meters
    a_t = (a_o + d_a)/1e3; % convert SMA to kilometers
    e_x_t = (d_e_x/a_o) + e_x_o;
    e_y_t = (d_e_y/a_o) + e_y_o;
    i_t = rad2deg((d_i_x/a_o)) + i_o;
    RAAN_t = rad2deg(d_i_y/(a_o*sind(i_o))) + RAAN_o;
    u_t = rad2deg((d_lambda/a_o)) + u_o - (RAAN_t-RAAN_o)*cosd(i_o);
    NOE = [a_t, e_x_t, e_y_t, i_t, RAAN_t, u_t];
end
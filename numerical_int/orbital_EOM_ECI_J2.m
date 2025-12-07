function statedot = orbital_EOM_ECI_J2(t, state)
    % Returns derivative of earth satellite state with J2 perturbation
    % state = [x; y; z; vx; vy; vz] in km, km/s

    % Earth physical constants from Vallado
    mu = 3.986 * 10^5;              % [km^3/sec^2]
    Re =   6378.1363;               % [km] mean equatorial radius
    J2 =      0.0010826266;        % J2

    r = state(1:3);
    v = state(4:6);
    x = r(1);  y = r(2);  z = r(3);
    r_norm = norm(r);

    acc_grav = -mu / r_norm^3 * r;

    % J2 perturbation - check these some conflicting
    factor = 3*J2*mu*Re^2 / (2 * r_norm^5);
    ax_J2 = factor * x * (5*(z^2)/(r_norm^2) - 1);
    ay_J2 = factor * y * (5*(z^2)/(r_norm^2) - 1);
    az_J2 = factor * z * (5*(z^2)/(r_norm^2) - 3);
    acc_J2 = [ax_J2; ay_J2; az_J2];

    acc = acc_grav + acc_J2;

    statedot = [v; acc];
end

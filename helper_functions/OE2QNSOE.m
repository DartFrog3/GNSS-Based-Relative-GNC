function oe = OE2QNSOE(a, e, inc, RAAN, w, f)
    %OE2QNSOE Convert classical orbital elements to quasi-nonsingular orbital elements
    %
    % Inputs:
    %   a     - semi-major axis (km)
    %   e     - eccentricity (unitless)
    %   inc   - inclination (degrees)
    %   RAAN  - right ascension of ascending node (degrees)
    %   w     - argument of periapsis (degrees)
    %   f     - true anomaly (degrees)
    %
    % Outputs:
    %   oe    - struct containing quasi-nonsingular orbital elements:
    %           .a     - semi-major axis (km)
    %           .ex    - x-component of eccentricity vector (unitless)
    %           .ey    - y-component of eccentricity vector (unitless)
    %           .inc   - inclination (degrees)
    %           .RAAN  - right ascension of ascending node (degrees)
    %           .u     - mean argument of latitude (degrees)

    oe.a = a;
    oe.ex = e*cosd(w);
    oe.ey = e*sind(w);
    oe.inc = inc;
    oe.RAAN = RAAN;
    % Convert true anomaly to mean anomaly
    M = true2mean(deg2rad(f), e); 
    oe.u = rad2deg(M) + w;
end
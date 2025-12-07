
function [a, e, inclination, RAAN, w, f, M] = QNSOE2OE(a, ex, ey, inclination, RAAN, u)
    % QNSOE2OE converts quasi-nonsingular orbital elements (QNSOE) to classical orbital elements (COE).
    %
    % Inputs:
    %   a          - Semi-major axis (km)
    %   ex         - Eccentricity component along x-axis (dimensionless)
    %   ey         - Eccentricity component along y-axis (dimensionless)
    %   inclination - Orbital inclination (degrees)
    %   RAAN       - Right Ascension of the Ascending Node (degrees)
    %   u          - Argument of latitude (degrees)
    %
    % Outputs:
    %   a          - Semi-major axis (km)
    %   e          - Eccentricity (dimensionless)
    %   inclination - Orbital inclination (degrees)
    %   RAAN       - Right Ascension of the Ascending Node (degrees)
    %   w          - Argument of periapsis (degrees)
    %   f          - True anomaly (degrees)
    %   M          - Mean anomaly (degrees)
    % Notes:
    %   - The conversion from mean anomaly to true anomaly is performed using an iterative method
    %     with a tolerance of 1e-10.

    tol = 1e-10;
    % Eccentricity
    e = sqrt(ex^2 + ey^2);

    % Argument of periapsis
    w = atan2d(ey, ex);

    % Mean anomaly
    M = u - w;

    % True anomaly
    f = rad2deg(mean2true(deg2rad(M), e, tol));
end
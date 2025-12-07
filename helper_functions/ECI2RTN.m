function [rot_RTN, w_RTN] = ECI2RTN(r_ECI, v_ECI, a, e, f)
% createRTNFrame
% Returns the rotation matrix and angular velocity vector (in RTN) that
% transform vectors from ECI to RTN coordinates
%
% Inputs:
%   r_ECI  - 3x1 position vector in ECI [km]
%   v_ECI  - 3x1 velocity vector in ECI [km/s]
%   a      - semi-major axis [km]
%   e      - eccentricity
%   f_deg  - true anomaly [deg]
%
% Outputs:
%   rot_RTN - 3x3 rotation matrix from ECI to RTN
%   w_RTN   - 3x1 angular velocity vector of RTN frame w.r.t. ECI, in RTN frame

    mu = 3.986004418e5;  % [km^3/s^2] Earth's gravitational parameter

    % Angular Velocity Vector in RTN
    h = cross(r_ECI, v_ECI);

    % Find RTN unit vectors in xyz frame
    R = r_ECI / norm(r_ECI);
    N = h / norm(h);
    T = cross(N, R);
    
    % xyz frame
    X = [1 0 0]';
    Y = [0 1 0]';
    Z = [0 0 1]';

    % find rotation matrix from ECI to RTN
    rot_RTN = [dot(R,X) dot(R,Y) dot(R,Z); dot(T,X) dot(T,Y) dot(T,Z); dot(N,X) dot(N,Y) dot(N,Z)];

    % determine angular velocity of ECI wrt. RTN 
    f_dot = sqrt(mu/a^3/(1 - e^2)^3)*(1 + e*cosd(f))^2;
    w_RTN = [0 0 f_dot]';
end

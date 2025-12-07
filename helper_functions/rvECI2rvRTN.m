function [r_RTN, v_RTN] = rvECI2rvRTN(r_ECI, v_ECI, rot_RTN, w_RTN)
% convertECItoRTN
% Converts ECI vectors into RTN frame using precomputed rotation matrix and angular velocity.
%
% Inputs:
%   r_ECI   - 3x1 position in ECI [km]
%   v_ECI   - 3x1 velocity in ECI [km/s]
%   rot_RTN - 3x3 rotation matrix from ECI to reference RTN frame
%   w_RTN   - 3x1 angular velocity of RTN frame [rad/s]
%
% Outputs:
%   r_RTN - 3x1 position in RTN frame [km]
%   v_RTN - 3x1 velocity in RTN frame [km/s]

    % from lec2 slide 12

    r_RTN = rot_RTN * r_ECI;
    v_RTN = rot_RTN * v_ECI - cross(w_RTN, r_RTN);
end

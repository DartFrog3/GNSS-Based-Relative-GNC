function [statedot] = orbital_EOM_ECI(t, state)
    % Returns derivative of earth satellite state
    % State is [rx ry rz vx vy vz]' in [km] and [km/sec]
    %
    % AA 279A Winter 2024
    % Andrew K. Barrows
    % Earth physical constants from Vallado
    muearth = 3.986 * 10^5;              % [km^3/sec^2]
    
    r = state(1:3);
    v = state(4:6);
    
    normr = norm(r);
    rhat = r/normr;
    
    % Acceleration due to 1/r^2 gravity
    acc = - (muearth/(normr*normr))*rhat;
           
    statedot = zeros(6, 1);
    statedot(1:3) = v;
    statedot(4:6) = acc;

end % terminates MATLAB function



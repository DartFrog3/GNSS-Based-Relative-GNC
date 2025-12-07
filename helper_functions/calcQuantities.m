function [e, h, sme] = calcQuantities(y)
    % calcQuantities returns the eccentricity vector, the angular momentum vector, and the specific mechanical energy given a state[r; v] 
    % Inputs :
    % y = [r; v]
    % Outputs :
    % e - eccentricity vector
    % h - angular momentum vector  [km^2 / sec]
    % sme - specific mechanical energy  [km^2 / sec^2]

    % Earth physical constants from Vallado
    muearth = 3.986 * 10^5;              % [km^3/sec^2]
    

    rVec = [y(:,1) y(:,2) y(:,3)];
    vVec = [y(:,4) y(:,5) y(:,6)];
    rNorm = sqrt(y(:,1).^2 + y(:,2).^2 + y(:,3).^2); % [km]
    
    h = cross(rVec, vVec, 2);
    
    e = cross(vVec, h, 2) / muearth - rVec ./ rNorm;
    
    
    speed = sqrt(y(:,4).^2 + y(:,5).^2 + y(:,6).^2); % [km/sec] 
    ske = speed.^2/2;  % specific kinetic energy [km^2/sec^2]
    spe = -muearth./rNorm; % specific potential energy [km^2/sec^2]
    sme = ske + spe;   % specific mechanical energy [km^2/sec^2]

end
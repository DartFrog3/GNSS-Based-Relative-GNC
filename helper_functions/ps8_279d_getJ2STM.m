function [STM] = ps8_279d_getJ2STM(tau, a, ex, ey, inc)
    %  propagation of QNROE under J2 using STM
    %
    %  Inputs:
    %  tau   – time [ sec ]
    %  a   – semi-major axis of the chief [ km ]
    %  ex  – chief eccentricity vector, x-comp  
    %  ey  – chief eccentricity vector, y-comp
    %  inc  – inclination of the chief [ deg ]
    %
    %  Output:
    %  STM - stm under J2
    
    %% constants and setup quantities
    % Earth physical constants from Vallado
    mu = 3.986 * 10^5;              % [km^3/sec^2]
    Re =   6378.1363;               % [km] mean equatorial radius
    J2 =      0.0010826266;        % J2
    
    % negative with retrograde??? --> POTENTIALLY FIX THIS
    w = atan2d(ey, ex);
    n = sqrt(mu / a^3);
    e = sqrt(ex^2 + ey^2);
    
    eta = sqrt(1 - e^2);
    kap = 0.75 * J2 * (Re^2) * sqrt(mu) / (a^(7/2) * eta^4);
    Epar = 1 + eta;
    Fpar = 4 + 3*eta;
    Gpar = 1/eta^2;
    
    Ppar = 3*cosd(inc)^2 - 1;
    Qpar = 5*cosd(inc)^2 - 1;
    Rpar = cosd(inc);
    Spar = sind(2*inc);
    Tpar = sind(inc)^2;
    
    wdot = kap * Qpar * 180 / pi;
    wf = w + wdot * tau;
    
    ex_f = e * cosd(wf);
    ey_f = e * sind(wf);
    ex_i = ex;
    ey_i = ey;
    %% build STM
    STM = zeros(6);
    
    % dela row
    STM(1,1) = 1;
    
    % delLambda row
    STM(2,1) = -(1.5*n + 3.5*kap*Epar*Ppar)*tau;
    STM(2,2) = 1;
    STM(2,3) = kap*ex_i*Fpar*Gpar*Ppar*tau;
    STM(2,4) = kap*ey_i*Fpar*Gpar*Ppar*tau;
    STM(2,5) = -kap*Fpar*Spar*tau;
    
    % delex row
    STM(3,1) = 3.5*kap*ey_f*Qpar*tau;
    STM(3,3) = cosd(wdot*tau)-4*kap*ex_i*ey_f*Gpar*Qpar*tau;
    STM(3,4) = -sind(wdot*tau)-4*kap*ey_i*ey_f*Gpar*Qpar*tau;
    STM(3,5) = 5*kap*ey_f*Spar*tau;
    
    % deley row
    STM(4,1) = -3.5*kap*ex_f*Qpar*tau;
    STM(4,3) = sind(wdot*tau)+4*kap*ex_i*ex_f*Gpar*Qpar*tau;
    STM(4,4) = cosd(wdot*tau) + 4*kap*ey_i*ex_f*Gpar*Qpar*tau;
    STM(4,5) = -5*kap*ex_f*Spar*tau;
    
    % delix
    STM(5,5) = 1;
    
    % deliy row
    STM(6,1)       =  3.5*kap*Spar*tau;
    STM(6,3)       = -4*kap*ex_i*Gpar*Spar*tau;
    STM(6,4)       = -4*kap*ey_i*Gpar*Spar*tau;
    STM(6,5)       =  2*kap*Tpar*tau;
    STM(6,6)       =  1;
end
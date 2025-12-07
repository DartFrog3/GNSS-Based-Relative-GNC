function [Ac] = ps6_279d_getReducedSTM(a, ex, ey, inc, w)
    % Returns reduced STM for given OEs
    % Inputs :
    % a = sma  [ km ]
    % ex = eccentricity x initial
    % ey = eccentricity y initial
    % inc = inclination  [ deg ]
    % w = aop [ deg ]
    % Outputs :
    % Ac = reduced STM for given OEs
    mu = 3.986 * 10^5;
    J2 = 0.0010826266;
    Re = 6378.1363;
    n = sqrt(mu / a^3);

    % setup quantities
    
    e_mag = hypot(ex, ey);
    eta = sqrt(1 - e_mag^2);
    G = 1/eta^2;
    
    gamma = (3/4)*J2*Re^2*sqrt(mu);
    kappa = gamma / (a^(7/2)*eta^4);
    P = 3*cosd(inc)^2 - 1;
    Q = 5*cosd(inc)^2 - 1;
    T = sind(inc)^2;
    S = sind(2*inc);
    C = sind(w);
    D = sind(w);
    
    %% Set phi entries
    Phi = zeros(5, 5);


    % KEEP WITH 5x5 for unconstrained case
    % Row 1
    %Phi(1, 6) = 1/kappa;
    
    % Row 2
    Phi(2,1) = 3.5*ey*Q;
    Phi(2, 2) = -(4*ex*ey*G + C)*Q;
    Phi(2,3) = -(1 + 4*ey^2*G - D)*Q;
    Phi(2,4) =  5*ey*S;
    %Phi(2, 6) = (D-ex)/kappa;
    
    % Row 3
    Phi(3,1) =  -3.5*ex*Q;
    Phi(3, 2) = (1 + 4*ex^2*G - D)*Q;
    Phi(3,3) =  (4*ex*ey*G - C)*Q;
    Phi(3,4) = -5*ex*S;
    %Phi(3, 6) = (C-ey)/kappa;
 
    % Row 5
    Phi(5,1) =  3.5*S;
    Phi(5,2) = -4*ex*G*S;
    Phi(5,3) = -4*ey*G*S;
    Phi(5,4) =  2*T;

    Ac = kappa*Phi;
end
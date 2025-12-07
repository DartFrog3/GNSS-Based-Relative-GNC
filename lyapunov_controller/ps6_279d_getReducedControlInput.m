function [gamma] = ps6_279d_getReducedControlInput(a, ex, ey, inc, w, f)
    % Returns Control Input Matrix for HCW Model
    % Inputs :
    % a = sma  [ km ]
    % ex = ecc x
    % ey = ecc y
    % inc = inclinatin [ deg ]
    % w = aop [ deg ]
    % f = true anom [ deg ]
    % Outputs :
    % gamma = reduced dynamics control input matrix at f

    muearth = 3.986 * 10^5;
    n = sqrt(muearth / a^3);
    ecc = hypot(ex, ey);
    eta = sqrt(1 - ecc^2);
    den = 1 + ecc*cosd(f);

    
    % define structure
    B = zeros(5,2);
    
    B(1,1) = (2/eta) * den;
    
    B(2,1) = eta * ((2 + ecc*cosd(f))*cosd(w+f) + ex) / den;
    B(2, 2) = (eta * ey / tand(inc)) * sind(w+f)/den;
    
    B(3,1) = eta * ((2 + ecc*cosd(f))*sind(w+f) + ey) / den;
    B(3,2) = (-eta * ex / tand(inc)) * sind(w+f)/den;
    
    B(4,2) = eta * cosd(w + f) / den;

    B(5,2) = eta * sind(w + f) / den;
     
    gamma = B / (n*a);
end
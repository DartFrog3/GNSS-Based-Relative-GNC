function [a , e , i , Om , w , f]=QN2OE(bet, gam, del, eps, zet, eta)
    % QN2OE Converts quasi-nonsingular elements to standard orbital elements %
    % Inputs :
    % bet = a  [ km ]
    % gam = e cos(w)
    % del = e sin(w)
    % eps = i [ deg ]
    % zet = Om [ deg ]
    % eta = w + M [ deg ]
    % Outputs :
    % a - semi - major axis of orbit [ km ]
    % e - eccentricity of orbit
    % i - inclination of orbit [ deg ]
    % Om - right ascension of the ascending node [ deg ]
    % w - argument of periapsis [ deg ]
    % v - true anomaly [ deg ] %

    a = bet;
    i = eps;
    Om = zet;
    e = gam^2 + del^2;
    w = atan2d(del, gam);
    M = eta - w;

    % External call takes radians
    Mrad = deg2rad(M);
    addpath('./mean_osc/');

    % CHECK THIS TOL IS APPROPRIATE
    f = rad2deg(mean2true(Mrad, e, 10^-5));
    
end

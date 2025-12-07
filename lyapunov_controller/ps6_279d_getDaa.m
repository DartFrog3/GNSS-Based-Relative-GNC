function [daa] = ps6_279d_getDaa(a, ecc, M, deltaDlambda, ddeNorm, ddiNorm)
    % Returns reduced dda for given dlambda
    % Inputs :
    % a = sma  [ km ]
    % ecc = eccentricity
    % M = mean anomaly [ deg ]
    % deltaDlambda = delta change mean longitude
    % controlInp = tangential control input as hyperparam for control input size
    % Outputs :
    % daa
    
    mu = 3.986 * 10^5;
    n = sqrt(mu / a^3);
    eta = sqrt(1 - ecc^2);
    T = 2*pi*sqrt(a^3/mu);

    % params
    tau = T/4; % THIS PARAMETER IS LITERALLY THE ANSWER - higher is lower delta v but cruder maneuvers while lower is higher dv but more accurate
    
    %% get u_d
    v_opt_ip = 0.5 * a*n * ddeNorm / eta;
    v_opt_oop = a*n * (1 - ecc) * ddiNorm / eta;

    nOrbits = 1; % TAKING 6 TO BE NUMBER OF ORBITS like in Tau BUT WHAT IS THIS
    nOrbits_ip = ddeNorm / (2*(1-ecc)*ddiNorm + ddeNorm) * nOrbits; 
    nOrbits_oop = 2*(1-ecc)*ddiNorm / (2*(1-ecc)*ddiNorm + ddeNorm) * nOrbits;

    u_d_ip = 2*v_opt_ip / (T*nOrbits_ip) * (8/3);
    u_d_oop = 2*v_opt_oop / (T*nOrbits_oop) * (8/3);

    u_d = u_d_ip + u_d_oop;
    %% solve
    % setup quantities
    f = rad2deg(mean2true(deg2rad(M), ecc));
    
    dvtan = u_d*T/(4) * (3/8); %for N>4 then for i=1:4 prod = prod * product_i^4 (q-1)/q - WHY DOES IT SAY STOP AT 4????? % need to take N parm as input?

    deltaDatan = (2/(a*n*eta))*(1+ecc*cosd(f))*dvtan;
    da_ref = (deltaDatan)/2;
    dlambdaDotref = 1.5*n*abs(da_ref);

    if deltaDlambda >= 0
        dlambdaDot = -min(abs(deltaDlambda)/tau, dlambdaDotref);
    else
        dlambdaDot = min(abs(deltaDlambda)/tau, dlambdaDotref);
    end
    
    daa = -(2/3)*dlambdaDot / n;
end
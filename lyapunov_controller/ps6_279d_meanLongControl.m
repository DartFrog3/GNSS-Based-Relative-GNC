function [statedot, u, trackErr, da_des] = ps6_279d_meanLongControl(t, state)
    % Returns derivative of the relative state for reduced continuous
    % control with mean longitude control
    % state = [chief_oe; roe; state_ref] where here
    % state_ref is the desired state [da; dlambda; dex; dey; dix; diy]

    % Earth physical constants from Vallado
    mu = 3.986 * 10^5;              % [km^3/sec^2]
    Re =   6378.1363;               % [km] mean equatorial radius
    J2 = 0.0010826266;
    
    %% Setup quantities
    chief_oe = state(1:6);
    state_ref = [state(13); state(15:18)];
    dAlpha = [state(7); state(9:12)];
    
    % define OEs
    a = chief_oe(1);
    e = chief_oe(2);
    inc = chief_oe(3);
    RAAN = chief_oe(4);
    w = chief_oe(5);
    M = chief_oe(6);
    ex = e*cosd(w);
    ey = e*sind(w);

    % get f
    f = rad2deg(mean2true(deg2rad(M), e));

    % define ROEs
    da = dAlpha(1);
    dlambda = state(8);
    dex = dAlpha(2);
    dey = dAlpha(3);
    dix = dAlpha(4);
    diy = dAlpha(5);

    % get dlambda change
    deltaDlambda =  dlambda - state(14);
    
    % get de, di changes
    ddeNorm = hypot(state_ref(2)-dex, state_ref(3)-dey);
    ddiNorm = hypot(state_ref(4)-dix, state_ref(5)-diy);

    da_des = ps6_279d_getDaa(a, e, M, deltaDlambda, ddeNorm, ddiNorm);
    state_ref(1) = da_des;
    
    % get Ac, Bc
    Ac = ps6_279d_getReducedSTM(a, ex, ey, inc, w); 
    Bc = ps6_279d_getReducedControlInput(a, ex, ey, inc, w, f); % these might be a bit small, but produce reasonable control input so prob ok

    % get P - 
    u_ip = rad2deg(atan2(dey, dex));  
    u_oop = rad2deg(atan2(diy, dix));
    uarg = w + M;

    J = rad2deg(wrapTo2Pi(uarg-u_ip));
    H = rad2deg(wrapTo2Pi(uarg-u_oop));

    P = ps6_279d_getControlGain(J, H);

    %% calculate control input and apply
    u = -Bc\(Ac * dAlpha + P*(dAlpha - state_ref)); 
    dAlphaDot = Ac*dAlpha + Bc*u;

    trackErr = [state_ref(1); state(14); state_ref(2:5)]  - [dAlpha(1); state(8); dAlpha(2:5)];

    %% propagate chief_oes with J2
    % Calculate the drift rates
    % Here we want to calculate the drift rates of RAAN, AOP, and M
    n   = sqrt(mu/a^3);                  % mean motion
    eta   = sqrt(1-e^2);
    gamma = 3/4 * J2*Re^2*sqrt(mu);
    kappa = gamma/(a^(7/2)*eta^4);
    Ppar     = 3*cosd(inc)^2 - 1;
    Q     = 5*cosd(inc)^2 - 1;

    % this computionally wasteful? should do outside this function in
    % theory and not calculate every step
    RAAN_dot = -2*cosd(inc)*kappa*(180/pi);
    aop_dot  = kappa * Q*(180/pi);
    M_dot    = (n + kappa*eta*Ppar)*(180/pi);

    chief_oeDot = [0;0;0;RAAN_dot;aop_dot;M_dot];


    prefactor = 1/2 * J2*Re^2/(a^2*eta^4);
    dld1 = prefactor*cosd(2*inc)*dix;
    dld2 = da/7;
    dlambdaDot = -1.5*n*da*180/pi;
    %dlambdaDot = -10.5*(dld1 + dld2)*180/pi;

    %% return statedot
    statedot = [chief_oeDot; dAlphaDot(1); dlambdaDot; dAlphaDot(2:5); zeros(size([state_ref;0]))]; % as no change to state_ref and directly pop out new dda
end
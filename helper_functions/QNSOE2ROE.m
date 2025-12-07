function delta_alpha = QNSOE2ROE(chief_QNSOE, deputy_QNSOE)
    % QNSOE2ROE Computes the relative orbit elements (ROE) from two QNSOE states
    %
    % Inputs:
    %   chief_QNSOE  - 6x1 vector [a, ex, ey, i, RAAN, u] of the chief (reference) in km and degrees
    %   deputy_QNSOE - 6x1 vector [a, ex, ey, i, RAAN, u] of the deputy in km and degrees
    %
    % Output:
    %   delta_alpha  - 6x1 vector of relative orbit elements [da, dlambda, dex, dey, dix, diy] in meters
    
    % Extract chief elements
    [a_c, ex_c, ey_c, i_c, RAAN_c, u_c] = deal(chief_QNSOE(1), chief_QNSOE(2), chief_QNSOE(3), chief_QNSOE(4), chief_QNSOE(5), chief_QNSOE(6));
    
    % Extract deputy elements
    [a_d, ex_d, ey_d, i_d, RAAN_d, u_d] = deal(deputy_QNSOE(1), deputy_QNSOE(2), deputy_QNSOE(3), deputy_QNSOE(4), deputy_QNSOE(5), deputy_QNSOE(6));
    
    % Compute relative orbit elements
    da = (a_d - a_c) / a_c;
    dlambda = deg2rad(wrapTo180(u_d - u_c) + deg2rad(RAAN_d - RAAN_c) * cosd(i_c)); % ARE THE PARENS HERE RIGHT?
    dex = (ex_d - ex_c);
    dey = (ey_d - ey_c);
    dix = deg2rad(i_d - i_c);
    diy = deg2rad((RAAN_d - RAAN_c) * sind(i_c)); % ARE THE PARENS HERE RIGHT?
    
    % Output
    delta_alpha = a_c*1000*[da; dlambda; dex; dey; dix; diy]';
    
end
    
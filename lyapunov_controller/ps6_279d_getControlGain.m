function [P] = ps6_279d_getControlGain(J, H)
    % Returns Control Gain Matrix
    % Inputs :
    % J = \phi - \phi_opt_ip [ deg ]
    % H = \phi - \phi_opt_oop [ deg ]
    % Outputs :
    % P = control gain matrix
    N = 10;
    k = 1000;
    
    % define structure
    preP = [cosd(J)^N, 0, 0, 0, 0; ...
                0, cosd(J)^N, 0, 0, 0; ...
                0, 0, cosd(J)^N, 0, 0; ...
                0, 0, 0, cosd(H)^N, 0; ...
                0, 0, 0, 0, cosd(H)^N];

    P = preP / k;
end
function rot_RTN = ECI2RTN_mtrx(r_ECI, v_ECI)
% Find Matrix that converts a vector in ECI to a vector in RTN
    h = cross(r_ECI, v_ECI);

    % Find RTN unit vectors in xyz frame
    R = r_ECI / norm(r_ECI);
    N = h / norm(h);
    T = cross(N, R);
    
    % xyz frame
    X = [1 0 0]';
    Y = [0 1 0]';
    Z = [0 0 1]';

    rot_RTN = [dot(R,X) dot(R,Y) dot(R,Z); dot(T,X) dot(T,Y) dot(T,Z); dot(N,X) dot(N,Y) dot(N,Z)];
end
function [ r_eci , v_eci ] = OE2ECI2(a , e , i , O , w , v)
    % OE2ECI Converts orbital elements to r, v in ECI frame
    
    % Inputs :
    % a - semi - major axis of orbit [km]
    % e - eccentricity of orbit
    % i - inclination of orbit [ deg ]
    % O - right ascension of the ascending node [ deg ]
    % w - argument of periapsis [ deg ]
    % v - true anomaly [ deg ]
    
    % Outputs :
    % r_eci - 3x1 vector of radius in ECI frame [km]
    % v_eci - 3x1 vector of velocity in ECI frame [km/s]
    
    % Notes :
    % 1) This function assumes the central body is Earth
    % 2) In case of equatorial orbits , the sum `O + w`
    % represents the longitude of periapsis
    % 3) In case of circular orbits , the sum `w + v`
    % represents the argument of latitude
    % 4) In case of equatorial circular orbits , the sum
    % `O + w + v` represents the true longitude

    mu = 3.986e5 ; % [km ^3/ s ^2]
    p = a * (1 - e * e ); % semi - latus rectum [m]
    % Radius and velocity of orbit in perifocal coordinates
    r = p / (1 + e * cosd(v) ) ;% orbit radius [m]
    r_peri = r * [ cosd(v ); sind(v) ; 0];
    v_peri = sqrt ( mu / p) * [- sind(v ); e + cosd(v ); 0];
    % Rotate vectors into ECI frame
    R_peri_eci = rotz (O ) * rotx ( i) * rotz (w );
    r_eci = R_peri_eci * r_peri ;
    v_eci = R_peri_eci * v_peri ;

end
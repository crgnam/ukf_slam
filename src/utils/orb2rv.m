function [r, v] = orb2rv(a,e,inc,w,RAAN,M0,mu)
    % This function converts orbital elements to position and velocity.
    n = sqrt(mu/(a^3)); %Mean motion
    dt = 1;
    M = M0 + n*dt;

    % Calculate eccentric anomaly:
    KeplerConverged = 0;
    convergencePercentage = 0.05;
    E = 1;
    while KeplerConverged == 0
        % Set up Newton-Raphson method
        f = E - e*sind(E) - M;
        df = 1 - e*cosd(E);
        E_new = E - f / df;

        % Check for convergence
        relativeDifference = abs(E_new - E) / E * 100;
        if relativeDifference < convergencePercentage
            KeplerConverged = 1;
        end
        E = E_new;
    end

    r_mag = a*(1 - e*cosd(E));
    x = a*(cosd(E) - e);
    y = a*sqrt(1 - e^2)*sind(E);
    xDot = -sqrt(mu*a) / r_mag * sind(E);
    yDot = sqrt(mu*a*(1-e^2) ) /r_mag * cosd(E);

    a11 = cosd(RAAN)*cosd(w) -sind(RAAN)*sind(w)*cosd(inc);
    a12 = sind(RAAN)*cosd(w) + cosd(RAAN)*sind(w)*cosd(inc);
    a13 = sind(w)*sind(inc);

    a21 = -cosd(RAAN)*sind(w) - sind(RAAN)*cosd(w)*cosd(inc);
    a22 = -sind(RAAN)*sind(w) + cosd(RAAN)*cosd(w)*cosd(inc);
    a23 =  cosd(w) * sind(inc);

    a31 = sind(RAAN)*sind(inc);
    a32 = -cosd(RAAN)*sind(inc);
    a33 = cosd(inc);

    A = [a11, a12, a13; a21, a22, a23; a31, a32, a33];

    r = A'*[x; y; 0];
    v = A'*[xDot; yDot; 0];
end
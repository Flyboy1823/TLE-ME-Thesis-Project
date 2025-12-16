function [r_eci, v_eci] = OE2ECI_Test(a, e, i, RAAN, omega, M, rp)
    % Constants
    mu = 398600.4418; % Earth's gravitational parameter [km^3/s^2]
    E=zeros(size(a));
    nu=zeros(size(a));

    % Handle special cases
    for count=1:length(a)

        if e(count) == 0
            omega(count) = 0; % Argument of perigee undefined
        end
        if i(count) == 0
            RAAN(count) = 0; % RAAN undefined
        end
    
        % Solve Kepler's Equation for Eccentric Anomaly E
        E(count) = solveKeplersEquation(M(count), e(count));
        
        nu(count) = acos((cos(E(count)-e(count)))/(1-cos(E(count)*e(count))));

    end

%%%%%%%%%%%% Stopped here



    % Compute True Anomaly (nu)


    % Compute distance
    r = a .* (1 - e .* cos(E));

    % Perifocal position and velocity
    r_pf = r .* [cos(nu); sin(nu); 0];
    v_pf = sqrt(mu * a) / r * [-sin(E); sqrt(1 - e^2) * cos(E); 0];

    % Rotation matrices
    R3_W = [cos(-RAAN), -sin(-RAAN), 0;
            sin(-RAAN),  cos(-RAAN), 0;
                 0     ,     0     , 1];

    R1_i = [1,      0       ,       0     ;
            0, cos(-i), -sin(-i);
            0, sin(-i),  cos(-i)];

    R3_w = [cos(-omega), -sin(-omega), 0;
            sin(-omega),  cos(-omega), 0;
                0     ,      0      , 1];

    % Complete rotation from perifocal to ECI
    Q_pX = R3_W * R1_i * R3_w;

    % Position and velocity in ECI
    r_eci = Q_pX * r_pf;
    v_eci = Q_pX * v_pf;
end

function E = solveKeplersEquation(M, e)
    % Iterative solution using Newton-Raphson
    tol = 1e-10;
    E = M; % Initial guess
    ratio = 1;
    while abs(ratio) > tol
        ratio = (E - e * sin(E) - M) / (1 - e * cos(E));
        E = E - ratio;
    end
end

function [r_t1, v_t1] = propagate_orbit_j2_drag(r0, v0, t0, t1, Bstar)
% Propagate orbit from t0 to t1 using J2 and atmospheric drag (via B*)
%
% Inputs:
%   r0    - Initial position vector [km]
%   v0    - Initial velocity vector [km/s]
%   t0    - Initial time [s]
%   t1    - Final time [s]
%   Bstar - B* drag term from TLE [1/earth radii]
%
% Outputs:
%   r_t1  - Final position vector [km]
%   v_t1  - Final velocity vector [km/s]

    % Constants
    mu = 398600.4418;             % Earth gravitational parameter [km^3/s^2]
    Re = 6378.137;                % Earth radius [km]
    J2 = 1.08262668e-3;           % J2 coefficient
    CdA_m = Bstar * Re;           % Convert B* to [km^2/kg]

    % Simple atmospheric model (scale height approximation)
    rho0 = 3.614e-13;             % Reference density [kg/km^3]
    H = 88.667;                   % Scale height [km]

    % Initial state
    y0 = [r0(:); v0(:)];
    dt = t1 - t0;

    % Integrate using ode45
    opts = odeset('RelTol',1e-9,'AbsTol',1e-9);
    [~, Y] = ode45(@(t, y) ode_dynamics(t, y, mu, Re, J2, CdA_m, rho0, H), [0 dt], y0, opts);

    % Final state
    r_t1 = Y(end,1:3)';
    v_t1 = Y(end,4:6)';
end

function dydt = ode_dynamics(~, y, mu, Re, J2, CdA_m, rho0, H)
    r = y(1:3);
    v = y(4:6);

    r_norm = norm(r);
    x = r(1); y_ = r(2); z = r(3);

    % Gravitational acceleration
    a_grav = -mu / r_norm^3 * r;

    % J2 perturbation
    z2 = z^2;
    r2 = r_norm^2;
    factor = 1.5 * J2 * mu * Re^2 / r_norm^5;
    ax_J2 = factor * x * (5 * z2 / r2 - 1);
    ay_J2 = factor * y_ * (5 * z2 / r2 - 1);
    az_J2 = factor * z * (5 * z2 / r2 - 3);
    a_J2 = [ax_J2; ay_J2; az_J2];

    % Atmospheric drag
    v_rel = v;                   % Assume no wind in atmosphere
    v_rel_norm = norm(v_rel);
    altitude = r_norm - Re;
    rho = rho0 * exp(-altitude / H);
    a_drag = -0.5 * CdA_m * rho * v_rel_norm * v_rel;

    % Total acceleration
    a_total = a_grav + a_J2 + a_drag;

    dydt = [v; a_total];
end

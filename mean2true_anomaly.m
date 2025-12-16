function nu = mean2true_anomaly(M, e, tol)
    % mean2true_anomaly: Convert mean anomaly to true anomaly
    % Works for vectors or matrices of M and e (element-wise)
    %
    % Inputs:
    %   M   - mean anomaly [degree], array (same size as e or scalar)
    %   e   - eccentricity, array (same size as M or scalar)
    %   tol - (optional) convergence tolerance for Kepler iteration
    %
    % Output:
    %   nu  - true anomaly [deg], same size as M

    if nargin < 3
        tol = 1e-8;
    end

    % Ensure M and e are same size or one is scalar
    if ~isequal(size(M), size(e))
        if isscalar(M)
            M = repmat(M, size(e));
        elseif isscalar(e)
            e = repmat(e, size(M));
        else
            error('M and e must be the same size or one must be scalar.');
        end
    end

    % Transform M from Degrees to Rad
    M = deg2rad(M);

    % Normalize M to [0, 2*pi)
    M = mod(M, 2*pi);

    % Initial guess: E = M
    E = M;

    % Newton-Raphson iteration
    max_iter = 100;
    for k = 1:max_iter
        f = E - e .* sin(E) - M;
        fp = 1 - e .* cos(E);
        delta = f ./ fp;
        E = E - delta;

        if all(abs(delta(:)) < tol)
            break;
        end
    end

    % Compute true anomaly
    nu = 2 * atan2( sqrt(1+e) .* sin(E/2), sqrt(1-e) .* cos(E/2) );
    nu = mod(nu, 2*pi);  % Normalize to [0, 2*pi)
    nu = rad2deg(nu);
end

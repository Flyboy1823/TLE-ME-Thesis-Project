function [delta_r_lvlh_all, delta_v_lvlh_all] = eciToLvlhRelativeStates(r1_eci, v1_eci, r2_eci, v2_eci)
%ECITOLVLHRELATIVESTATES Transforms relative position and velocity from ECI to LVLH frame.
%
%   [delta_r_lvlh_all, delta_v_lvlh_all] = eciToLvlhRelativeStates(r1_eci, v1_eci, r2_eci, v2_eci)
%
%   Inputs:
%       r1_eci (Nx3 double): Position vectors of the Chief satellite in ECI frame (km).
%                            Each row is [x, y, z] at a specific time step.
%       v1_eci (Nx3 double): Velocity vectors of the Chief satellite in ECI frame (km/s).
%                            Each row is [vx, vy, vz] at a specific time step.
%       r2_eci (Nx3 double): Position vectors of the Follower satellite in ECI frame (km).
%                            Each row is [x, y, z] at a specific time step.
%       v2_eci (Nx3 double): Velocity vectors of the Follower satellite in ECI frame (km/s).
%                            Each row is [vx, vy, vz] at a specific time step.
%
%   Outputs:
%       delta_r_lvlh_all (Nx3 double): Relative position vectors of the Follower
%                                      with respect to the Chief in the LVLH frame (km).
%                                      Rows: X (Radial), Y (Along-track), Z (Cross-track).
%                                      (Output format remains 3xN for consistency with 3D vectors as columns)
%       delta_v_lvlh_all (Nx3 double): Relative velocity vectors of the Follower
%                                      with respect to the Chief in the LVLH frame (km/s).
%                                      Rows: Vx (Radial), Vy (Along-track), Vz (Cross-track).
%                                      (Output format remains 3xN for consistency with 3D vectors as columns)
%
%   The LVLH frame is defined with respect to the Chief satellite:
%   - X-axis (Radial): Points outward along the Chief's position vector.
%   - Y-axis (Along-track/Transverse): Points in the direction of the Chief's
%     velocity, completing a right-handed system (perpendicular to X and Z).
%   - Z-axis (Cross-track/Normal): Points perpendicular to the orbital plane,
%     opposite to the Chief's orbital angular momentum vector.

    % Ensure input dimensions are consistent (Nx3 format)
    if size(r1_eci, 2) ~= 3 || size(v1_eci, 2) ~= 3 || ...
       size(r2_eci, 2) ~= 3 || size(v2_eci, 2) ~= 3 || ...
       size(r1_eci, 1) ~= size(v1_eci, 1) || ...
       size(r1_eci, 1) ~= size(r2_eci, 1) || ...
       size(r1_eci, 1) ~= size(v2_eci, 1)
        error('Input matrices must be Nx3 and have the same number of rows (time steps).');
    end

    num_timesteps = size(r1_eci, 1); % Number of rows is now the number of time steps

    % Pre-allocate arrays for results (output format remains 3xN for 3D vectors)
    delta_r_lvlh_all = zeros(3, num_timesteps);
    delta_v_lvlh_all = zeros(3, num_timesteps);

    % Loop through each time step
    for k = 1:num_timesteps

        % Get current state vectors for chief and follower (transpose to column vectors for calculations)
        current_r1_eci = r1_eci(k, :)'; % Take k-th row and transpose to column
        current_v1_eci = v1_eci(k, :)'; % Take k-th row and transpose to column
        current_r2_eci = r2_eci(k, :)'; % Take k-th row and transpose to column
        current_v2_eci = v2_eci(k, :)'; % Take k-th row and transpose to column

        % --- Step 1: Calculate the Relative Position and Velocity Vectors in ECI ---
        % Position of Follower relative to Chief, expressed in ECI.
        delta_r_eci = current_r2_eci - current_r1_eci;
        % Velocity of Follower relative to Chief, expressed in ECI.
        delta_v_eci = current_v2_eci - current_v1_eci;

        % --- Step 2: Define the LVLH Frame Basis Vectors (relative to the Chief) ---

        % 2.1. Radial (R) axis - LVLH X-axis
        % Unit vector along the position vector of the chief satellite.
        x_lvlh_eci = current_r1_eci / norm(current_r1_eci);

        % 2.2. Normal (H) axis - LVLH Z-axis
        % Unit vector perpendicular to the orbital plane, opposite to the angular momentum.
        % Angular momentum vector of the chief: h = r x v
        h1_eci = cross(current_r1_eci, current_v1_eci);
        
        % Handle potential zero angular momentum (e.g., if r and v are collinear)
        if norm(h1_eci) < eps % Use a small epsilon to check for near-zero norm
            error('Chief satellite position and velocity vectors are nearly collinear, cannot define LVLH frame.');
        end
        z_lvlh_eci = -h1_eci / norm(h1_eci); % Negative because Z is typically anti-angular momentum

        % 2.3. Transverse (T) axis - LVLH Y-axis
        % Completes the right-handed system: Y = Z x X
        y_lvlh_eci = cross(z_lvlh_eci, x_lvlh_eci);

        % --- Step 3: Create the Transformation Matrix (ECI to LVLH) ---
        % The rotation matrix from LVLH to ECI has the LVLH basis vectors as its columns.
        % R_LVLH_to_ECI = [x_lvlh_eci, y_lvlh_eci, z_lvlh_eci]
        % The rotation matrix from ECI to LVLH is the transpose of R_LVLH_to_ECI.
        R_ECI_to_LVLH = [x_lvlh_eci';
                         y_lvlh_eci';
                         z_lvlh_eci'];

        % --- Step 4: Transform the Relative Position and Velocity Vectors to LVLH Frame ---
        delta_r_lvlh_all(k, :) = (R_ECI_to_LVLH * delta_r_eci)';
        delta_v_lvlh_all(k, :) = (R_ECI_to_LVLH * delta_v_eci)';

    end

end

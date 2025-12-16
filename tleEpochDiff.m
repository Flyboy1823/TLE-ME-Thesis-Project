function dt_seconds = tleEpochDiff(epoch1, epoch2)
    % tleEpochDiff - Calculate seconds between two TLE epoch times.
    %
    % Inputs:
    %   epoch1_str - First epoch string (e.g. '24153.54791435')
    %   epoch2_str - Second epoch string (e.g. '24153.64891435')
    %
    % Output:
    %   dt_seconds - Time difference in seconds (epoch2 - epoch1)
    %
    % Note: Assumes epochs are in format YYDDD.DDDDDDDD
    
    % Extract year and day of year
    year1 = 2000 + floor(epoch1 / 1000);
    year2 = 2000 + floor(epoch2 / 1000);
    
    doy1 = mod(epoch1, 1000);
    doy2 = mod(epoch2, 1000);
    
    % Convert to datetime objects
    dt1 = datetime(year1,1,1) + days(doy1 - 1);
    dt2 = datetime(year2,1,1) + days(doy2 - 1);
    
    % Calculate difference in seconds
    dt_seconds = seconds(dt2 - dt1);
end

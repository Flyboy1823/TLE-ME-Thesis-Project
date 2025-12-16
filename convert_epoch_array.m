function [start_datetime, end_datetime, t] = convert_epoch_array(epoch_array)
% convert_epoch_array: Converts an array of TLE epoch times into
% datetime values and seconds since first epoch.
%
% Input:
%   epoch_array - vector of TLE epoch values in YYDDD.DDDDDDDD format
%
% Outputs:
%   start_datetime - datetime of the first epoch
%   end_datetime   - datetime of the last epoch
%   t              - array of seconds since the first epoch

    % Ensure column vector
    epoch_array = epoch_array(:);

    % Extract year and DOY
    YY = floor(epoch_array / 1000);          % Two-digit year
    DOY = mod(epoch_array, 1000);            % Day of year with fraction

    % Expand YY to full 4-digit year
    year = zeros(size(YY));
    year(YY < 57) = 2000 + YY(YY < 57);      % 2000–2056
    year(YY >= 57) = 1900 + YY(YY >= 57);    % 1957–1999

    % Convert to datetime array
    epoch_datetime = datetime(year, 1, 1) + days(DOY - 1);

    % Sort in case input isn't ordered
    [epoch_datetime, sort_idx] = sort(epoch_datetime);

    % First and last timestamps
    start_datetime = epoch_datetime(1);
    end_datetime = epoch_datetime(end);

    % Time in seconds since first entry
    t = seconds(epoch_datetime - start_datetime);

    % Reorder t to match original input order
    t(sort_idx) = t;
end

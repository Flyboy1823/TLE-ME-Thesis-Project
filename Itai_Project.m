clc;
close all;
clear all

%% Calculating Orbital Parameters of Desired Satellites
% all angles in radians and distance in [Km]
[a,e,i,Om,om,n,M,Epoch,rp] = prs_TLE ('46807'); 

% Convert Epoch Times to seconds
[start_datetime, end_datetime, t] = convert_epoch_array(Epoch);

%% Calculate ECI coordinates

% Calcualte tru anomolys (degrees)
nu = rad2deg(mean2true_anomaly(M, e)); % degrees

% Calculate r and v cooridnates
r=[];
v=[];
for counter=1:length(t)
    [r(counter, :), v(counter, :)]=keplerian2ijk(a(counter), e(counter), i(counter), Om(counter), om(counter), nu(counter));
end

earth_sphere('km')
hold on

plot3(r(:, 1), r(:, 2), r(:, 3));
title('Oribt of Satelite Around Earth');
xlabel('x in ECI [km]');
ylabel('y in ECI [km]');
zlabel('z in ECI [km]');

grid on

hold off


%%  Estimation vs Observed

% Interpolations Orbit

time_tol=60;  %no need to interpolate if the time difference between the two observations is less than 1 minute

r_interoplated=[];
v_interpolated=[];
for counter=1:500%(length(t)-1)
    if (t(counter+1)-t(counter)) > time_tol
        interpolate_time=seconds(t(counter)+seconds(0:60:t(counter+1)));
        [temp_prop_positions,temp_prop_velocities] = propagateOrbit(interpolate_time,a(counter),e(counter),i(counter),Om(counter),om(counter),nu(counter));
        r_interoplated=[r_interoplated; temp_prop_positions(:,1:(end-1))'];
        v_interpolated=[v_interpolated; temp_prop_velocities(:,1:(end-1))'];
    else
        r_interoplated=[r_interoplated; r(counter,:)];
        v_interpolated=[v_interpolated; v(counter,:)];    
    end


end


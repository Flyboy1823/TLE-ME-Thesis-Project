clc;
close all;
clear all

%% Calculating Orbital Parameters of Desired Satellites
% all angles in radians and distance in [Km]
SatelliteID='46807';
[a,e,i,Om,om,n,M,Epoch,rp] = prs_TLE (SatelliteID); 
% a in km
% e without units
% i in rad
% n in 1/sec
% Epoch in datetime UTC
% rp in km

% Convert Angles to Degrees, semi major axis to m:
i=rad2deg(i);
Om = rad2deg(Om);
om = rad2deg(om);
M = rad2deg(M);
a = a*1000;

% Convert Epoch Times to seconds
[start_datetime, end_datetime, t] = convert_epoch_array(Epoch);
dates=start_datetime+seconds(t);  %time in datetime format
dates.TimeZone='UTC';

tleData_a = tleread('test.txt');




%% Calculate ECI coordinates

% Calcualte tru anomolys (degrees)
tru = mean2true_anomaly(M, e); % [deg]

% Calculate r and v cooridnates
r=[];
v=[];
for counter=1:length(t)
    [r(counter, :), v(counter, :)]=keplerian2ijk(a(counter), e(counter), i(counter), Om(counter), om(counter), tru(counter));
end

earth_sphere('km')
hold on

plot3(r(:, 1)/1000, r(:, 2)/1000, r(:, 3)/1000);
title('Oribt of Satelite Around Earth');
xlabel('x in ECI [km]');
ylabel('y in ECI [km]');
zlabel('z in ECI [km]');

grid on

hold off


%%  Estimation vs Observed

% Interpolations Orbit
time_tol=600;  % [sec] - no need to interpolate if the time difference between the two observations is less than 1 minute
r_interoplated=[];
v_interpolated=[];
t_interpolated=[];
timer=start_datetime;
fileLines = readlines(fullfile('TLE Histories', [SatelliteID,'.txt']));
fid = fopen(fullfile('TLE Histories', SatelliteID, 'temp.txt'), 'w');

t_tles = [];

end_run = 500;; %length(t);
for counter=2:end_run
    tleLine1_a = fileLines{(((counter-1)*2)-1)};
    tleLine2_a = fileLines{(((counter-1)*2))};
    tleToRead_a = {tleLine1_a; tleLine2_a};

    tleLine1_b = fileLines{(((counter)*2)-1)};
    tleLine2_b = fileLines{(((counter)*2))};
    tleToRead_b = {tleLine1_b; tleLine2_b};

    writelines(tleToRead_a, 'temp_a.txt');
    tleData_a = tleread('temp_a.txt');

    writelines(tleToRead_b, 'temp_b.txt');
    tleData_b = tleread('temp_b.txt');

    t_tles=[t_tles, tleData_a.Epoch];

    if counter == end_run
        t_tles=[t_tles, tleData_b.Epoch];
    end

    diff=t(counter)-t(counter-1);
    if (diff) > time_tol
        tleLine1_a = fileLines{(((counter-1)*2)-1)};
        tleLine2_a = fileLines{(((counter-1)*2))};
        tleToRead_a = {tleLine1_a; tleLine2_a};

        tleLine1_b = fileLines{(((counter)*2)-1)};
        tleLine2_b = fileLines{(((counter)*2))};
        tleToRead_b = {tleLine1_b; tleLine2_b};

        writelines(tleToRead_a, 'temp_a.txt')
        tleData_a = tleread('temp_a.txt');

        writelines(tleToRead_b, 'temp_b.txt')
        tleData_b = tleread('temp_b.txt');

        interpolate_time = [tleData_a.Epoch:seconds(time_tol):tleData_b.Epoch];

        if interpolate_time ~= tleData_b.Epoch
            interpolate_time= [interpolate_time, tleData_b.Epoch];
        end

        if isempty(t_interpolated)
            t_interpolated = interpolate_time;
        else
        t_interpolated=[t_interpolated(1:end-1), interpolate_time];
        end

        [temp_prop_positions,temp_prop_velocities] = propagateOrbit(interpolate_time,tleData_a);
        r_interoplated=[r_interoplated(1:end-1,:); temp_prop_positions'];
        v_interpolated=[v_interpolated(1:end-1,:); temp_prop_velocities'];
        % interpolate_time=[timer+seconds(0:time_tol:diff), dates(counter)];
        % tleLine1 = fileLines{(((counter-1)*2)-1)};
        % tleLine2 = fileLines{(((counter-1)*2))};
        % tleToRead = {tleLine1; tleLine2};
        % writelines(tleToRead, 'temp.txt')
        % tleData = tleread('temp.txt');
        % t_interpolated=[t_interpolated, interpolate_time(2:(end))];
        % [temp_prop_positions,temp_prop_velocities] = propagateOrbit(interpolate_time,tleData);
        % r_interoplated=[r_interoplated; temp_prop_positions(:,2:(end))'];
        % v_interpolated=[v_interpolated; temp_prop_velocities(:,2:(end))'];
        % timer=t_interpolated(end);
    else
        r_interoplated=[r_interoplated; r(counter,:)];
        v_interpolated=[v_interpolated; v(counter,:)]; 
        t_interpolated=[t_interpolated, dates(counter)];
    end

    %%%%% GET a tle's oribtal elemtns as new variables (a_tle=tle.a....)
    %%% Fix initial offset (why?)



end


a_interpolated=[a(1)];
e_interpolated=[e(1)];
i_interpolated=[i(1)];
Om_interpolated = [Om(1)];
om_interpolated = [om(1)];
tru_interpolated = [tru(1)];

%% Convert back to orbital elements
for counter=2:length(t_interpolated)
    [a_temp, e_temp, i_temp, Om_temp, om_temp, tru_temp] = ijk2keplerian(r_interoplated(counter,:), v_interpolated(counter,:)); 
    a_interpolated=[a_interpolated, a_temp];  %km
    e_interpolated=[e_interpolated, e_temp];    
    i_interpolated=[i_interpolated, i_temp];
    Om_interpolated = [Om_interpolated, Om_temp];
    om_interpolated = [om_interpolated, om_temp];
    tru_interpolated = [tru_interpolated, tru_temp];
end

% get indecies of interpolated array for tle data times
idx = ismember(t_interpolated, t_tles);

a_interpolated_filtered=a_interpolated(idx);
e_interpolated_filtered=e_interpolated(idx);
i_interpolated_filtered=i_interpolated(idx);
Om_interpolated_filtered=Om_interpolated(idx);
om_interpolated_filtered=om_interpolated(idx);
tru_interpolated_filtered=tru_interpolated(idx);



figure
plot(t_tles, a(1:length(t_tles)), 'b')
hold on
plot(t_tles, a_interpolated_filtered, 'r')
hold off
xlabel('$Epoch\ [UTC]$','interpreter','latex','FontSize',12')
ylabel('$Semi Majoro Axis\ [km]$','interpreter','latex','FontSize',12')
title('$Semi Major Axis\ Vs.\ Time$','interpreter','latex','FontSize',12')
subtitle(['Satellite ID - ', SatelliteID] ,'interpreter','latex','FontSize',9)
grid on

figure
plot(t_tles, e(1:length(t_tles)), 'b')
hold on
plot(t_tles, e_interpolated_filtered, 'r')
hold off
xlabel('$Epoch\ [UTC]$','interpreter','latex','FontSize',12')
ylabel('$Eccentricity\$','interpreter','latex','FontSize',12')
title('$Eccentricity\ Vs.\ Time$','interpreter','latex','FontSize',12')
subtitle(['Satellite ID - ', SatelliteID] ,'interpreter','latex','FontSize',9)
grid on

figure
plot(t_tles, i(1:length(t_tles)), 'b')
hold on
plot(t_tles, i_interpolated_filtered, 'r')
hold off
xlabel('$Epoch\ [UTC]$','interpreter','latex','FontSize',12')
ylabel('$i\ [deg]$','interpreter','latex','FontSize',12')
title('$Inclination\ Vs.\ Time$','interpreter','latex','FontSize',12')
subtitle(['Satellite ID - ', SatelliteID] ,'interpreter','latex','FontSize',9)
grid on


%% Steps:
% 1. Review polots for tolerences
% 2. set tolerances
% 3. set if conditions to pinpoint when diffreence between interpolated and
% and tle measurments are high enough then it is flagged as a potential
% change in orbit
% 4. Change a, e, i, Om, om, tru to a_tle, e_tle....
% 5. Smooth over the plots for a_interpolated, e_interpolated, i_interpolated

index_potential_a_change=[];
index_potential_i_change=[];
index_potential_e_change=[];

% % for counter=1:length(t)
% % 	%%%%% if a_interporeted_filtered(counter) < a(counter) ....
% % end
% 
% 





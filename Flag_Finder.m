function [flag, time_tle, r_output, v_output] = Flag_Finder(SatelliteID)


    %% Calculating Orbital Parameters of Desired Satellites
    if isnumeric(SatelliteID)
        SatelliteID= num2str(SatelliteID);
    end
    currentFolder = fileparts(mfilename('fullpath')); % Get the directory of the running script
    filepath = fullfile('TLE Histories', [SatelliteID, '.txt']);  % Full file name
    fileContent = readlines(filepath);
    nonEmptyLines = sum((strlength(fileContent) > 0));
    tle_len = nonEmptyLines/2;
    
    time_tol=minutes(1);  % [duration] - no need to interpolate if the time difference between the two observations is less than 1 minute
    a_tol=3000; % [meters]
    i_tol=2; %[deg]
    e_tol=.05;
    a_tle=[];
    e_tle=[];
    i_tle=[];
    RAAN_tle=[];
    om_tle=[];
    M_tle=[];
    tru_tle=[];
    time_tle=[];
    r_tle=[];
    v_tle=[];
    r_interoplated=[];
    v_interpolated=[];
    time_interpolated=[];
    a_interpolated_filtered=[];
    e_interpolated_filtered=[];
    i_interpolated_filtered=[];
    flag={};
    
    MU = 3.986004418e14; %[m^3/s^2]
    
    for counter=2:round(tle_len)  
        tleLine1_a = fileContent{(((counter-1)*2)-1)};
        tleLine2_a = fileContent{(((counter-1)*2))};
        tleToRead_a = {tleLine1_a; tleLine2_a};
        writelines(tleToRead_a, 'temp_a.txt');
        tleLine1_b = fileContent{(((counter)*2)-1)};
        tleLine2_b = fileContent{(((counter)*2))};
        tleToRead_b = {tleLine1_b; tleLine2_b};
        writelines(tleToRead_b, 'temp_b.txt');
        
        tleData_a = tleread('temp_a.txt');
    
        n = (((tleData_a.MeanMotion)/0.0042)./(24*60*60))*2*pi; %mean motion [1/sec]
        a = (MU./(n.^2)).^(1/3); %semi major axis [m]
        e = tleData_a.Eccentricity; %
        M = tleData_a.MeanAnomaly; % [deg]
        i = tleData_a.Inclination; % [deg]
        RAAN = tleData_a.RightAscensionOfAscendingNode; %[deg]
        om = tleData_a.ArgumentOfPeriapsis; %[deg]
        tru = mean2true_anomaly(M, e); % [deg]
        Epoch_a = tleData_a.Epoch; %[UTC]
        [r, v]=keplerian2ijk(a, e, i, RAAN, om, tru);
    
        tleData_b = tleread('temp_b.txt');
    
        n_b = (((tleData_b.MeanMotion)/0.0042)./(24*60*60))*2*pi; %mean motion [1/sec]
        a_b = (MU./(n_b.^2)).^(1/3); %semi major axis [m]
        e_b = tleData_b.Eccentricity; %
        i_b = tleData_b.Inclination; % [deg]
        Epoch_b = tleData_b.Epoch; %[UTC]
    
        a_tle=[a_tle, a];
        e_tle=[e_tle, e];
        i_tle=[i_tle, i];
        RAAN_tle=[RAAN_tle, RAAN];
        om_tle=[om_tle, om];
        M_tle=[M_tle, M];
        tru_tle=[tru_tle, tru];
        r_tle=[r_tle, r];
        %v_tle=[v_tle, v];
        time_tle=[time_tle, Epoch_a];
    
        if counter == 2
            a_interpolated_filtered=[a];
            i_interpolated_filtered=[i];
            e_interpolated_filtered=[e];
        end
    
        if (Epoch_b-Epoch_a) > time_tol
            interpolate_time = [Epoch_a:time_tol:Epoch_b];
    
            if interpolate_time(end) ~= tleData_b.Epoch
                interpolate_time= [interpolate_time, Epoch_b];
            end
    
            time_interpolated=[time_interpolated, interpolate_time];
            [temp_prop_positions,temp_prop_velocities] = OrbitProp(r,v,Epoch_a, Epoch_b);
            %[temp_prop_positions, temp_prop_velocities] = Orbit_Prop_Test(r'/1000,v'/1000,seconds(Epoch_b-Epoch_a));
            
            r_interoplated=[r_interoplated, temp_prop_positions(:, 1:(end-1))];
            %v_interpolated=[v_interpolated, temp_prop_velocities(:, 1:(end-1))];
            
            [a_b_inter, e_b_inter, i_b_inter, RAAN_b_inter, om_b_inter, tru_b_inter] = ijk2keplerian(temp_prop_positions(:,end), temp_prop_velocities(:,end)); 
    
            a_interpolated_filtered=[a_interpolated_filtered, a_b_inter];
            i_interpolated_filtered=[i_interpolated_filtered, i_b_inter];
            e_interpolated_filtered=[e_interpolated_filtered, e_b_inter];
            
            % Flag potential orbit change indecies in the tle array
            if (abs(a_b_inter-a_b)) > a_tol
                flag(:, end+1)={'Semi-Major Axis tolarance exceeded'; Epoch_b; ['Difference of ' sprintf('%.1f', (a_b_inter - a_b)/1000) ' KM']; ['Tolerence is ' sprintf('%.1f', a_tol/1000) ' KM']};
            end
            
            if (abs(e_b_inter-e_b)) > e_tol
                flag(:, end+1)={'Eccentricity tolarance exceeded'; Epoch_b; ['Difference of ' sprintf('%.2f', (e_b_inter - e_b))]; ['Tolerence is ' sprintf('%.1f', e_tol)]};
            end
             
            if (abs(i_b_inter-i_b)) > i_tol
                flag(:, end+1)={'Inclination tolarance exceeded'; Epoch_b; ['Difference of ' sprintf('%.2f',(i_b_inter - i_b)) ' Degrees']; ['Tolerence is ' sprintf('%.1f', i_tol) ' KM']};
            end
    
        else
            time_interpolated=[time_interpolated, Epoch_b];
            r_interoplated=[r_interoplated, r];
            %v_interpolated=[v_interpolated, v];
            a_interpolated_filtered=[a_interpolated_filtered, a_b];
            e_interpolated_filtered=[e_interpolated_filtered, e_b];
            i_interpolated_filtered=[i_interpolated_filtered, i_b];
        end
    end
    
    a_tle=[a_tle, a_b];
    e_tle=[e_tle, e_b];
    i_tle=[i_tle, i_b];
    RAAN_tle=[RAAN_tle, tleData_b.RightAscensionOfAscendingNode];
    om_tle=[om_tle, tleData_b.ArgumentOfPeriapsis];
    M_tle=[M_tle, tleData_b.MeanAnomaly];
    tru_tle=[tru_tle, mean2true_anomaly(M(end), e(end))];
    [r_end, v_end]=keplerian2ijk(a_tle(end), e_tle(end), i_tle(end), RAAN_tle(end), om_tle(end), tru_tle(end));
    r_tle=[r_tle, r_end];
    v_tle=[v_tle, v_end];
    time_tle=[time_tle, Epoch_b];
    
    
    %% Plots
    
    earth_sphere('km')
    hold on
    plot3(r_tle(1, :)/1000, r_tle(2, :)/1000, r_tle(3, :)/1000);
    title('Oribt of Satelite Around Earth');
    xlabel('x in ECI [km]');
    ylabel('y in ECI [km]');
    zlabel('z in ECI [km]');
    grid on
    hold off
    
    figure
    plot(time_tle, a_tle/1000, 'b')
    hold on
    plot(time_tle, a_interpolated_filtered/1000, 'r')
    hold off
    xlabel('Epoch [UTC]')
    ylabel('Semi Major Axis [km]')
    title('Semi Major Axis Vs. Time')
    legend('TLE measurments', 'Interpolations')
    subtitle(['Satellite ID - ', SatelliteID] ,'interpreter','latex','FontSize',9)
    grid on
    
    figure
    plot(time_tle, e_tle, 'b')
    hold on
    plot(time_tle, e_interpolated_filtered, 'r')
    hold off
    xlabel('Epoch [UTC]$')
    ylabel('Eccentricity')
    title('Eccentricity Vs. Time')
    legend('TLE measurments', 'Interpolations')
    subtitle(['Satellite ID - ', SatelliteID] ,'FontSize',9)
    grid on
    
    figure
    plot(time_tle, i_tle, 'b')
    hold on
    plot(time_tle, i_interpolated_filtered, 'r')
    hold off
    xlabel('Epoch [UTC]')
    ylabel('i [deg]')
    legend('TLE measurments', 'Interpolations');
    title('Inclination Vs. Time')
    subtitle(['Satellite ID - ', SatelliteID] ,'FontSize',9)
    grid on

    r_output = r_tle';
    v_output = v_tle';
end

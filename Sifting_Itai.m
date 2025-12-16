%% Sifting through TLE file containing different satellites in every row. 
%% Input: Range of inclination, eccentricity and altitude.
%% Output: The matching satellite catalogs. 

clc
clear all
close all

% arranging the TLE data in three different vectors:
% 1. Satellite Name
% 2. First line of TLE:
% - Line Number
% - Satellite Catalog and Classification
% - International Designator
% - Epoch
% - 1st Derivative of mean motion
% - 2nd Derivative of mean motion
% - BSTAR Drag Term,Ephemeris type
% - Element number
% - Checksum 
% 3. Second line of TLE:
% - Line Number
% - Satellite Catalog
% - Inclination, i [deg]
% - RAAN, Omega [deg]
% - eccentricity, e
% - Argument of Perigee, omega [deg]
% - Mean Anomaly, M [deg]
% - Mean Motion, n [rev/day]

%constants
mu = 3.986e5; %[km^3/s^2]
Re = 6378; %[km]

TLE = "celes_tle.json"; %this file needs to be downlaoded as json using a script in python
txt = fileread(TLE);
data = jsondecode(txt); 

%initializing vectors
Sat_Name = string(NaN(length(data),1));
TLE_1 = string(NaN(length(data),1));
TLE_2 =  string(NaN(length(data),1));


%converting the TLE data to vector of type double 
for i=1:length(data)
    Sat_Name(i,:) =  data(i, 1).OBJECT_NAME;
    %Sat_Name(i,:) =  data{i, 1}.satellite_name;
    TLE_1(i,:) = data(i, 1).TLE_LINE1;
    TLE_2(i,:) = data(i, 1).TLE_LINE2;
    temp_string = strsplit(TLE_2(i,:));
    split_temp(i,:) = temp_string(1:8);
    split_temp(i,5) = "0." + split_temp(i,5);
    TLE2_param(i,:) = str2double(split_temp(i,:)); 
    n(i) = (TLE2_param(i,8))*2*pi/(24*60*60); %mean motion [1/sec]
    e(i) = TLE2_param(i,5);
    a(i) = (mu/(n(i)^2))^(1/3); %semi major axis [Km]
    h(i) = a(i)- Re; %altitude [Km]
    rp(i) = a(i)*(1-e(i)); %[Km]
end

i_min = input("Please enter min value of the inclination (in degrees) : ");
i_max = input("Please enter max value of the inclination (in degrees) : ");

e_min = input("Please enter min value of the eccentricity: ");
e_max = input("Please enter max value of the eccentricity: ");

h_min = input("Please enter min value of the altitude (in Km) : ");
h_max = input("Please enter max value of the altitude (in Km) : ");

flag = zeros(length(data),1);

% trimmed_TLE2 = TLE2_param(TLE2_param>i_min & TLE2_param<i_max & TLE2_param<e_min TLE2>e_max)
for i=1:length(data)
    if TLE2_param(i,5)>e_min && TLE2_param(i,5)<e_max && TLE2_param(i,3)>i_min && TLE2_param(i,3)<i_max && h_min<h(i) && h(i)<h_max
        flag(i) = 1;
    end
end

TLE2_param_sifted = TLE2_param(flag==1,:);
Satellite_Name = Sat_Name(flag==1);
Satellite_Catalog = TLE2_param_sifted(:,2);
Sifted_Sats = [Satellite_Name,Satellite_Catalog];

if isempty(Satellite_Catalog)
    disp('No satellites meet the specifies input conditions')
else
    fprintf('There are %d satellites that match the given parameters and they are:\n ',length(Satellite_Catalog));
    disp(Sifted_Sats)
end

%% Save the search results
% Create the folder "Sat Searches" if it doesn't exist
folderName = 'Sat Searches';
if ~exist(folderName, 'dir')
    mkdir(folderName); % Create the folder if it doesn't exist
end

% Create the file name based on the user inputs
fileName = sprintf('satellite_output_Inc%.2f-%.2f_Ecc%.2f-%.2f_Alt%.2f-%.2f.txt', ...
    i_min, i_max, e_min, e_max, h_min, h_max);

% Construct the full path for the file in the "Sat Searches" folder
filePath = fullfile(folderName, fileName);

% Open the text file for writing
fileID = fopen(filePath, 'w');

% Write the user input data to the text file
fprintf(fileID, 'User Inputs:\n');
fprintf(fileID, 'Inclination Range: %.2f to %.2f degrees\n', i_min, i_max);
fprintf(fileID, 'Eccentricity Range: %.2f to %.2f\n', e_min, e_max);
fprintf(fileID, 'Altitude Range: %.2f to %.2f km\n', h_min, h_max);
fprintf(fileID, '\n');

% Check if there are any satellites in the filtered list
if isempty(Sifted_Sats)
    fprintf(fileID, 'No satellites match the given parameters.\n');
else
    % Write the filtered list of satellites to the text file
    fprintf(fileID, 'Filtered Satellites:\n');
    for i = 1:size(Sifted_Sats, 1)
        fprintf(fileID, 'Satellite: %s, Catalog: %d\n', Sifted_Sats{i, 1}, Sifted_Sats{i, 2});
    end
end

% Close the file
fclose(fileID);

% Notify user that the file has been created
disp(['Output has been saved to ', filePath]);


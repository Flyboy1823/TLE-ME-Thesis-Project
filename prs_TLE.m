function [a,e,i,Om,om,n,M,Epoch_Vec,rp] = prs_TLE (TLE_Designator)

%% Input: TLE txt file designator number as string. For example, '26360'
%% Output: Orbital Elements and epoch vector

currentFolder = fileparts(mfilename('fullpath')); % Get the directory of the running script
file = fullfile(currentFolder, 'TLE Histories', TLE_Designator);  % Full file name


%file = strcat(TLE_Designator,'.txt');

mu = 3.986e5; %[km^3/s^2]
Re = 6378; %[km]

%% Importing the TLE txt file, dividing it into two rows, and labeling the variables
TLE = readtable(file,'Delimiter',' ', 'MultipleDelimsAsOne', true, 'Format','%d%s%s%f%s%s%s%f%f');

TLE1 = TLE(1:2:end,:);
TLE_1 = renamevars(TLE1,["Var2","Var3", "Var4", "Var5", "Var6", "Var7","Var8","Var9"],...
    ["Satellite_Number", "International_Designator","Epoch", "1st_Derivative_of_Mean_Motion","2nd_Derivative_of_Mean_Motion","Drag_Term_or_Radiation_Pressure_Coefficient","Ephemeris_Type", "Element_Number&Check_Sum"]);
TLE2 = TLE(2:2:end,:);
TLE_2 = renamevars(TLE2,["Var2","Var3", "Var4", "Var5", "Var6", "Var7","Var8","Var9"],...
    ["Satellite_Number", "Inclination","RAAN", "Eccentricity","Argument_of_Perigee","Mean_Anomaly",...
    "Mean_Motion", "Revolution#_at_Epoch&Check_Sum"]);

%% Fixing the Eccentricity value and converting all relevent values to be double/string
TLE_2 = convertvars(TLE_2,{'Satellite_Number','Inclination','Eccentricity','Argument_of_Perigee','Mean_Anomaly'},'string');
for i=1:height(TLE_2)
    TLE_2{i,5} = strcat('0.',TLE_2{i,5});
end
TLE_2 = convertvars(TLE_2,{'Satellite_Number','Inclination','RAAN','Eccentricity','Argument_of_Perigee','Mean_Anomaly',...
    'Mean_Motion','Revolution#_at_Epoch&Check_Sum'}, 'double');

TLE_1 = convertvars(TLE_1,{'Satellite_Number', 'International_Designator','Epoch', '1st_Derivative_of_Mean_Motion',...
    '2nd_Derivative_of_Mean_Motion','Drag_Term_or_Radiation_Pressure_Coefficient','Ephemeris_Type', 'Element_Number&Check_Sum'},'string');
for i=1:height(TLE_1)
    TLE_1{i,5} = strcat('0',TLE_1{i,5});
end
TLE_1 = convertvars(TLE_1,{'Epoch', '1st_Derivative_of_Mean_Motion','2nd_Derivative_of_Mean_Motion',...
    'Drag_Term_or_Radiation_Pressure_Coefficient','Ephemeris_Type', 'Element_Number&Check_Sum'},'double');

n = ((TLE_2.Mean_Motion)./(24*60*60))*2*pi; %mean motion [1/sec]
a = (mu./(n.^2)).^(1/3); %semi major axis [Km]
%     h(i) = a(i)- Re; %altitude [Km]
rp = a.*(1-TLE_2.Eccentricity); %[Km]
M = deg2rad(TLE_2.Mean_Anomaly); %[rad]
i = deg2rad(TLE_2.Inclination); %[rad]
Om = deg2rad(TLE_2.RAAN); %[rad]
om = deg2rad(TLE_2.Argument_of_Perigee); %[rad]
e = TLE_2.Eccentricity;
Epoch_Vec = TLE_1.Epoch; %[UTC]
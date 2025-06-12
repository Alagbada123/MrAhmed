% create_experimental_data.m
% This script creates experimental data structure for concrete samples
% Data includes Oven-Dried (OD), Saturated Surface-Dried (SSD) measurements
% for Young's Modulus (E) and Compressive Strength (CS).
% Measurements aligned with ages: 2, 7, 14, 28, 56, and 91 days

clear experimental_data; % Clear existing variable if script is re-run
experimental_data = struct();

% --- Common Data ---
% Hydration degrees corresponding to 2, 7, 14, 28, 56, 91 days respectively
experimental_data.xi_marker_vals = [0.60, 0.70, 0.80, 0.86, 0.90, 0.92];
experimental_data.days = [2, 7, 14, 28, 56, 91];
% Sample descriptions used consistently
experimental_data.sample_descriptions = {
    'NAC     - Natural Aggregate Concrete', ...
    'RPB     - Recycled Perforated Brick Aggregate', ...
    'RSB     - Recycled Solid Brick Aggregate', ...
    'RC      - Recycled Concrete Aggregate', ...
    'RCa     - Recycled Abraded Concrete Aggregate', ...
    'RCc     - Recycled Carbonated Concrete Aggregate', ...
    'RCh     - Recycled Hydrophobic Concrete Aggregate'
};

%% Oven-Dried (OD) Samples Data
experimental_data.OD = struct();
experimental_data.OD.condition = 'Oven-Dried';

% OD Young's modulus measurements [GPa]
% Rows: Samples in order of experimental_data.sample_descriptions
% Columns: Measurements at 2, 7, 14, 28, 56, 91 days
experimental_data.OD.youngs_modulus = [
    34.0  40.1  42.0  43.1  45.8  46.4;  % NAC
    22.2  26.2  26.9  28.7  29.9  30.1;  % RPB-OD
    26.2  30.5  31.8  32.9  33.4  34.5;  % RSB-OD
    29.7  31.7  35.9  38.8  40.5  40.2;  % RC-OD (Value for day 7 corrected from previous version based on image)
    29.8  35.5  36.4  38.3  38.5  40.4;  % RCa-OD
    30.5  35.5  37.2  37.9  39.4  40.6;  % RCc-OD
    29.8  34.5  35.7  36.7  37.7  39.7   % RCh-OD
];
% OD Standard deviation for Young's Modulus [GPa]
experimental_data.OD.youngs_modulus_stddev = [
    0.6  0.7  1.3  0.8  0.7  1.1;  % NAC
    0.3  0.5  0.5  0.2  0.2  0.5;  % RPB-OD
    0.8  0.3  0.2  0.2  0.2  0.2;  % RSB-OD
    0.1  1.1  1.0  1.3  0.5  0.2;  % RC-OD
    0.8  0.3  1.8  0.3  0.6  0.7;  % RCa-OD
    0.4  1.9  1.6  0.2  0.1  0.8;  % RCc-OD
    0.9  1.0  1.2  1.1  1.0  1.0   % RCh-OD
];

% OD Compressive strength measurements [MPa] --- NEW ---
% Columns: Measurements at 2, 7, 14, 28, 56, 91 days (Using NaN for missing days 2, 14)
experimental_data.OD.strength_MPa = [
%    2     7     14    28    56    91  days
    NaN  38.6   NaN  47.2  49.8  53.9;  % NAC
    NaN  25.7   NaN  30.5  35.9  32.9;  % RPB-OD
    NaN  32.3   NaN  38.0  41.6  42.5;  % RSB-OD
    NaN  39.2   NaN  44.5  48.4  50.1;  % RC-OD
    NaN  34.5   NaN  45.1  50.5  48.7;  % RCa-OD
    NaN  36.3   NaN  46.7  48.8  49.9;  % RCc-OD
    NaN  35.8   NaN  43.1  47.4  46.6   % RCh-OD
];
% OD Standard deviation for Compressive Strength [MPa] --- NEW ---
% Columns: Measurements at 2, 7, 14, 28, 56, 91 days (Using NaN for missing days 2, 14)
experimental_data.OD.strength_MPa_stddev = [
%    2     7     14    28    56    91  days
    NaN  0.4    NaN   0.3   0.8   1.9;  % NAC
    NaN  0.0    NaN   1.1   1.9   0.9;  % RPB-OD
    NaN  1.4    NaN   1.7   1.2   0.8;  % RSB-OD
    NaN  1.7    NaN   0.4   0.4   0.5;  % RC-OD
    NaN  0.3    NaN   1.0   0.3   0.3;  % RCa-OD
    NaN  0.9    NaN   0.2   0.2   2.0;  % RCc-OD
    NaN  0.4    NaN   1.1   1.2   0.7   % RCh-OD
];


%% Saturated Surface-Dried (SSD) Samples Data
experimental_data.SSD = struct();
experimental_data.SSD.condition = 'Saturated Surface-Dried';

% SSD Young's modulus measurements [GPa]
% Rows: Samples in order of experimental_data.sample_descriptions
% Columns: Measurements at 2, 7, 14, 28, 56, 91 days
experimental_data.SSD.youngs_modulus = [
    34.0  40.1  42.0  43.1  45.8  46.4;  % NAC
    21.3  25.6  27.4  28.4  29.3  30.5;  % RPB-SSD
    25.8  31.1  32.4  33.8  34.7  35.3;  % RSB-SSD
    30.3  35.3  37.2  38.3  39.3  40.3;  % RC-SSD
    31.3  36.4  38.4  38.9  40.0  41.9;  % RCa-SSD
    30.7  35.1  37.1  38.1  39.3  39.4;  % RCc-SSD
    29.3  34.8  36.7  36.8  39.6  39.9   % RCh-SSD
];
% SSD Standard deviation for Young's Modulus [GPa]
experimental_data.SSD.youngs_modulus_stddev = [
    0.6  0.7  1.3  0.8  0.7  1.1;  % NAC
    1.4  1.2  1.4  1.5  1.6  1.6;  % RPB-SSD
    0.4  0.5  0.8  1.0  0.9  0.5;  % RSB-SSD
    0.9  1.2  1.6  1.0  0.5  0.8;  % RC-SSD
    0.5  1.2  0.8  0.8  0.6  2.1;  % RCa-SSD
    0.6  0.6  0.5  0.6  0.8  0.3;  % RCc-SSD
    0.3  1.0  0.8  0.8  0.9  0.2   % RCh-SSD
];

% SSD Compressive strength measurements [MPa] --- UPDATED ---
% Columns: Measurements at 2, 7, 14, 28, 56, 91 days (Using NaN for missing days 2, 14)
experimental_data.SSD.strength_MPa = [
%    2     7     14    28    56    91  days
    NaN  38.6   NaN  47.2  49.8  53.9;  % NAC
    NaN  26.4   NaN  30.9  32.7  34.1;  % RPB-SSD
    NaN  30.9   NaN  38.7  41.8  42.4;  % RSB-SSD
    NaN  35.2   NaN  44.2  46.0  48.8;  % RC-SSD
    NaN  35.2   NaN  46.6  50.4  49.5;  % RCa-SSD
    NaN  35.6   NaN  45.6  49.9  49.4;  % RCc-SSD
    NaN  35.6   NaN  42.4  45.2  47.6   % RCh-SSD
];
% SSD Standard deviation for Compressive Strength [MPa] --- UPDATED ---
% Columns: Measurements at 2, 7, 14, 28, 56, 91 days (Using NaN for missing days 2, 14)
experimental_data.SSD.strength_MPa_stddev = [
%    2     7     14    28    56    91  days
    NaN  0.4    NaN   0.3   0.8   1.9;  % NAC
    NaN  2.7    NaN   1.1   2.1   1.5;  % RPB-SSD
    NaN  0.7    NaN   0.9   0.9   0.6;  % RSB-SSD
    NaN  1.2    NaN   0.5   1.0   0.2;  % RC-SSD
    NaN  0.4    NaN   0.9   1.3   0.5;  % RCa-SSD
    NaN  0.5    NaN   0.1   1.3   1.1;  % RCc-SSD
    NaN  0.2    NaN   0.8   0.4   0.2   % RCh-SSD
];

% Display message
disp('Experimental data structure created/updated.');
disp('Includes Young''s Modulus and Compressive Strength for OD and SSD conditions.');
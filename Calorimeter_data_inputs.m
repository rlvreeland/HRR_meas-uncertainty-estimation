clc, clear, close all

% This codes uses the user specified fuel properties, calorimeter 
% specifications and an array of HRR to create model calorimeter data that 
% can be used to estimate the uncertainty for a proposed calorimeter
% design. This code is specifically for a paragamatic O2 sensor and accounts
% for the effects of CO2 on the measurements.

% Define HRR range 
spacing = 150;
HRR_kw = linspace(75, 750, spacing).';
% HRR_kw = 75;

%% Propane Burner
%%%%% Fire Model
% Specify the fuel properties in the model
Propane_burner = Fire_model;
Propane_burner.Mw_fuel_g_per_mol=44.09;
Propane_burner.delta_hc_kJ_per_kg=46.0*10^3;
Propane_burner.yield_CO_g_per_g_fuel=0.005;
Propane_burner.yield_CO2_g_per_g_fuel=2.85;

% Add in the calorimeter specifications
Propane_burner.D_duct_m = 0.450;
Propane_burner.h_hood_m = 7*12*2.54/100;
Propane_burner.flow_fraction = 1.2;

% For this square burner the the equivalent circular fire diameter by area was
% used
L_burner_m = 1.0;
Propane_burner = Burner_diameter(Propane_burner, L_burner_m);

% Create model data given input model parameters for the temperature,
% differential pressure, and concentrations of CO, CO2, and O2 in the duct.
Propane_burner = model_HRR_measurements(Propane_burner, HRR_kw);

%%%%% Uncertainty
% Specify the uncertainty in each measurement
Propane_burner_uncert = measurement_uncert;

% Uncertainties in flow measurement device, assuming no change in ambient
% pressure
Propane_burner_uncert.Delta_P_Pa = 0.0015*Propane_burner.Delta_P_Pa; 
       % from pressure transducer data sheet
% add in uncertainty in duct temperature measurement given that TC is a
% regular k-type TC
Propane_burner_uncert = uncert_Te_K_k_type(Propane_burner_uncert, ...
    Propane_burner.Te_K);

% Uncertainties in duct specification
Propane_burner_uncert.D_duct_m = 0.0005;

% Calibration bottles with uncertainties
Propane_burner.x_CO_bottle_span = 4/100;
Propane_burner.x_CO2_bottle_span = 8/100;
Propane_burner.x_O2_bottle_span = 21/100;
% Propane_burner.x_O2_bottle_span = 100/100;
Propane_burner.x_i_bottle_zero = 0; %N2 is used for zero, in this 
% case 99.99% N2 was used. 

Propane_burner_uncert.x_CO_bottle_span = 0.02*Propane_burner.x_CO_bottle_span...
    /(3^0.5); % convert to standard error
Propane_burner_uncert.x_CO2_bottle_span = 0.02*Propane_burner.x_CO2_bottle_span...
    /(3^0.5); % convert to standard error
Propane_burner_uncert.x_O2_bottle_span = 0.02*Propane_burner.x_O2_bottle_span...
    /(3^0.5); % convert to standard error
%Propane_burner_uncert.x_O2_bottle_span = 0.1/100/(3^0.5); 
              % convert to standard error
% Calculate the uncertainty for the zero given nitrogen concentration
x_N2_bottle = 99.9/100;
Propane_burner_uncert.x_i_bottle_zero = (1-x_N2_bottle)/(3^0.5); 
            % converts to standard error

% Run functions to calculate the uncertainty in concentration and flow 
% measurements. 
% Please note that these will not be used for HRR are just for reference.
% Instead calculation will be repeated when calculating the HRR.
Propane_burner_uncert = uncert_conc_measure(Propane_burner_uncert, ...
    Propane_burner);
Propane_burner_uncert = uncert_mass_flow_measure(Propane_burner_uncert, ...
    Propane_burner);

% Calculate HRR and uncertainty given model data not considering CO and CO2
[Propane_burner_uncert_without, Propane_burner_HRR_calc_kW_without, ~] = ...
    uncert_HRR(Propane_burner_uncert, Propane_burner, 'False');

% Calculate HRR and uncertainty given model data considering CO and CO2
[Propane_burner_uncert_with, Propane_burner_HRR_calc_kW_with, ...
    Propane_burner_uncert_frac] = uncert_HRR(Propane_burner_uncert, ...
    Propane_burner, 'True');



%% Heptane Pool Fire
%%%%% Fire Model
% Specify the fuel properties in the model
Heptane_pool = Fire_model;
Heptane_pool.Mw_fuel_g_per_mol=100.20;
Heptane_pool.delta_hc_kJ_per_kg=44.6*10^3;
Heptane_pool.yield_CO_g_per_g_fuel=0.010;
Heptane_pool.yield_CO2_g_per_g_fuel=2.85;

% Add in the calorimeter specifications
Heptane_pool.D_duct_m = 0.450;
Heptane_pool.h_hood_m = 7*12*2.54/100;
Heptane_pool.flow_fraction = 1.2;

% Calculate the diameter for a circular pool fire to obtain the given HRR
% using correlation from Zabetakis and Burgess (1961).
m_dot_inf_kg_per_m2 = 0.101;
k_beta = 1.1;
Heptane_pool = Pool_fire_diameter(Heptane_pool, HRR_kw, m_dot_inf_kg_per_m2,...
    k_beta);

% Create model data given input model parameters for the temperature,
% differential pressure, and concentrations of CO, CO2, and O2 in the duct.
Heptane_pool = model_HRR_measurements(Heptane_pool, HRR_kw);


%%%%% Uncertainty
% Specify the uncertainty in each measurement
Heptane_pool_uncert = measurement_uncert;

% Uncertainties in flow measurement device, assuming no change in ambient
% pressure
Heptane_pool_uncert.Delta_P_Pa = 0.0015*Heptane_pool.Delta_P_Pa; 
       % from pressure transducer data sheet
% add in uncertainty in duct temperature measurement given that TC is a
% regular k-type TC
Heptane_pool_uncert = uncert_Te_K_k_type(Heptane_pool_uncert, ...
    Heptane_pool.Te_K);

% Uncertainties in duct specification
Heptane_pool_uncert.D_duct_m = 0.0005;

% Calibration bottles with uncertainties
Heptane_pool.x_CO_bottle_span = 5/100;
Heptane_pool.x_CO2_bottle_span = 10/100;
Heptane_pool.x_O2_bottle_span = 21/100;
Heptane_pool.x_i_bottle_zero = 0; %N2 is used for zero, in this 
% case 99.9% N2 was used. 

Heptane_pool_uncert.x_CO_bottle_span = 0.02*Heptane_pool.x_CO_bottle_span...
    /(3^0.5); % convert to standard error
Heptane_pool_uncert.x_CO2_bottle_span = 0.02*Heptane_pool.x_CO2_bottle_span...
    /(3^0.5); % convert to standard error
Heptane_pool_uncert.x_O2_bottle_span = 0.02*Heptane_pool.x_O2_bottle_span...
    /(3^0.5); % convert to standard error
% Calculate the uncertainty for the zero given nitrogen concentration
x_N2_bottle = 99.9/100;
Heptane_pool_uncert.x_i_bottle_zero = (1-x_N2_bottle)/(3^0.5); 
            % converts to standard error

% Run functions to calculate the uncertainty in concentration and flow 
% measurements. 
% Please note that these will not be used for HRR are just for reference.
% Instead calculation will be repeated when calculating the HRR.
Heptane_pool_uncert = uncert_conc_measure(Heptane_pool_uncert, ...
    Heptane_pool);
Heptane_pool_uncert = uncert_mass_flow_measure(Heptane_pool_uncert, ...
    Heptane_pool);


% Calculate HRR and uncertainty given model data not considering CO and CO2
[Heptane_pool_uncert_without, Heptane_pool_HRR_calc_kW_without, ~] = ...
    uncert_HRR(Heptane_pool_uncert, Heptane_pool, 'False');

% Calculate HRR and uncertainty given model data considering CO and CO2
[Heptane_pool_uncert_with, Heptane_pool_HRR_calc_kW_with, ...
    Heptane_pool_uncert_frac] = uncert_HRR(Heptane_pool_uncert, ...
    Heptane_pool, 'True');


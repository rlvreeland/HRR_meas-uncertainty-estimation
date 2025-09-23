clc, clear, close all

% This codes uses the user specfied fuel properties, calorimeter 
% specifications and an array of HRR to create model calorimeter data that 
% can be used to estimate the uncertainty for a proposed calorimeter
% design. This code is specfically for a paragamatic O2 sensor and accounts
% for the effects of CO and CO2 on the measurmements

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

% For this square burner the the equivent circular fire diamter by area was
% used
L_burner_m = 1.0;
Propane_burner = Burner_diameter(Propane_burner, L_burner_m);

% Create model data given input model parameters for the temperature,
% differential pressure, and concentrations of CO, CO2, and O2 in the duct.
Propane_burner = model_HRR_measurements(Propane_burner, HRR_kw);

%%%%% Uncertainty
% Specify the uncertainty in each measurement
Propane_burner_uncert = measurement_uncert;

% Uncertainities in flow measurement device, assuming no change in ambient
% pressure
Propane_burner_uncert.Delta_P_Pa = 0.0015*Propane_burner.Delta_P_Pa; 
       % from pressure transducter data sheet
% add in uncertainty in duct temperature measurement given that TC is a
% regular k-type TC
Propane_burner_uncert = uncert_Te_K_k_type(Propane_burner_uncert, ...
    Propane_burner.Te_K);

% Uncertainties in duct specfication
Propane_burner_uncert.D_duct_m = 0.0005;

% Calibration bottles with uncertainties
Propane_burner.x_CO_bottle_span = 4/100;
Propane_burner.x_CO2_bottle_span = 8/100;
Propane_burner.x_O2_bottle_span = 21/100;
% Propane_burner.x_O2_bottle_span = 100/100;
Propane_burner.x_i_bottle_zero = 0; %N2 is used for zero, in this 
% case 99.99% N2 was used. 

Propane_burner_uncert.x_CO_bottle_span = 0.02*Propane_burner.x_CO_bottle_span...
    /(3^0.5); % convert to stanard error
Propane_burner_uncert.x_CO2_bottle_span = 0.02*Propane_burner.x_CO2_bottle_span...
    /(3^0.5); % convert to stanard error
Propane_burner_uncert.x_O2_bottle_span = 0.02*Propane_burner.x_O2_bottle_span...
    /(3^0.5); % convert to stanard error
%Propane_burner_uncert.x_O2_bottle_span = 0.1/100/(3^0.5); 
              % convert to stanard error
% Calulate the uncertainty for the zero given nitrogen concentration
x_N2_bottle = 99.9/100;
Propane_burner_uncert.x_i_bottle_zero = (1-x_N2_bottle)/(3^0.5); 
            % converts to stanard error

% Run functions to calculate the uncertainty in concentration and flow 
% measurements. 
% Please note that these will not be used for HRR are just for reference.
% Instead calculation will be reapeted when calcuating the HRR.
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

% alter check only O2 old code
% % x0 = dlarray([Propane_burner.x_O2_amb, Propane_burner.x_O2_measure, Propane_burner.Delta_P_Pa, Propane_burner.Te_K, Propane_burner.D_duct_m, Propane_burner.f_Re, Propane_burner.alpha, ...
% %     Propane_burner.E_kJ_per_kg]);
% % [HRR_kW_dl, partials_dl] = dlfeval(@HRR_O2,x0);



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

% Calculate the diamter for a circular pool fire to obtain the given HRR
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

% Uncertainities in flow measurement device, assuming no change in ambient
% pressure
Heptane_pool_uncert.Delta_P_Pa = 0.0015*Heptane_pool.Delta_P_Pa; 
       % from pressure transducter data sheet
% add in uncertainty in duct temperature measurement given that TC is a
% regular k-type TC
Heptane_pool_uncert = uncert_Te_K_k_type(Heptane_pool_uncert, ...
    Heptane_pool.Te_K);

% Uncertainties in duct specfication
Heptane_pool_uncert.D_duct_m = 0.0005;

% Calibration bottles with uncertainties
Heptane_pool.x_CO_bottle_span = 5/100;
Heptane_pool.x_CO2_bottle_span = 10/100;
Heptane_pool.x_O2_bottle_span = 21/100;
Heptane_pool.x_i_bottle_zero = 0; %N2 is used for zero, in this 
% case 99.99% N2 was used. 

Heptane_pool_uncert.x_CO_bottle_span = 0.02*Heptane_pool.x_CO_bottle_span...
    /(3^0.5); % convert to stanard error
Heptane_pool_uncert.x_CO2_bottle_span = 0.02*Heptane_pool.x_CO2_bottle_span...
    /(3^0.5); % convert to stanard error
Heptane_pool_uncert.x_O2_bottle_span = 0.02*Heptane_pool.x_O2_bottle_span...
    /(3^0.5); % convert to stanard error
% Calulate the uncertainty for the zero given nitrogen concentration
x_N2_bottle = 99.9/100;
Heptane_pool_uncert.x_i_bottle_zero = (1-x_N2_bottle)/(3^0.5); 
            % converts to stanard error

% Run functions to calculate the uncertainty in concentration and flow 
% measurements. 
% Please note that these will not be used for HRR are just for reference.
% Instead calculation will be reapeted when calcuating the HRR.
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



%% Create plot for HRR Models for fire diamter and flame height
figure
tiledlayout(1,2, 'TileSpacing','compact', 'Padding','compact')

nexttile
hold on
box on

xlabel('HRR [kW]', 'FontName', 'Times', 'FontSize', 8,'FontWeight', 'bold')
ylabel('D_f [m]', 'FontName', 'Times', 'FontSize', 8, ...
    'FontWeight', 'bold')

%limit the axis
xlim([0 750])
ylim([0 1])

% Plot pool diamter
% plot(HRR_kw, Propane_burner.D_fire_m, 'linewidth', 0.75, 'Color', ...
%    '#5ec962', 'LineStyle','-')
plot(HRR_kw, Heptane_pool.D_fire_m, 'linewidth', 0.75, 'Color', ...
    '#3b528b', 'LineStyle','--')

% label
title('a)', 'FontSize', 8, 'FontWeight','bold')

%specifs axis numbers size and font
ax = gca;
ax.FontSize = 8;
ax.FontName = 'Times';
ax.TitleHorizontalAlignment = 'left';
ax.YAxis.TickValues = [0 0.25 .5 0.75 1];

nexttile
hold on
box on

xlabel('HRR [kW]', 'FontName', 'Times', 'FontSize', 8,'FontWeight', 'bold')
ylabel('l_f [m]', 'FontName', 'Times', 'FontSize', 8, ...
    'FontWeight', 'bold')

%limit the x-axis
xlim([0 750])
ylim([0 3])

% Plot hood height
plot_space = linspace(0, 750, spacing);
h_hood_m = ones(size(plot_space));
h_hood_m = Propane_burner.h_hood_m*h_hood_m;
h_hood_m = reshape(h_hood_m, [5, spacing/5]);
h_hood_m(5, :) = NaN;
h_hood_m = reshape(h_hood_m, [spacing, 1]);
plot(plot_space, h_hood_m, ...
    'linewidth', 0.5, 'HandleVisibility','off',LineStyle='-', Color='k')
% using hezkestad virtual orgin correction 
plot(HRR_kw, Propane_burner.flame_height_m, 'linewidth', 0.75, 'Color', ...
    '#5ec962', 'LineStyle','-')
plot(HRR_kw, Heptane_pool.flame_height_m, 'linewidth', 0.75, 'Color', ...
    '#3b528b', 'LineStyle','--')

% label
title('b)', 'FontSize', 8, 'FontWeight','bold')
text(20, 2.32, 'z_{hood}', 'FontSize', 8)

%specifs axis numbers size and font
ax = gca;
ax.FontSize = 8;
ax.FontName = 'Times';
ax.TitleHorizontalAlignment = 'left';

lg = legend({ 'Propane', 'Heptane'}, FontSize=8, Box="off");
% lg.Layout.Tile = 'South';
lg.Location = 'southeast';
% lg.NumColumns = 2;

%Save the figure as a pdf
set(gcf,'PaperUnits','inches','PaperPosition', [0 0 4 2], 'PaperSize', [4 2]);
print -dpdf Figures/D_l_f.pdf -r150

%% Create plot for HRR Models results (Te, deltaP, volume fraction 
%   measurements)
figure
tiledlayout(2,6, 'TileSpacing','compact', 'Padding','compact')

nexttile([1 2])
hold on
box on

xlabel('HRR [kW]', 'FontName', 'Times', 'FontSize', 8,'FontWeight', 'bold')
ylabel('X_{O_2,e} [-]', 'FontName', 'Times', 'FontSize', 8, ...
    'FontWeight', 'bold')

%limit the axis
xlim([0 750])
ylim([0.19 0.21])

% x_O2
plot(HRR_kw, Propane_burner.x_O2_measure, 'linewidth', 0.75, 'Color', ...
   '#5ec962', 'LineStyle','-')
plot(HRR_kw, Heptane_pool.x_O2_measure, 'linewidth', 0.75, 'Color', ...
    '#3b528b', 'LineStyle','--')

% label
title('a)', 'FontSize', 8, 'FontWeight','bold')

%specifs axis numbers size and font
ax = gca;
ax.FontSize = 8;
ax.FontName = 'Times';
ax.TitleHorizontalAlignment = 'left';


nexttile([1 2])
hold on
box on

xlabel('HRR [kW]', 'FontName', 'Times', 'FontSize', 8,'FontWeight', 'bold')
ylabel('X_{CO_2,e} [-]', 'FontName', 'Times', 'FontSize', 8, ...
    'FontWeight', 'bold')

%limit the x-axis
xlim([0 750])
% ylim([0 3])


% Plot X_CO2
plot(HRR_kw, Propane_burner.x_CO2_measure, 'linewidth', 0.75, 'Color', ...
    '#5ec962', 'LineStyle','-')
plot(HRR_kw, Heptane_pool.x_CO2_measure, 'linewidth', 0.75, 'Color', ...
    '#3b528b', 'LineStyle','--')

% label
title('b)', 'FontSize', 8, 'FontWeight','bold')

%specifs axis numbers size and font
ax = gca;
ax.FontSize = 8;
ax.FontName = 'Times';
ax.TitleHorizontalAlignment = 'left';

% x_CO
nexttile([1 2])
hold on
box on

xlabel('HRR [kW]', 'FontName', 'Times', 'FontSize', 8,'FontWeight', 'bold')
ylabel('X_{CO,e} [-]', 'FontName', 'Times', 'FontSize', 8, ...
    'FontWeight', 'bold')

%limit the x-axis
xlim([0 750])
% ylim([0 3])

plot(HRR_kw, Propane_burner.x_CO_measure, 'linewidth', 0.75, 'Color', ...
    '#5ec962', 'LineStyle','-')
plot(HRR_kw, Heptane_pool.x_CO_measure, 'linewidth', 0.75, 'Color', ...
    '#3b528b', 'LineStyle','--')

% label
title('c)', 'FontSize', 8, 'FontWeight','bold')

%specifs axis numbers size and font
ax = gca;
ax.FontSize = 8;
ax.FontName = 'Times';
ax.TitleHorizontalAlignment = 'left';

% Temperature
nexttile(8, [1 2])
hold on
box on

xlabel('HRR [kW]', 'FontName', 'Times', 'FontSize', 8,'FontWeight', 'bold')
ylabel('T_e [K]', 'FontName', 'Times', 'FontSize', 8, ...
    'FontWeight', 'bold')

%limit the x-axis
xlim([0 750])
% ylim([0 3])

plot(HRR_kw, Propane_burner.Te_K, 'linewidth', 0.75, 'Color', ...
    '#5ec962', 'LineStyle','-')
plot(HRR_kw, Heptane_pool.Te_K, 'linewidth', 0.75, 'Color', ...
    '#3b528b', 'LineStyle','--')

% label
title('d)', 'FontSize', 8, 'FontWeight','bold')

%specifs axis numbers size and font
ax = gca;
ax.FontSize = 8;
ax.FontName = 'Times';
ax.TitleHorizontalAlignment = 'left';

lg = legend({ 'Propane', 'Heptane'}, FontSize=8, Box="off");
% lg.Layout.Tile = 'East';
lg.Location = 'southeast';
% lg.NumColumns = 2;

%%% Delta P
nexttile([1 2])
hold on
box on

xlabel('HRR [kW]', 'FontName', 'Times', 'FontSize', 8,'FontWeight', 'bold')
ylabel('\Delta P [Pa]', 'FontName', 'Times', 'FontSize', 8, ...
    'FontWeight', 'bold')

%limit the x-axis
xlim([0 750])
% ylim([0 3])

plot(HRR_kw, Propane_burner.Delta_P_Pa, 'linewidth', 0.75, 'Color', ...
    '#5ec962', 'LineStyle','-')
plot(HRR_kw, Heptane_pool.Delta_P_Pa, 'linewidth', 0.75, 'Color', ...
    '#3b528b', 'LineStyle','--')

% label
title('e)', 'FontSize', 8, 'FontWeight','bold')

%specifs axis numbers size and font
ax = gca;
ax.FontSize = 8;
ax.FontName = 'Times';
ax.TitleHorizontalAlignment = 'left';


%Save the figure as a pdf
set(gcf,'PaperUnits','inches','PaperPosition', [0 0 6 4], 'PaperSize', [6 4]);
print -dpdf Figures/model_results.pdf -r150

%% Plot Calculated HRR comparsion and Uncertainty

figure
tiledlayout(1,2, 'TileSpacing','compact', 'Padding','compact')

nexttile
hold on
box on

xlabel('HRR [kW]', 'FontName', 'Times', 'FontSize', 8,'FontWeight', 'bold')
% ylabel('$(\dot{Q}_{w/o}-\dot{Q}_{w/})/\dot{Q}_{w/}$ \bf{(\%)}', ...
%     'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 8, ...
%     'FontWeight', 'bold')
ylabel('$(\dot{Q}_{1}-\dot{Q}_{2})/\dot{Q}_{2}$ \bf{[\%]}', ...
    'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 8, ...
    'FontWeight', 'bold')

%limit the axis
xlim([0 750])
ylim([0 20])

% Plot HRR differences
plot(HRR_kw, (Propane_burner_HRR_calc_kW_without - ...
    Propane_burner_HRR_calc_kW_with)./Propane_burner_HRR_calc_kW_with*100, ...
    'linewidth', 0.75, 'Color', '#5ec962', 'LineStyle','-')
plot(HRR_kw, (Heptane_pool_HRR_calc_kW_without - ...
    Heptane_pool_HRR_calc_kW_with)./Heptane_pool_HRR_calc_kW_with*100, ...
    'linewidth', 0.75, 'Color', '#3b528b', 'LineStyle','--')

% label
title('a)', 'FontSize', 8, 'FontWeight','bold')

%specifs axis numbers size and font
ax = gca;
ax.FontSize = 8;
ax.FontName = 'Times';
ax.TitleHorizontalAlignment = 'left';
%ax.YAxis.TickValues = [0 0.25 .5 0.75 1];

lg = legend({ 'Propane', 'Heptane'}, FontSize=8, Box="off");
% lg.Layout.Tile = 'South';
lg.Location = 'southeast';
% lg.NumColumns = 2;

nexttile
hold on
box on

xlabel('HRR [kW]', 'FontName', 'Times', 'FontSize', 8,'FontWeight', 'bold')
ylabel('$U(\dot{Q})/\dot{Q}$ \bf{[\%]}', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 8, ...
    'FontWeight', 'bold')

%limit the x-axis
xlim([0 750])
ylim([0 20])


% O2 only
plot(HRR_kw(1:10:end), Propane_burner_uncert_without.HRR_kW(1:10:end)./...
    Propane_burner_HRR_calc_kW_without(1:10:end)*100, 'linewidth', 0.75, 'Color', ...
    '#5ec962', 'LineStyle','-', 'Marker', 'o', 'MarkerSize', 3)
% The fraction of the heat realease rate 
plot(HRR_kw, Propane_burner_uncert_with.HRR_kW./HRR_kw*100, 'linewidth', 0.75, 'Color', ...
    '#5ec962', 'LineStyle','-')


plot(HRR_kw(1:10:end), Heptane_pool_uncert_without.HRR_kW(1:10:end)./...
    Heptane_pool_HRR_calc_kW_without(1:10:end)*100, 'linewidth', 0.75, 'Color', ...
    '#3b528b', 'LineStyle','--', 'Marker', '^', 'MarkerSize', 3)
plot(HRR_kw, Heptane_pool_uncert_with.HRR_kW./HRR_kw*100, 'linewidth', 0.75, 'Color', ...
    '#3b528b', 'LineStyle','--')

% label
title('b)', 'FontSize', 8, 'FontWeight','bold')
%text(20, 2.32, 'z_{hood}', 'FontSize', 8)

%specifs axis numbers size and font
ax = gca;
ax.FontSize = 8;
ax.FontName = 'Times';
ax.TitleHorizontalAlignment = 'left';

lg = legend({ 'Propane ($\dot{Q}_1$)', 'Propane ($\dot{Q}_2$)', ...
    'Heptane ($\dot{Q}_1$)', 'Heptane ($\dot{Q}_2$)'}, ...
    'Interpreter', 'latex', FontSize=8, Box="off");
% lg.Layout.Tile = 'South';
lg.Location = 'southeast';
% lg.NumColumns = 2;

%Save the figure as a pdf
set(gcf,'PaperUnits','inches','PaperPosition', [0 0 4 2], 'PaperSize', [4 2]);
print -dpdf Figures/HRR_uncert_comp.pdf -r150

%% Plot HRR Uncertainty contibution fractions

figure
tiledlayout(1,2, 'TileSpacing','none', 'Padding','compact')

nexttile
hold on
box on

xlabel('HRR [kW]', 'FontName', 'Times', 'FontSize', 8,'FontWeight', 'bold')
ylabel('$(u_c(\dot{x}))^2/(u_c(\dot{Q}))^2$ \bf{[\%]}', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 8, ...
    'FontWeight', 'bold')

%limit the axis
xlim([0 750])
ylim([1 100])
yscale("log")

% Plot uncert fracs
% O2 span gas
plot(HRR_kw(1:length(Propane_burner_uncert_frac.cal_gas(:,3))),...
    Propane_burner_uncert_frac.cal_gas(:,3)*100, ...
    'linewidth', 0.75, 'Color', '#fde725', 'LineStyle','-')
scatter(HRR_kw(1:10:length(Propane_burner_uncert_frac.cal_gas(:,3))),...
    Propane_burner_uncert_frac.cal_gas(1:10:end,3)*100, 6, ...
    'MarkerEdgeColor', '#fde725', 'MarkerFaceColor', '#fde725', 'Marker', ...
    'o')
% Nitrogen calibration gas
plot(HRR_kw(1:length(Propane_burner_uncert_frac.cal_gas(:,3))),...
    Propane_burner_uncert_frac.cal_gas(:,4)*100, ...
    'linewidth', 0.75, 'Color', '#c8e020', 'LineStyle','-')
scatter(HRR_kw(1:10:length(Propane_burner_uncert_frac.cal_gas(:,3))),...
    Propane_burner_uncert_frac.cal_gas(1:10:end,4)*100, 6, ...
    'MarkerEdgeColor', '#c8e020', 'MarkerFaceColor', '#c8e020', 'Marker', ...
    'square')
% e_int CO2
plot(HRR_kw(1:length(Propane_burner_uncert_frac.gas_cons(:,3))),...
    Propane_burner_uncert_frac.misc_CO2(:,1)*100, ...
    'linewidth', 0.75, 'Color', '#90d743', 'LineStyle','-')
scatter(HRR_kw(1:10:length(Propane_burner_uncert_frac.cal_gas(:,3))),...
    Propane_burner_uncert_frac.misc_CO2(1:10:end,1)*100, 6, ...
    'MarkerEdgeColor', '#90d743', 'MarkerFaceColor', '#90d743', 'Marker', ...
    'v')
% e_lin CO2
plot(HRR_kw(1:length(Propane_burner_uncert_frac.gas_cons(:,3))),...
    Propane_burner_uncert_frac.misc_CO2(:,2)*100, ...
    'linewidth', 0.75, 'Color', '#5ec962', 'LineStyle','-')
scatter(HRR_kw(1:10:length(Propane_burner_uncert_frac.cal_gas(:,3))),...
    Propane_burner_uncert_frac.misc_CO2(1:10:end,2)*100, 6, ...
    'MarkerEdgeColor', '#5ec962', 'MarkerFaceColor', '#5ec962', 'Marker', ...
    '<')
% e_rep CO2
plot(HRR_kw(1:length(Propane_burner_uncert_frac.gas_cons(:,3))),...
    Propane_burner_uncert_frac.misc_CO2(:,3)*100, ...
    'linewidth', 0.75, 'Color', '#35b779', 'LineStyle','-')
scatter(HRR_kw(1:10:length(Propane_burner_uncert_frac.cal_gas(:,3))),...
    Propane_burner_uncert_frac.misc_CO2(1:10:end,3)*100, 6, ...
    'MarkerEdgeColor', '#35b779', 'MarkerFaceColor', '#35b779', 'Marker', ...
    '>')
% e_int CO
plot(HRR_kw(1:length(Propane_burner_uncert_frac.gas_cons(:,3))),...
    Propane_burner_uncert_frac.misc_CO(:,1)*100, ...
    'linewidth', 0.75, 'Color', '#20a486', 'LineStyle','-')
scatter(HRR_kw(1:10:length(Propane_burner_uncert_frac.cal_gas(:,3))),...
    Propane_burner_uncert_frac.misc_CO(1:10:end,1)*100, 6, ...
    'MarkerEdgeColor', '#20a486', 'MarkerFaceColor', 'none', 'Marker', ...
    'v')
% e_lin CO
plot(HRR_kw(1:length(Propane_burner_uncert_frac.gas_cons(:,3))),...
    Propane_burner_uncert_frac.misc_CO(:,2)*100, ...
    'linewidth', 0.75, 'Color', '#21918c', 'LineStyle','-')
scatter(HRR_kw(1:10:length(Propane_burner_uncert_frac.cal_gas(:,3))),...
    Propane_burner_uncert_frac.misc_CO(1:10:end,2)*100, 6, ...
    'MarkerEdgeColor', '#21918c', 'MarkerFaceColor', 'none', 'Marker', ...
    '<')
% e_rep CO
plot(HRR_kw(1:length(Propane_burner_uncert_frac.gas_cons(:,3))),...
    Propane_burner_uncert_frac.misc_CO(:,3)*100, ...
    'linewidth', 0.75, 'Color', '#287c8e', 'LineStyle','-')
scatter(HRR_kw(1:10:length(Propane_burner_uncert_frac.cal_gas(:,3))),...
    Propane_burner_uncert_frac.misc_CO(1:10:end,3)*100, 6, ...
    'MarkerEdgeColor', '#287c8e', 'MarkerFaceColor', 'none', 'Marker', ...
    '>')
%f(Re)
plot(HRR_kw(1:length(Propane_burner_uncert_frac.gas_cons(:,3))),...
    Propane_burner_uncert_frac.flow_vars(:,1)*100, ...
    'linewidth', 0.75, 'Color', '#31688e', 'LineStyle','-')
scatter(HRR_kw(1:10:length(Propane_burner_uncert_frac.cal_gas(:,3))),...
    Propane_burner_uncert_frac.flow_vars(1:10:end,1)*100, 6, ...
    'MarkerEdgeColor', '#31688e', 'MarkerFaceColor', '#31688e', 'Marker', ...
    'pentagram')
%E
plot(HRR_kw(1:length(Propane_burner_uncert_frac.gas_cons(:,3))),...
    Propane_burner_uncert_frac.OCC_cons(:,2)*100, ...
    'linewidth', 0.75, 'Color', '#3b528b', 'LineStyle','-')
scatter(HRR_kw(1:10:length(Propane_burner_uncert_frac.cal_gas(:,3))),...
    Propane_burner_uncert_frac.OCC_cons(1:10:end,2)*100, 6, ...
    'MarkerEdgeColor', '#3b528b', 'MarkerFaceColor', 'none', 'Marker', ...
    'x')
%alpha
plot(HRR_kw(1:length(Propane_burner_uncert_frac.gas_cons(:,3))),...
    Propane_burner_uncert_frac.OCC_cons(:,1)*100, ...
    'linewidth', 0.75, 'Color', '#443983', 'LineStyle','-')
scatter(HRR_kw(1:10:length(Propane_burner_uncert_frac.cal_gas(:,3))),...
    Propane_burner_uncert_frac.OCC_cons(1:10:end,1)*100, 6, ...
    'MarkerEdgeColor', '#443983', 'MarkerFaceColor', '#443983', 'Marker', ...
    '^')
%P_amb
plot(HRR_kw(1:length(Propane_burner_uncert_frac.gas_cons(:,3))),...
    Propane_burner_uncert_frac.gas_cons(:,2)*100, ...
    'linewidth', 0.75, 'Color', '#481f70', 'LineStyle','-')
scatter(HRR_kw(1:10:length(Propane_burner_uncert_frac.cal_gas(:,3))),...
    Propane_burner_uncert_frac.gas_cons(1:10:end,2)*100, 6, ...
    'MarkerEdgeColor', '#481f70', 'MarkerFaceColor', '#481f70', 'Marker', ...
    'd')
%Mw_air
plot(HRR_kw(1:length(Propane_burner_uncert_frac.gas_cons(:,3))),...
    Propane_burner_uncert_frac.gas_cons(:,3)*100, ...
    'linewidth', 0.75, 'Color', '#440154', 'LineStyle','-')
scatter(HRR_kw(1:10:length(Propane_burner_uncert_frac.cal_gas(:,3))),...
    Propane_burner_uncert_frac.gas_cons(1:10:end,3)*100, 6, ...
    'MarkerEdgeColor', '#440154', 'MarkerFaceColor', 'none', 'Marker', ...
    'o')

% label
title('Propane', 'FontSize', 8, 'FontWeight','bold')

%specifs axis numbers size and font
ax = gca;
ax.FontSize = 8;
ax.FontName = 'Times';
ax.TitleHorizontalAlignment = 'center';
ax.YAxis.TickValues = [10e-5 10e-4 10e-3 10e-2 10e-1 10e-0, 10e1, 10e2];


nexttile
hold on
box on

xlabel('HRR (kW)', 'FontName', 'Times', 'FontSize', 8,'FontWeight', 'bold')
% ylabel('$(u_c(\dot{x}))^2/(u_c(\dot{Q}))^2$ \bf{(\%)}', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 8, ...
%     'FontWeight', 'bold')


%limit the x-axis
xlim([0 750])
ylim([0 100])
%yscale("log")

% Heptane pool
% O2 span gas
plot(HRR_kw(1:length( Heptane_pool_uncert_frac.cal_gas(:,3))),...
     Heptane_pool_uncert_frac.cal_gas(:,3)*100, ...
    'linewidth', 0.75, 'Color', '#fde725', 'LineStyle','-')
scatter(HRR_kw(1:10:length( Heptane_pool_uncert_frac.cal_gas(:,3))),...
    Heptane_pool_uncert_frac.cal_gas(1:10:end,3)*100, 6, ...
    'MarkerEdgeColor', '#fde725', 'MarkerFaceColor', '#fde725', 'Marker', ...
    'o')
% Nitrogen calibration gas
plot(HRR_kw(1:length( Heptane_pool_uncert_frac.cal_gas(:,3))),...
     Heptane_pool_uncert_frac.cal_gas(:,4)*100, ...
    'linewidth', 0.75, 'Color', '#c8e020', 'LineStyle','-')
scatter(HRR_kw(1:10:length( Heptane_pool_uncert_frac.cal_gas(:,3))),...
     Heptane_pool_uncert_frac.cal_gas(1:10:end,4)*100, 6, ...
    'MarkerEdgeColor', '#c8e020', 'MarkerFaceColor', '#c8e020', 'Marker', ...
    'square')
% e_int CO2
plot(HRR_kw(1:length( Heptane_pool_uncert_frac.gas_cons(:,3))),...
     Heptane_pool_uncert_frac.misc_CO2(:,1)*100, ...
    'linewidth', 0.75, 'Color', '#c8e020', 'LineStyle','-')
scatter(HRR_kw(1:10:length( Heptane_pool_uncert_frac.cal_gas(:,3))),...
     Heptane_pool_uncert_frac.misc_CO2(1:10:end,1)*100, 6, ...
    'MarkerEdgeColor', '#c8e020', 'MarkerFaceColor', '#c8e020', 'Marker', ...
    'v')
% e_lin CO2
plot(HRR_kw(1:length( Heptane_pool_uncert_frac.gas_cons(:,3))),...
     Heptane_pool_uncert_frac.misc_CO2(:,2)*100, ...
    'linewidth', 0.75, 'Color', '#5ec962', 'LineStyle','-')
scatter(HRR_kw(1:10:length( Heptane_pool_uncert_frac.cal_gas(:,3))),...
     Heptane_pool_uncert_frac.misc_CO2(1:10:end,2)*100, 6, ...
    'MarkerEdgeColor', '#5ec962', 'MarkerFaceColor', '#5ec962', 'Marker', ...
    '<')
% e_rep CO2
plot(HRR_kw(1:length( Heptane_pool_uncert_frac.gas_cons(:,3))),...
     Heptane_pool_uncert_frac.misc_CO2(:,3)*100, ...
    'linewidth', 0.75, 'Color', '#35b779', 'LineStyle','-')
scatter(HRR_kw(1:10:length( Heptane_pool_uncert_frac.cal_gas(:,3))),...
     Heptane_pool_uncert_frac.misc_CO2(1:10:end,3)*100, 6, ...
    'MarkerEdgeColor', '#35b779', 'MarkerFaceColor', '#35b779', 'Marker', ...
    '>')
% e_int CO
plot(HRR_kw(1:length( Heptane_pool_uncert_frac.gas_cons(:,3))),...
     Heptane_pool_uncert_frac.misc_CO(:,1)*100, ...
    'linewidth', 0.75, 'Color', '#20a486', 'LineStyle','-')
scatter(HRR_kw(1:10:length( Heptane_pool_uncert_frac.cal_gas(:,3))),...
     Heptane_pool_uncert_frac.misc_CO(1:10:end,1)*100, 6, ...
    'MarkerEdgeColor', '#20a486', 'MarkerFaceColor', 'none', 'Marker', ...
    'v')
% e_lin CO
plot(HRR_kw(1:length( Heptane_pool_uncert_frac.gas_cons(:,3))),...
     Heptane_pool_uncert_frac.misc_CO(:,2)*100, ...
    'linewidth', 0.75, 'Color', '#21918c', 'LineStyle','-')
scatter(HRR_kw(1:10:length( Heptane_pool_uncert_frac.cal_gas(:,3))),...
     Heptane_pool_uncert_frac.misc_CO(1:10:end,2)*100, 6, ...
    'MarkerEdgeColor', '#21918c', 'MarkerFaceColor', 'none', 'Marker', ...
    '<')
% e_rep CO
plot(HRR_kw(1:length( Heptane_pool_uncert_frac.gas_cons(:,3))),...
     Heptane_pool_uncert_frac.misc_CO(:,3)*100, ...
    'linewidth', 0.75, 'Color', '#287c8e', 'LineStyle','-')
scatter(HRR_kw(1:10:length( Heptane_pool_uncert_frac.cal_gas(:,3))),...
     Heptane_pool_uncert_frac.misc_CO(1:10:end,3)*100, 6, ...
    'MarkerEdgeColor', '#287c8e', 'MarkerFaceColor', 'none', 'Marker', ...
    '>')
%f(Re)
plot(HRR_kw(1:length( Heptane_pool_uncert_frac.gas_cons(:,3))),...
     Heptane_pool_uncert_frac.flow_vars(:,1)*100, ...
    'linewidth', 0.75, 'Color', '#287c8e', 'LineStyle','-')
scatter(HRR_kw(1:10:length( Heptane_pool_uncert_frac.cal_gas(:,3))),...
     Heptane_pool_uncert_frac.flow_vars(1:10:end,1)*100, 6, ...
    'MarkerEdgeColor', '#287c8e', 'MarkerFaceColor', '#287c8e', 'Marker', ...
    'pentagram')
%E
plot(HRR_kw(1:length( Heptane_pool_uncert_frac.gas_cons(:,3))),...
     Heptane_pool_uncert_frac.OCC_cons(:,2)*100, ...
    'linewidth', 0.75, 'Color', '#3b528b', 'LineStyle','-')
scatter(HRR_kw(1:10:length( Heptane_pool_uncert_frac.cal_gas(:,3))),...
     Heptane_pool_uncert_frac.OCC_cons(1:10:end,2)*100, 6, ...
    'MarkerEdgeColor', '#3b528b', 'MarkerFaceColor', 'none', 'Marker', ...
    'x')
%alpha
plot(HRR_kw(1:length( Heptane_pool_uncert_frac.gas_cons(:,3))),...
     Heptane_pool_uncert_frac.OCC_cons(:,1)*100, ...
    'linewidth', 0.75, 'Color', '#3b528b', 'LineStyle','-')
scatter(HRR_kw(1:10:length( Heptane_pool_uncert_frac.cal_gas(:,3))),...
     Heptane_pool_uncert_frac.OCC_cons(1:10:end,1)*100, 6, ...
    'MarkerEdgeColor', '#3b528b', 'MarkerFaceColor', '#3b528b', 'Marker', ...
    '^')
%P_amb
plot(HRR_kw(1:length(Heptane_pool_uncert_frac.gas_cons(:,3))),...
    Heptane_pool_uncert_frac.gas_cons(:,2)*100, ...
    'linewidth', 0.75, 'Color', '#481f70', 'LineStyle','-')
scatter(HRR_kw(1:10:length(Heptane_pool_uncert_frac.cal_gas(:,3))),...
    Heptane_pool_uncert_frac.gas_cons(1:10:end,2)*100, 6, ...
    'MarkerEdgeColor', '#481f70', 'MarkerFaceColor', '#481f70', 'Marker', ...
    'd')
%Mw_air
plot(HRR_kw(1:length(Heptane_pool_uncert_frac.gas_cons(:,3))),...
    Heptane_pool_uncert_frac.gas_cons(:,3)*100, ...
    'linewidth', 0.75, 'Color', '#440154', 'LineStyle','-')
scatter(HRR_kw(1:10:length(Heptane_pool_uncert_frac.cal_gas(:,3))),...
    Heptane_pool_uncert_frac.gas_cons(1:10:end,3)*100, 6, ...
    'MarkerEdgeColor', '#440154', 'MarkerFaceColor', 'none', 'Marker', ...
    'o')

% label
title('Heptane', 'FontSize', 8, 'FontWeight','bold')
%text(20, 2.32, 'z_{hood}', 'FontSize', 8)

%specifs axis numbers size and font
ax = gca;
ax.FontSize = 8;
ax.FontName = 'Times';
ax.TitleHorizontalAlignment = 'left';
% ax.YAxis.TickValues = [10e-5 10e-4 10e-3 10e-2 10e-1 10e-0, 10e1, 10e2];
ax.YAxis.TickLabel = [];

lg = legend({ '', 'X_{O_2}^{span}', '', 'X_{i}^{zero}', '', ...
    '\epsilon _{CO_2,int}', '', '\epsilon _{CO_2,lin}', '', ...
    '\epsilon _{CO_2,rep}','', '\epsilon _{CO,int}', '', '\epsilon _{CO,lin}', ...
    '', '\epsilon _{CO,rep}', '', 'f(Re)', '', 'E', '', '\alpha', '', ...
    'P_{amb}', '', 'Mw_{air}'}, FontSize=8, Box="off");
lg.Layout.Tile = 'East';
%lg.Location = 'northeast';
lg.NumColumns = 2;

%Save the figure as a pdf
set(gcf,'PaperUnits','inches','PaperPosition', [0 0 5.5 2], 'PaperSize', [5.5 2]);
print -dpdf Figures/Uncertainty_fraction.pdf -r150

%% Plot HRR Uncertainty contibution fractions only vars greater than 1% contribution

figure
tiledlayout(1,2, 'TileSpacing','none', 'Padding','compact')

nexttile
hold on
box on

xlabel('HRR [kW]', 'FontName', 'Times', 'FontSize', 8,'FontWeight', 'bold')
ylabel({'\bf{Uncertainty Budget}' ; '$u_c^2(x)/u_c^2(\dot{Q})$ \bf{[\%]}'}, 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 8, ...
    'FontWeight', 'bold')


%limit the axis
xlim([0 750])
ylim([0 100])
%yscale("log")

% Plot uncert fracs
% O2 span gas
plot(HRR_kw(1:length(Propane_burner_uncert_frac.cal_gas(:,3))),...
    Propane_burner_uncert_frac.cal_gas(:,3)*100, ...
    'linewidth', 0.75, 'Color', '#fde725', 'LineStyle','-')
scatter(HRR_kw(1:10:length(Propane_burner_uncert_frac.cal_gas(:,3))),...
    Propane_burner_uncert_frac.cal_gas(1:10:end,3)*100, 6, ...
    'MarkerEdgeColor', '#fde725', 'MarkerFaceColor', '#fde725', 'Marker', ...
    'o')
%f(Re)
plot(HRR_kw(1:length(Propane_burner_uncert_frac.gas_cons(:,3))),...
    Propane_burner_uncert_frac.flow_vars(:,1)*100, ...
    'linewidth', 0.75, 'Color', '#35b779', 'LineStyle','-')
scatter(HRR_kw(1:10:length(Propane_burner_uncert_frac.cal_gas(:,3))),...
    Propane_burner_uncert_frac.flow_vars(1:10:end,1)*100, 6, ...
    'MarkerEdgeColor', '#35b779', 'MarkerFaceColor', '#35b779', 'Marker', ...
    'pentagram')
%E
plot(HRR_kw(1:length(Propane_burner_uncert_frac.gas_cons(:,3))),...
    Propane_burner_uncert_frac.OCC_cons(:,2)*100, ...
    'linewidth', 0.75, 'Color', '#31688e', 'LineStyle','-')
scatter(HRR_kw(1:10:length(Propane_burner_uncert_frac.cal_gas(:,3))),...
    Propane_burner_uncert_frac.OCC_cons(1:10:end,2)*100, 6, ...
    'MarkerEdgeColor', '#31688e', 'MarkerFaceColor', '#31688e', 'Marker', ...
    's')
% e CO
plot(HRR_kw(1:length(Propane_burner_uncert_frac.gas_cons(:,3))),...
    Propane_burner_uncert_frac.misc_CO(:,1)*100*3, ...
    'linewidth', 0.75, 'Color', '#440154', 'LineStyle','-')
scatter(HRR_kw(1:10:length(Propane_burner_uncert_frac.cal_gas(:,3))),...
    Propane_burner_uncert_frac.misc_CO(1:10:end,1)*100*3, 6, ...
    'MarkerEdgeColor', '#440154', 'MarkerFaceColor', '#440154', 'Marker', ...
    '^')
% label
title('Propane', 'FontSize', 8, 'FontWeight','bold')

%specifs axis numbers size and font
ax = gca;
ax.FontSize = 8;
ax.FontName = 'Times';
ax.TitleHorizontalAlignment = 'center';
%ax.YAxis.TickValues = [10e-5 10e-4 10e-3 10e-2 10e-1 10e-0, 10e1, 10e2];


nexttile
hold on
box on

xlabel('HRR [kW]', 'FontName', 'Times', 'FontSize', 8,'FontWeight', 'bold')
% ylabel('$(u_c(\dot{x}))^2/(u_c(\dot{Q}))^2$ \bf{(\%)}', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 8, ...
%     'FontWeight', 'bold')


%limit the x-axis
xlim([0 750])
ylim([0 100])
%yscale("log")

% Heptane pool
% O2 span gas
% Plot uncert fracs
% O2 span gas
plot(HRR_kw(1:length(Heptane_pool_uncert_frac.cal_gas(:,3))),...
    Heptane_pool_uncert_frac.cal_gas(:,3)*100, ...
    'linewidth', 0.75, 'Color', '#fde725', 'LineStyle','-')
scatter(HRR_kw(1:10:length(Heptane_pool_uncert_frac.cal_gas(:,3))),...
    Heptane_pool_uncert_frac.cal_gas(1:10:end,3)*100, 6, ...
    'MarkerEdgeColor', '#fde725', 'MarkerFaceColor', '#fde725', 'Marker', ...
    'o')
%f(Re)
plot(HRR_kw(1:length(Heptane_pool_uncert_frac.gas_cons(:,3))),...
    Heptane_pool_uncert_frac.flow_vars(:,1)*100, ...
    'linewidth', 0.75, 'Color', '#35b779', 'LineStyle','-')
scatter(HRR_kw(1:10:length(Heptane_pool_uncert_frac.cal_gas(:,3))),...
    Heptane_pool_uncert_frac.flow_vars(1:10:end,1)*100, 6, ...
    'MarkerEdgeColor', '#35b779', 'MarkerFaceColor', '#35b779', 'Marker', ...
    'pentagram')
%E
plot(HRR_kw(1:length(Heptane_pool_uncert_frac.gas_cons(:,3))),...
    Heptane_pool_uncert_frac.OCC_cons(:,2)*100, ...
    'linewidth', 0.75, 'Color', '#31688e', 'LineStyle','-')
scatter(HRR_kw(1:10:length(Heptane_pool_uncert_frac.cal_gas(:,3))),...
    Heptane_pool_uncert_frac.OCC_cons(1:10:end,2)*100, 6, ...
    'MarkerEdgeColor', '#31688e', 'MarkerFaceColor', '#31688e', 'Marker', ...
    's')
% e CO
plot(HRR_kw(1:length(Heptane_pool_uncert_frac.gas_cons(:,3))),...
    Heptane_pool_uncert_frac.misc_CO(:,1)*100*3, ...
    'linewidth', 0.75, 'Color', '#440154', 'LineStyle','-')
scatter(HRR_kw(1:10:length(Heptane_pool_uncert_frac.cal_gas(:,3))),...
    Heptane_pool_uncert_frac.misc_CO(1:10:end,1)*100*3, 6, ...
    'MarkerEdgeColor', '#440154', 'MarkerFaceColor', '#440154', 'Marker', ...
    '^')
% label
title('Heptane', 'FontSize', 8, 'FontWeight','bold')
%text(20, 2.32, 'z_{hood}', 'FontSize', 8)

%specifs axis numbers size and font
ax = gca;
ax.FontSize = 8;
ax.FontName = 'Times';
ax.TitleHorizontalAlignment = 'center';
%ax.YAxis.TickValues = [10e-5 10e-4 10e-3 10e-2 10e-1 10e-0, 10e1, 10e2];
ax.YAxis.TickLabel = [];

lg = legend({ '', 'X_{O_2}^{span}', '', 'f(Re)', '', 'E', '', ...
    '\epsilon _{CO, comb}'}, FontSize=8, Box="off");
lg.Layout.Tile = 'East';
%lg.Location = 'northeast';
%lg.NumColumns = 2;

%Save the figure as a pdf
set(gcf,'PaperUnits','inches','PaperPosition', [0 0 4.75 2], 'PaperSize', [4.75 2]);
print -dpdf Figures/Uncertainty_fraction.pdf -r150

% % %% Do analysis without CO2 correction
% % Propane_burner.CO2_zero_offset_percent = 0;
% % Heptane_pool.CO2_zero_offset_percent = 0;
% % % run the anlysis
% % % Calculate HRR and uncertainty given model data considering CO and CO2
% % [Propane_burner_uncert_cross_sens_false, Propane_burner_HRR_calc_kW_cross_sens_false, ...
% %     ~] = uncert_HRR(Propane_burner_uncert, ...
% %     Propane_burner, 'True');
% % % Calculate HRR and uncertainty given model data considering CO and CO2
% % [Heptane_pool_uncert_cross_sens_false, Heptane_pool_HRR_calc_kW_cross_sens_false, ...
% %     ~] = uncert_HRR(Heptane_pool_uncert, ...
% %     Heptane_pool, 'True');
% % 
% % % make the plot
% % figure
% % hold on
% % box on
% % 
% % xlabel('HRR (kW)', 'FontName', 'Times', 'FontSize', 8,'FontWeight', 'bold')
% % ylabel('$U(\dot{Q})/\dot{Q}$ \bf{(\%)}', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 8, ...
% %     'FontWeight', 'bold')
% % 
% % %limit the x-axis
% % xlim([0 750])
% % ylim([0 20])
% % 
% % 
% % % O2 only
% % plot(HRR_kw(1:10:end), Propane_burner_uncert_cross_sens_false.HRR_kW(1:10:end)./...
% %     Propane_burner_HRR_calc_kW_cross_sens_false(1:10:end)*100, 'linewidth', 0.75, 'Color', ...
% %     '#5ec962', 'LineStyle','-', 'Marker', 'o', 'MarkerSize', 3)
% % % The fraction of the heat realease rate 
% % plot(HRR_kw, Propane_burner_uncert_with.HRR_kW./HRR_kw*100, 'linewidth', 0.75, 'Color', ...
% %     '#5ec962', 'LineStyle','-')
% % 
% % 
% % plot(HRR_kw(1:10:end), Heptane_pool_uncert_cross_sens_false.HRR_kW(1:10:end)./...
% %     Heptane_pool_HRR_calc_kW_cross_sens_false(1:10:end)*100, 'linewidth', 0.75, 'Color', ...
% %     '#3b528b', 'LineStyle','--', 'Marker', '^', 'MarkerSize', 3)
% % plot(HRR_kw, Heptane_pool_uncert_with.HRR_kW./HRR_kw*100, 'linewidth', 0.75, 'Color', ...
% %     '#3b528b', 'LineStyle','--')
% % 
% % % label
% % title('b)', 'FontSize', 8, 'FontWeight','bold')
% % %text(20, 2.32, 'z_{hood}', 'FontSize', 8)
% % 
% % %specifs axis numbers size and font
% % ax = gca;
% % ax.FontSize = 8;
% % ax.FontName = 'Times';
% % ax.TitleHorizontalAlignment = 'left';
% % 
% % lg = legend({ 'Propane (no correction)', 'Propane ($\dot{Q}_2$)', ...
% %     'Heptane (no correction)', 'Heptane ($\dot{Q}_2$)'}, ...
% %     'Interpreter', 'latex', FontSize=8, Box="off");
% % % lg.Layout.Tile = 'South';
% % lg.Location = 'southeast';
% % % lg.NumColumns = 2;
% % 
% % %Save the figure as a pdf
% % set(gcf,'PaperUnits','inches','PaperPosition', [0 0 4 2], 'PaperSize', [4 2]);
% % print -dpdf Figures/HRR_uncert_comp_CO2_cross_sens.pdf -r150
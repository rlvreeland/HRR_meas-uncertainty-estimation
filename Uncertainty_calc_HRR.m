function [HRR, grad_cal_gas, grad_misc_CO, grad_misc_CO2, grad_misc_O2,...
    grad_gas_cons, grad_OCC_cons, grad_flow_vars] = ...
    Uncertainty_calc_HRR(cal_gas, misc_CO, misc_CO2, misc_O2, gas_cons, ...
    OCC_cons, flow_vars, conc_extra, include_CO_CO2)
% This function calcuculates the HRR and the partial derivates with respect
% to each varable. The CO and CO2 deterimines if CO and CO2 are inculded in
% the HRR calculation

% Please note that the variables for the function are combined in mulitple
% dl arrays. The arrays are as follows:
% calibration_gas = [x_CO_bottle_span, x_CO2_bottle_span, x_O2_bottle_span,
%                    x_i_bottle_zero]
% misc_CO = [intresic, linearity, repeatability];
% misc_CO2 = [intresic, linearity, repeatability];
% misc_O2 = [flow, line drop, T, P_baro, drift, absorption, HC_cross_sens, 
%         CO2_zero_offset_percent]
% gas_cons = [R, P_amb, Mw_air, Mw_O2 x_CO2_amb, x_O2_amb]
% OCC_cons = [alpha, E_O2, E_CO2]
% flow_vars = [f(Re) for the BDP, Delta_P, Te, D_duct]
% conc_extra = given values of CO, CO2, CO2, then the span range for CO,
% CO2, O2, and the zero value.

%% Extract each parameter out of each array to inprove code redability
% Calibration Gases
x_CO_bottle_span = cal_gas(1);
x_CO2_bottle_span = cal_gas(2);
x_O2_bottle_span = cal_gas(3);
x_i_bottle_zero = cal_gas(4);

% CO concentraion miscienous parameters
intresic_CO = misc_CO(1);
linearity_CO = misc_CO(2);
repeatability_CO = misc_CO(3);

% CO2 concentraion miscienous parameters
intresic_CO2 = misc_CO2(1);
linearity_CO2 = misc_CO2(2);
repeatability_CO2 = misc_CO2(3);

% O2 concentration miscienous parameters
%   please note that all except the CO2 offset are just added at the end of
%   the O2 calculation so to not have to write out each paramter the first
%   seven uncertainty sources
misc_uncert_sorces_O2 = sum(misc_O2(1:7), 'all');
CO2_zero_offset_percent = misc_O2(8);

% Gas constants
R_m3_Pa_per_K_mol = gas_cons(1);
P_amb_Pa = gas_cons(2);
Mw_air_g_per_mol = gas_cons(3);
Mw_O2_g_per_mol = gas_cons(4);
x_CO2_amb_given = gas_cons(5);
x_O2_amb_given = gas_cons(6);

% Oxygen consumption calorimtery constatnts
alpha = OCC_cons(1);
E_kJ_per_kg = OCC_cons(2);
E_CO_kJ_per_kg = OCC_cons(3);

% Flow variables
f_Re = flow_vars(1);
Delta_P_Pa = flow_vars(2);
Te_K = flow_vars(3);
D_duct_m = flow_vars(4);

% Extra concentration paramters
x_CO_given = conc_extra(1);
x_CO2_given = conc_extra(2);
x_O2_given = conc_extra(3);
span_range_CO = conc_extra(4);
span_range_CO2 = conc_extra(5);
span_range_O2 = conc_extra(6);
zero = conc_extra(7);

%% Calcuate concentration of gas species
x_CO_calc = (x_CO_bottle_span-x_i_bottle_zero)/(span_range_CO)*(x_CO_given-zero)...
    +intresic_CO + linearity_CO + repeatability_CO;

x_CO2_calc = (x_CO2_bottle_span-x_i_bottle_zero)/(span_range_CO2)*(x_CO2_given-zero)...
    +intresic_CO2 + linearity_CO2 + repeatability_CO2;
x_CO2_amb = (x_CO2_bottle_span-x_i_bottle_zero)/(span_range_CO2)*(x_CO2_amb_given-zero)...
    +intresic_CO2 + linearity_CO2 + repeatability_CO2;
% for x_O2 have to use the x_CO2_calc so that the uncertantiy in CO2 concentration
% is accounted for
x_O2_calc = (x_O2_bottle_span-x_i_bottle_zero)/(span_range_O2)*...
    (x_O2_given-zero) + misc_uncert_sorces_O2 - x_CO2_calc* ...
    CO2_zero_offset_percent/100;
x_O2_amb = (x_O2_bottle_span-x_i_bottle_zero)/(span_range_O2)*...
    (x_O2_amb_given-zero) + misc_uncert_sorces_O2 - x_CO2_amb* ...
    CO2_zero_offset_percent/100;

%% Calculate the mass flow in the duct
rho_kg_per_m3 = P_amb_Pa*Mw_air_g_per_mol*10^(-3)/(R_m3_Pa_per_K_mol*Te_K);
vel_m_per_s = ((2*Delta_P_Pa/rho_kg_per_m3)^0.5)/f_Re;
A_duct_m2 = pi*(D_duct_m/2)^2;
m_e_kg_per_s = vel_m_per_s*A_duct_m2*rho_kg_per_m3;

%% Calculate the HRR in kW
% this deterimines if CO and CO2 are used or not
if contains(include_CO_CO2, 'False')
    % % % O2 only
    % % phi = (x_O2_amb-x_O2_calc)/((1-x_O2_calc)*x_O2_amb);
    % % HRR = E_kJ_per_kg*phi/(1+phi*(alpha-1))*...
    % %         Mw_O2_g_per_mol/Mw_air_g_per_mol*m_e_kg_per_s*x_O2_amb;
      % CO2 and O2 only
    phi = (x_O2_amb-x_O2_calc)/((1-x_O2_calc)*x_O2_amb);
    HRR = E_kJ_per_kg*phi/(1+phi*(alpha-1))*...
            Mw_O2_g_per_mol/Mw_air_g_per_mol*m_e_kg_per_s*x_O2_amb;
else % signifies that string contains true
    % CO, CO2, O2
     phi = (x_O2_amb*(1-x_CO2_calc-x_CO_calc)-x_O2_calc*(1-x_CO2_amb))/...
    ((1-x_O2_calc-x_CO2_calc-x_CO_calc)*x_O2_amb);
    HRR = E_kJ_per_kg*phi/(1+phi*(alpha-1))*...
    Mw_O2_g_per_mol/Mw_air_g_per_mol*m_e_kg_per_s*x_O2_amb;
end

%% Calculate the partial derivates for each of the dlarrays used
grad_cal_gas = dlgradient(HRR, cal_gas);
grad_misc_CO = dlgradient(HRR, misc_CO);
grad_misc_CO2 = dlgradient(HRR, misc_CO2);
grad_misc_O2 = dlgradient(HRR, misc_O2);
grad_gas_cons = dlgradient(HRR, gas_cons);
grad_OCC_cons = dlgradient(HRR, OCC_cons);
grad_flow_vars = dlgradient(HRR, flow_vars);

end
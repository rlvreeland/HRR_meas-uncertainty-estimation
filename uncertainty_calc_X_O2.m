function [x_O2_calc, grad_cal_gas, grad_misc_O2, grad_misc_CO2] = ...
    uncertainty_calc_X_O2(cal_gas, misc_O2, misc_CO2, x_O2_given, ...
    x_CO2_given, span_range_O2, span_range_CO2, zero)
% This function is used to calculate partial derviatives if the input 
% parameters for O2, when voltages are not available given the measurement
% concentration and the calibration gases. This fuction also
% includes a term to correct for the effects of CO2 concentration.

% Please note that the variables for the function are combined in mulitple
% dl arrays. The arrays are as follows:
% calibration_gas = [x_CO_bottle_span, x_CO2_bottle_span, x_O2_bottle_span,
%                    x_i_bottle_zero]
% misc_O2 = [flow, line drop, T, P_baro, drift, absorption, HC_cross_sens, 
%         CO2_zero_offset_percent]
% misc_CO2 = [intresic, linearity, repeatability];


% Extract each parameter out of each array to inprove code redability
x_CO2_bottle_span = cal_gas(2);
x_O2_bottle_span = cal_gas(3);
x_i_bottle_zero = cal_gas(4);


intresic = misc_CO2(1);
linearity = misc_CO2(2);
repeatability = misc_CO2(3);

% please note that all except the CO2 offset are just added at the end of
% the O2 calculation so to not have to write out each paramter the first
% seven uncertainty sources
misc_uncert_sorces_O2 = sum(misc_O2(1:7), 'all');
CO2_zero_offset_percent = misc_O2(8);

% Calculate the concentration of CO2
x_CO2_calc = (x_CO2_bottle_span-x_i_bottle_zero)/(span_range_CO2)*(x_CO2_given-zero)...
    +intresic + linearity + repeatability;

%Calcualte the concentration of O2
% have to use the x_CO2_calc so that the uncertantiy in CO2 concentration
% is accounted for
x_O2_calc = (x_O2_bottle_span-x_i_bottle_zero)/(span_range_O2)*...
    (x_O2_given-zero) + misc_uncert_sorces_O2 - x_CO2_calc* ...
    CO2_zero_offset_percent/100;

grad_cal_gas = dlgradient(x_O2_calc, cal_gas);
grad_misc_O2 = dlgradient(x_O2_calc, misc_O2);
grad_misc_CO2 = dlgradient(x_O2_calc, misc_CO2);

end
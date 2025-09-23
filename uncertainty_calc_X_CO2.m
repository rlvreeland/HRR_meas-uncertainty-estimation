function [x_CO2_calc, grad_cal_gas, grad_misc] = ...
    uncertainty_calc_X_CO2(cal_gas, misc, x_CO2_given, span_range, zero)
% This function is used to calculate partial derviatives if the input 
% parameters for CO2, when voltages are not available given the measurement
% concentration and the calibration gases.

% Please note that the variables for the function are combined in mulitple
% dl arrays. The arrays are as follows:
% calibration_gas = [x_CO_bottle_span, x_CO2_bottle_span, x_O2_bottle_span,
%                    x_i_bottle_zero]
% misc = [intresic, linearity, repeatability];

% Extract each parameter out of each array to inprove code redability
x_CO2_bottle_span = cal_gas(2);
x_i_bottle_zero = cal_gas(4);


intresic = misc(1);
linearity = misc(2);
repeatability = misc(3);

% Calculate the concentration of CO2
x_CO2_calc = (x_CO2_bottle_span-x_i_bottle_zero)/(span_range)*(x_CO2_given-zero)...
    +intresic + linearity + repeatability;



grad_cal_gas = dlgradient(x_CO2_calc, cal_gas);
grad_misc = dlgradient(x_CO2_calc, misc);


end
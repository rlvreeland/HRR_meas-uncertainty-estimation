function [x_CO_calc, grad_cal_gas, grad_misc] = ...
    uncertainty_calc_X_CO(cal_gas, misc, x_CO_given, span_range, zero)
% This function is used to calculate partial derivatives of the input 
% parameters for CO, when voltages are not available given the measurement
% concentration and the calibration gases.

% Please note that the variables for the function are combined in multiple
% dl arrays. The arrays are as follows:
% calibration_gas = [x_CO_bottle_span, x_CO2_bottle_span, x_O2_bottle_span,
%                    x_i_bottle_zero]
% misc = [intresic, linearity, repeatability];

% Extract each parameter out of each array to improve code readability
x_CO_bottle_span = cal_gas(1);
x_i_bottle_zero = cal_gas(4);

intresic = misc(1);
linearity = misc(2);
repeatability = misc(3);

% Calculate the concentration of CO
x_CO_calc = (x_CO_bottle_span-x_i_bottle_zero)/(span_range)*(x_CO_given-zero)...
    +intresic + linearity + repeatability;


grad_cal_gas = dlgradient(x_CO_calc, cal_gas);
grad_misc = dlgradient(x_CO_calc, misc);


end
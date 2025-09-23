function [m_e_kg_per_s, grad_flow_vars, grad_gas_cons] = ...
    uncert_mass_flow_measure(flow_vars, gas_cons)
% This function is used to calculate the partial derivatives for mass flow
% measurment in the duct.

% Please note that the variables for the function are combined in mulitple
% dlarray. The arrays is as follows:
% flow_vars = [f(Re) for the BDP, Delta_P, Te, D_duct]
% gas_cons = [R, P_amb, Mw]

%% Extract needed parameters out of each array to inprove code redability
% flow varables
f_Re = flow_vars(1);
Delta_P_Pa = flow_vars(2);
Te_K = flow_vars(3);
D_duct_m = flow_vars(4);

% Gas constants
R_m3_Pa_per_K_mol = gas_cons(1);
P_amb_Pa = gas_cons(2);
Mw_air_g_per_mol = gas_cons(3);


%% Calculate mass flow in the duct
rho_kg_per_m3 = P_amb_Pa*Mw_air_g_per_mol*10^(-3)/(R_m3_Pa_per_K_mol*Te_K);
vel_m_per_s = (2*Delta_P_Pa/rho_kg_per_m3)^0.5/f_Re;
A_duct_m2 = pi*(D_duct_m/2)^2;
m_e_kg_per_s = vel_m_per_s*A_duct_m2*rho_kg_per_m3;

%% Calculate partial derivates
grad_flow_vars = dlgradient(m_e_kg_per_s, flow_vars);
grad_gas_cons = dlgradient(m_e_kg_per_s, gas_cons);

end
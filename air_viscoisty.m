function [visco_Pa_s] = air_viscoisty(rho_kg_per_m3,T_K)
%This function calculates the dynamic viscosisty (eta) of air given the temperature and
%density of the air.

% Define neccarray paramters for corralation
T_star_K = 132.5;
rho_star_kg_m3 = 314.3;
H_Pa_s = 6.1609*10^-6;
A_1 = 0.128517;
A_0_point_5 = 2.60661;

A_array = [-1, -0.709661, 0.662534, -0.197846, 0.00770147];
A_index = [0, -1, -2, -3, -4];

B_array = [0.465601, 1.26469, -0.511425, 0.2746]; 
B_index = [1, 2, 3, 4];

% apply the correlation
T_r = T_K/T_star_K;
rho_r = rho_kg_per_m3/rho_star_kg_m3;

delta_eta = sum(B_array.*rho_r.^B_index);
eta_0 = A_1*T_r+A_0_point_5*T_r^0.5+ sum(A_array.*T_r.^A_index);
visco_Pa_s = H_Pa_s*(eta_0 + delta_eta);

end
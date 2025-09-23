classdef measurement_uncert
    properties
        % Flow measurement
        f_Re = 0.015;
        Delta_P_Pa
        Te_K

        % Calibration bottles
        x_CO_bottle_span {mustBeNumeric}
        x_CO2_bottle_span {mustBeNumeric}
        x_O2_bottle_span {mustBeNumeric}
        x_i_bottle_zero {mustBeNumeric}

        %CO sensor parameters
        Intrinsic_error_CO = 0.01*5/100;
        Linearity_error_CO = 0.01*5/100;
        Reapeatability_CO = 0.01*5/100;

        %CO2 sensor parameters
        Intrinsic_error_CO2 = 0.01*10/100;
        Linearity_error_CO2 = 0.01*10/100;
        Reapeatability_CO2 = 0.01*10/100;

        %O2 sensor parameters that contribute to the uncertainty in the
        %measurement. Each show the effect in mol fraction.
        O2_sens_flow = 5*10^(-6);
        O2_sens_exhaust_line_P_drop = 8*10^(-6);
        O2_sens_T_short_term = 2*10^(-6);
        O2_sens_P_baro_short_term = 200*10^(-6);
        O2_sens_noise_drift = 20*10^(-6);
        O2_sens_adsorption = 10*10^(-6);
        O2_sesn_cross_sens_HC = 1*10^(-6);
        CO2_zero_offset_percent = 0; % This code does not account for the
        %                               uncertainty in the offest in
        %                               measured O2 concentration from CO2,
        %                               since uncertainty is set to zero,
        %                               but will be accounted for if this
        %                               parameter is modified.


        % Calorimetery
        alpha = 0.105/(3^0.5); % This converts to standard error
        % assuming a square distibution
        E_kJ_per_kg = 655;
        E_CO_kJ_per_kg = 17.6*10^3*0.00057;

        % Calorimeter specfications
        D_duct_m {mustBeNumeric}

        % Gas constant and ambient conditions
        R_m3_Pa_per_K_mol = 0.00026;
        P_amb_Pa = 101325*0.0015;
        T_amb_K = 2.2;
        x_O2_amb
        x_CO2_amb
        x_CO_amb

        Mw_air_g_per_mol = 0.01;
        Mw_O2_g_per_mol = 0;

        % Calculated uncertaintes
        x_CO
        x_CO2
        x_O2
        m_e_kg_per_s
        HRR_kW

    end

    methods
        function [obj] = uncert_Te_K_k_type(obj, meas_Te_K)
            % calculate the uncertainty in Temperature measurements
            % assuming a regular k-type TC is used
            %preallocate
            uncert_Te_k = NaN(size(meas_Te_K));
            for i = 1:length(meas_Te_K)
                if (meas_Te_K(i)-273.15)*0.0075 > 2.2
                    uncert_Te_k(i) = (meas_Te_K(i)-273.15)*0.0075;
                else
                    uncert_Te_k(i) = 2.2;
                end
            end
            % added to object
            obj.Te_K = uncert_Te_k;
        end

        function [obj_uncert] = ...
                uncert_conc_measure(obj_uncert, obj_fire_model)
            % Create dl array for calibration gases
            cal_gas_dl = dlarray([obj_fire_model.x_CO_bottle_span, ...
                obj_fire_model.x_CO2_bottle_span, obj_fire_model.x_O2_bottle_span, ...
                obj_fire_model.x_i_bottle_zero]);
            % Create matching array with the uncertainties of each
            % calibation gas
            uncert_cal_gas_array = [obj_uncert.x_CO_bottle_span, ...
                obj_uncert.x_CO2_bottle_span, obj_uncert.x_O2_bottle_span, ...
                obj_uncert.x_i_bottle_zero];

            % Create dl array for misc parameters for each species
            % Note these are place holder values with just uncertainty
            % assigned
            misc_CO_dl = dlarray([0, 0, 0]); %[intresic, linearity, repeatability]
            misc_CO2_dl = dlarray([0, 0, 0]); %[intresic, linearity, repeatability]
            misc_O2_dl = dlarray([0, 0, 0, 0, 0, 0, 0, obj_fire_model.CO2_zero_offset_percent]);
            % [flow, line drop, T, P_baro, drift, absorption,
            % HC_cross_sens, CO2_zero_offset_percent]
            % Create matching uncertainty arrays
            uncert_misc_CO_array = [obj_uncert.Intrinsic_error_CO, ...
                obj_uncert.Linearity_error_CO, obj_uncert.Reapeatability_CO];
            uncert_misc_CO2_array = [obj_uncert.Intrinsic_error_CO2, ...
                obj_uncert.Linearity_error_CO2, obj_uncert.Reapeatability_CO2];
            uncert_misc_O2_array = [obj_uncert.O2_sens_flow, ...
                obj_uncert.O2_sens_exhaust_line_P_drop, ...
                obj_uncert.O2_sens_T_short_term, ...
                obj_uncert.O2_sens_P_baro_short_term, ...
                obj_uncert.O2_sens_noise_drift, obj_uncert.O2_sens_adsorption,...
                obj_uncert.O2_sesn_cross_sens_HC, obj_uncert.CO2_zero_offset_percent];

            % apply this operation for each HRR indvidually
            %preallocte
            obj_uncert.x_CO = NaN(size(obj_fire_model.x_CO_measure));
            obj_uncert.x_CO2 = NaN(size(obj_fire_model.x_CO_measure));
            obj_uncert.x_O2 = NaN(size(obj_fire_model.x_CO_measure));

            for i = 1:(length(obj_fire_model.x_CO_measure)+1)
                if i == length(obj_fire_model.x_CO_measure)+1
                    conc_CO2 = obj_fire_model.x_CO2_amb;
                    conc_CO = obj_fire_model.x_CO_amb;
                    conc_O2 = obj_fire_model.x_O2_amb;
                else
                    conc_CO2 = obj_fire_model.x_CO2_measure(i);
                    conc_CO = obj_fire_model.x_CO_measure(i);
                    conc_O2 = obj_fire_model.x_O2_measure(i);
                end
                % check if NaNs
                if isnan(conc_O2)
                else
                    %%% X_CO Ucertainty
                    % Get partials as dl arrays
                    [~, partials_cal_gas_dl, partials_misc_dl] = ...
                        dlfeval(@uncertainty_calc_X_CO, cal_gas_dl, misc_CO_dl, ...
                        dlarray(conc_CO), ...
                        dlarray(obj_fire_model.x_CO_bottle_span...
                        -obj_fire_model.x_i_bottle_zero), ...
                        dlarray(obj_fire_model.x_i_bottle_zero));
                    % Convert to regular array
                    partials_cal_gas_CO = extractdata(partials_cal_gas_dl);
                    partials_misc_CO = extractdata(partials_misc_dl);
                    % Calculate the uncertainty
                    uncert_conc = (sum((partials_cal_gas_CO.*uncert_cal_gas_array)...
                        .^2, 'all')+ sum((partials_misc_CO.*uncert_misc_CO_array)...
                        .^2, 'all'))^0.5;
                    % Place result in correct location
                    if i == length(obj_fire_model.x_CO_measure)+1
                        obj_uncert.x_CO_amb = uncert_conc;
                    else
                        obj_uncert.x_CO(i) = uncert_conc;
                    end

                    %%% X_CO2 Ucertainty
                    % Get partials as dl arrays
                    [check, partials_cal_gas_dl, partials_misc_dl] = ...
                        dlfeval(@uncertainty_calc_X_CO2, cal_gas_dl, misc_CO2_dl, ...
                        dlarray(conc_CO2), ...
                        dlarray(obj_fire_model.x_CO2_bottle_span...
                        -obj_fire_model.x_i_bottle_zero), ...
                        dlarray(obj_fire_model.x_i_bottle_zero));
                    % Convert to regular array
                    partials_cal_gas_CO2 = extractdata(partials_cal_gas_dl);
                    partials_misc_CO2 = extractdata(partials_misc_dl);
                    % Calculate the uncertainty
                    uncert_conc = (sum((partials_cal_gas_CO2.*uncert_cal_gas_array)...
                        .^2, 'all')+ sum((partials_misc_CO2.*uncert_misc_CO2_array)...
                        .^2, 'all'))^0.5;
                    % Place result in correct location
                    if i == length(obj_fire_model.x_CO_measure)+1
                        obj_uncert.x_CO2_amb = uncert_conc;
                    else
                        obj_uncert.x_CO2(i) = uncert_conc;
                    end

                    %%% X_O2 Ucertainty Considering x_CO2
                    % Get partials as dl arrays
                    [check, partials_cal_gas_dl, partials_misc_O2_dl, ...
                        partials_misc_CO2_dl] = dlfeval(@uncertainty_calc_X_O2, ...
                        cal_gas_dl, misc_O2_dl, misc_CO2_dl, ...
                        dlarray(conc_O2), dlarray(conc_CO2),...
                        dlarray(obj_fire_model.x_O2_bottle_span...
                        -obj_fire_model.x_i_bottle_zero), ...
                        dlarray(obj_fire_model.x_CO2_bottle_span...
                        -obj_fire_model.x_i_bottle_zero), ...
                        dlarray(obj_fire_model.x_i_bottle_zero));
                    % Convert to regular array
                    partials_cal_gas = extractdata(partials_cal_gas_dl);
                    partials_misc_O2 = extractdata(partials_misc_O2_dl);
                    partials_misc_CO2 = extractdata(partials_misc_CO2_dl);
                    % Calculate the uncertainty
                    uncert_conc = (sum((partials_cal_gas.*uncert_cal_gas_array)...
                        .^2, 'all')+ sum((partials_misc_O2.*uncert_misc_O2_array)...
                        .^2, 'all')+ sum((partials_misc_CO2.*uncert_misc_CO2_array)...
                        .^2, 'all'))^0.5;

                    % Place result in correct location
                    if i == length(obj_fire_model.x_CO_measure)+1
                        obj_uncert.x_O2_amb = uncert_conc;
                    else
                        obj_uncert.x_O2(i) = uncert_conc;
                    end

                    % % %%% X_O2 Ucertainty without accounting for CO2
                    % %    % Get partials as dl arrays
                    % %    [~, partials_cal_gas_dl, partials_misc_O2_dl] = ...
                    % %        dlfeval(@uncertainty_calc_X_O2_no_CO2_account, ...
                    % %        cal_gas_dl, misc_O2_dl, ...
                    % %        dlarray(conc_O2), ...
                    % %        dlarray(obj_fire_model.x_O2_bottle_span...
                    % %        -obj_fire_model.x_i_bottle_zero), ...
                    % %        dlarray(obj_fire_model.x_i_bottle_zero));
                    % %    % Convert to regular array
                    % %    partials_cal_gas = extractdata(partials_cal_gas_dl);
                    % %    partials_misc_O2 = extractdata(partials_misc_O2_dl);
                    % %    % Calculate the uncertainty
                    % %    uncert_conc = (sum((partials_cal_gas.*uncert_cal_gas_array)...
                    % %        .^2, 'all')+ sum((partials_misc_O2.*uncert_misc_O2_array)...
                    % %        .^2, 'all'))^0.5;
                    % Place result in correct location
                    % % if i == length(obj_fire_model.x_CO_measure)+1
                    % %     obj_uncert.x_O2_amb = uncert_conc;
                    % % else
                    % %     obj_uncert.x_O2(i) = uncert_conc;
                    % % end
                end
            end
        end

        function [obj_uncert] = ...
                uncert_mass_flow_measure(obj_uncert, obj_fire_model)
            %Create dlarray for gas constants
            gas_cons_dl = dlarray([obj_fire_model.R_m3_Pa_per_K_mol,...
                obj_fire_model.P_amb_Pa, obj_fire_model.Mw_air_g_per_mol]);
            % Create matching uncertainty array
            uncert_gas_cons_array = [obj_uncert.R_m3_Pa_per_K_mol,...
                obj_uncert.P_amb_Pa, obj_uncert.Mw_air_g_per_mol];
            % Complete assment for each HRR indvidually
            % Preallocate
            obj_uncert.m_e_kg_per_s = NaN(size(obj_fire_model.x_CO_measure));
            for i = 1:(length(obj_fire_model.x_CO_measure))
                if isnan(obj_fire_model.x_O2_measure(i))
                else
                    % Create dlarray for flow varables
                    flow_vars_dl = dlarray([obj_fire_model.f_Re,...
                        obj_fire_model.Delta_P_Pa(i), ...
                        obj_fire_model.Te_K(i), obj_fire_model.D_duct_m]);
                    % Create matching uncertainty array
                    uncert_flow_vars_array = [obj_uncert.f_Re,...
                        obj_uncert.Delta_P_Pa(i), ...
                        obj_uncert.Te_K(i), obj_uncert.D_duct_m];
                    % run function to calculate partials for m_e
                    [~, partials_flow_vars_dl, partials_gas_cons_dl] = ...
                        dlfeval(@uncert_mass_flow_measure, flow_vars_dl,...
                        gas_cons_dl);
                    % Convert to regualar array
                    partials_flow_vars = extractdata(partials_flow_vars_dl);
                    partials_gas_cons = extractdata(partials_gas_cons_dl);

                    % calcualate the uncertainty
                    obj_uncert.m_e_kg_per_s(i) = (sum((partials_flow_vars...
                        .*uncert_flow_vars_array).^2, 'all')+ ...
                        sum((partials_gas_cons.*uncert_gas_cons_array)...
                        .^2, 'all'))^0.5;
                end
            end

        end

        function [obj_uncert, HRR_calc, ds_uncert_fraction] = uncert_HRR(obj_uncert, ...
                obj_fire_model, include_CO_CO2)
            arguments
                obj_uncert 
                obj_fire_model 
                include_CO_CO2 (1,:) char {mustBeMember(include_CO_CO2, ...
                    {'True', 'False'})} 
            end
            % Create dl array for calibration gases
            cal_gas_dl = dlarray([obj_fire_model.x_CO_bottle_span, ...
                obj_fire_model.x_CO2_bottle_span, obj_fire_model.x_O2_bottle_span, ...
                obj_fire_model.x_i_bottle_zero]);
            % Create matching array with the uncertainties of each
            % calibation gas
            uncert_cal_gas_array = [obj_uncert.x_CO_bottle_span, ...
                obj_uncert.x_CO2_bottle_span, obj_uncert.x_O2_bottle_span, ...
                obj_uncert.x_i_bottle_zero];

            % Create dl array for misc parameters for each species
            % Note these are place holder values with just uncertainty
            % assigned
            misc_CO_dl = dlarray([0, 0, 0]); %[intresic, linearity, repeatability]
            misc_CO2_dl = dlarray([0, 0, 0]); %[intresic, linearity, repeatability]
            misc_O2_dl = dlarray([0, 0, 0, 0, 0, 0, 0, obj_fire_model.CO2_zero_offset_percent]);
            % [flow, line drop, T, P_baro, drift, absorption,
            % HC_cross_sens, CO2_zero_offset_percent]
            % Create matching uncertainty arrays
            uncert_misc_CO_array = [obj_uncert.Intrinsic_error_CO, ...
                obj_uncert.Linearity_error_CO, obj_uncert.Reapeatability_CO];
            uncert_misc_CO2_array = [obj_uncert.Intrinsic_error_CO2, ...
                obj_uncert.Linearity_error_CO2, obj_uncert.Reapeatability_CO2];
            uncert_misc_O2_array = [obj_uncert.O2_sens_flow, ...
                obj_uncert.O2_sens_exhaust_line_P_drop, ...
                obj_uncert.O2_sens_T_short_term, ...
                obj_uncert.O2_sens_P_baro_short_term, ...
                obj_uncert.O2_sens_noise_drift, obj_uncert.O2_sens_adsorption,...
                obj_uncert.O2_sesn_cross_sens_HC, obj_uncert.CO2_zero_offset_percent];
            %Create dlarray for gas constants include ambient conditions
            gas_cons_dl = dlarray([obj_fire_model.R_m3_Pa_per_K_mol,...
                obj_fire_model.P_amb_Pa, obj_fire_model.Mw_air_g_per_mol, ...
                obj_fire_model.Mw_O2_g_per_mol, obj_fire_model.x_CO2_amb, ...
                obj_fire_model.x_O2_amb]);
            % Create matching uncertainty array
            uncert_gas_cons_array = [obj_uncert.R_m3_Pa_per_K_mol,...
                obj_uncert.P_amb_Pa, obj_uncert.Mw_air_g_per_mol,  ...
                obj_uncert.Mw_O2_g_per_mol, 0, ...
                0]; 
            % Create dlarray for Oxgyen cosumption calorimetery constants
            OCC_cons_dl = dlarray([obj_fire_model.alpha, ...
                obj_fire_model.E_kJ_per_kg, obj_fire_model.E_CO_kJ_per_kg]);
            % Create matching uncertainty array
            uncert_OCC_cons_array = [obj_uncert.alpha, ...
                obj_uncert.E_kJ_per_kg, obj_uncert.E_CO_kJ_per_kg];

            % Complete assment for each HRR indvidually
            % Preallocate
            obj_uncert.HRR_kW = NaN(size(obj_fire_model.x_CO_measure));
            HRR_calc = NaN(size(obj_fire_model.x_CO_measure));
            for i = 1:(length(obj_fire_model.x_CO_measure))
                if isnan(obj_fire_model.x_O2_measure(i))
                else
                    % Create dlarray for flow varables
                    flow_vars_dl = dlarray([obj_fire_model.f_Re,...
                        obj_fire_model.Delta_P_Pa(i), ...
                        obj_fire_model.Te_K(i), obj_fire_model.D_duct_m]);
                    % Create matching uncertainty array
                    uncert_flow_vars_array = [obj_uncert.f_Re,...
                        obj_uncert.Delta_P_Pa(i), ...
                        obj_uncert.Te_K(i), obj_uncert.D_duct_m];
                    % Make dlarray to hold extra paramters for
                    % concentration calculation that have no assigned
                    % uncertainty.
                    Conc_extra_dl = dlarray([obj_fire_model.x_CO_measure(i), ...
                        obj_fire_model.x_CO2_measure(i),...
                        obj_fire_model.x_O2_measure(i), ...
                        obj_fire_model.x_CO_bottle_span-obj_fire_model.x_i_bottle_zero, ...
                        obj_fire_model.x_CO2_bottle_span-obj_fire_model.x_i_bottle_zero, ...
                        obj_fire_model.x_O2_bottle_span-obj_fire_model.x_i_bottle_zero, ...
                        obj_fire_model.x_i_bottle_zero]);
                % get paritals and HRR
                [HRR_calc_dl, partial_cal_gas_dl, partial_misc_CO_dl, ...
                    partial_misc_CO2_dl, partial_misc_O2_dl, partial_gas_cons_dl, ...
                    partial_OCC_cons_dl, parital_flow_vars_dl] = ...
                    dlfeval(@Uncertainty_calc_HRR, cal_gas_dl, misc_CO_dl,...
                    misc_CO2_dl, misc_O2_dl, gas_cons_dl, OCC_cons_dl, ...
                    flow_vars_dl, Conc_extra_dl, include_CO_CO2);

                % convert HRR to regular array
                HRR_calc(i) = extractdata(HRR_calc_dl);
                % convert to regular arrays and calcuate uncertiant
                % components^2
                delta_part2_cal_gas = (extractdata(partial_cal_gas_dl).* ...
                    uncert_cal_gas_array).^2;
                delta_part2_misc_CO = (extractdata(partial_misc_CO_dl) .* ...
                    uncert_misc_CO_array).^2;
                delta_part2_misc_CO2 = (extractdata(partial_misc_CO2_dl).* ...
                    uncert_misc_CO2_array).^2;
                delta_part2_misc_O2 = (extractdata(partial_misc_O2_dl) .* ...
                    uncert_misc_O2_array).^2;
                delta_part2_gas_cons = (extractdata(partial_gas_cons_dl) .* ...
                    uncert_gas_cons_array).^2;
                delta_part2_OCC_cons = (extractdata(partial_OCC_cons_dl) .* ...
                    uncert_OCC_cons_array).^2;
                delta_part2_flow_vars = (extractdata(parital_flow_vars_dl) .* ...
                    uncert_flow_vars_array).^2;
               
                % Calculate uncertainty by summing all the components in
                % quadirature
                obj_uncert.HRR_kW(i) = (sum(delta_part2_cal_gas, "all")+ ...
                    sum(delta_part2_misc_CO, "all")+ ...
                    sum(delta_part2_misc_CO2, "all") + ...
                    sum(delta_part2_misc_O2, "all") + ...
                    sum(delta_part2_gas_cons, 'all') + ...
                    sum(delta_part2_OCC_cons, 'all') + ...
                    sum(delta_part2_flow_vars, "all") ).^0.5;

                 % Calcuate uncerty fractions to save for later
                ds_uncert_fraction.cal_gas(i,:) = delta_part2_cal_gas/...
                    obj_uncert.HRR_kW(i)^2;
                ds_uncert_fraction.misc_CO(i,:) = delta_part2_misc_CO./ ...
                    obj_uncert.HRR_kW(i)^2;
                ds_uncert_fraction.misc_CO2(i,:) = delta_part2_misc_CO2./ ...
                    obj_uncert.HRR_kW(i)^2;
                ds_uncert_fraction.misc_O2(i,:) = delta_part2_misc_O2./ ...
                    obj_uncert.HRR_kW(i)^2;
                ds_uncert_fraction.gas_cons(i,:) = delta_part2_gas_cons./ ...
                    obj_uncert.HRR_kW(i)^2;
                ds_uncert_fraction.OCC_cons(i,:) = delta_part2_OCC_cons./ ...
                    obj_uncert.HRR_kW(i)^2;
                ds_uncert_fraction.flow_vars(i,:) = delta_part2_flow_vars./ ...
                    obj_uncert.HRR_kW(i)^2;

                % Covert to expanned uncertaity (95% cofedence interval)
                obj_uncert.HRR_kW(i) = 2*obj_uncert.HRR_kW(i);

                end
            end

        end
    end


end
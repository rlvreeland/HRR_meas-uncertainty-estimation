classdef Fire_model
    properties
        % Define the molecular weights of the species of interest
        Mw_air_g_per_mol = 28.97;
        Mw_O2_g_per_mol = 32.00;
        Mw_CO_g_per_mol = 28.01;
        Mw_CO2_g_per_mol = 44.01;
        Mw_fuel_g_per_mol {mustBeNumeric}

        % Specify the gas constant and ambient conditions
        R_m3_Pa_per_K_mol = 8.314;
        Cp_kJ_per_kg_K = 1;
        P_amb_Pa = 101325;
        T_amb_K = 25+273.15;
        x_O2_amb = 0.2095;
        x_CO2_amb = 0.00044;
        x_CO_amb = 0;

        % Radiative fraction
        X_rad = 0.3;

        % Fuel heat Release and CO, CO2 yield
        delta_hc_kJ_per_kg {mustBeNumeric}
        yield_CO_g_per_g_fuel {mustBeNumeric}
        yield_CO2_g_per_g_fuel {mustBeNumeric}

        % Bidirrectional probe parameters
        f_Re = 1.07;

        % Calorimetery parameters
        alpha = 1.105; % dirrectly between pure C and H2
        E_kJ_per_kg = 13.1*10^3; % heat release per kg O2
        E_CO_kJ_per_kg = 17.6*10^3;

        % Paramagnetic O2 analyzer correctors
        CO_zero_offset_percent = 0.07; % need to apply mole fraction of the
        CO2_zero_offset_percent = -0.3;% species when conducted analysis

        % Calorimeter specifications
        D_duct_m {mustBeNumeric}
        h_hood_m {mustBeNumeric}
        flow_fraction {mustBeNumeric} %The fraction of the emtainment flow 
                                      % the hood is sucking

        % Calibration bottles
        x_CO_bottle_span {mustBeNumeric}
        x_CO2_bottle_span {mustBeNumeric}
        x_O2_bottle_span {mustBeNumeric}
        x_i_bottle_zero {mustBeNumeric}

        % Place holder for chacteristic fire values
        D_fire_m
        flame_height_m

        % Mass flow of fuel
        m_dot_f_kg_per_s

        % Coditions in the duct
        m_e_kg_per_s
        v_e_m_per_s
        Te_K
        Delta_P_Pa
        x_O2_measure
        x_CO2_measure
        x_CO_measure
        Re

      
    end

    methods
        % Use this method for a square gas burner
        function obj = Burner_diameter(obj, burner_length_m)
            obj.D_fire_m = (4*burner_length_m^2/pi)^0.5;
        end
        
        %Use this method to determine the diameter for a circular pool
        %fire, using correlation from Zabetakis and Burgess (1961) for pool
        %fires with diamters > 0.2
        function obj = Pool_fire_diameter(obj, HRR, m_dot_inf_kg_per_m2,...
                k_beta)
            %Preallocate array
            obj.D_fire_m = NaN(size(HRR));
            for i = 1:length(HRR)
                syms D
                m_dot_kg_per_m2 = m_dot_inf_kg_per_m2*(1-exp(-k_beta*D));
                eq = HRR(i) == m_dot_kg_per_m2*(D/2)^2*...
                    obj.delta_hc_kJ_per_kg*pi;
                Diameter_m = vpasolve(eq, D, .5);
                % Remove diamters that are outside of model (<= 0.2)
                if Diameter_m > 0.2
                    obj.D_fire_m(i) = Diameter_m;
                end
            end
        end

        % After determining the diamter use this function to obtain the
        % model HRR data (note the HRR has to be in kW)
        function obj = model_HRR_measurements(obj, HRR)
            % Specify number of HRR points
            length_HRR = length(HRR);
            % Check is Diamter of fire is an array with length HRR, if 
            % single value make into an array with all the same values
            if isscalar(obj.D_fire_m)
                obj.D_fire_m = ones(length_HRR,1)*obj.D_fire_m;
            end

            % Calculate the flame heights using hezkestad virtual orgin
            obj.flame_height_m = 0.235*HRR.^(2/5)-1.02*obj.D_fire_m;

            % Use Plume correlations to calulate entrainment of the plume
            z_o_m = obj.D_fire_m*(-1.02)+0.083*HRR.^(2/5);
            % Preallocate
            obj.m_e_kg_per_s = NaN(size(HRR));
            % calculate convective fraction
            X_conv = 1-obj.X_rad;
            for i = 1:length_HRR
                % remove case where flames reach the hood
                if obj.flame_height_m(i) < obj.h_hood_m
                    % Calculate the flow in the duct, some flow faction X 
                    % the entrainment rate
                    obj.m_e_kg_per_s(i) = obj.flow_fraction*(0.071*(HRR(i)...
                        *X_conv)^(1/3)*(obj.h_hood_m-z_o_m(i))^(5/3)*(1+...
                        0.027*(HRR(i)*X_conv)^(2/3)*(obj.h_hood_m-z_o_m(i))...
                        ^(-5/3)));
                end
            end

            % Calulate mass flow of fuel supplied
            obj.m_dot_f_kg_per_s = HRR/obj.delta_hc_kJ_per_kg;
            % Calculate the concentration of CO and CO2 in the duct give
            % the yield of each species
            m_CO2_fire_kg_per_s = obj.yield_CO2_g_per_g_fuel.*...
                obj.m_dot_f_kg_per_s;
            obj.x_CO2_measure = ((obj.m_e_kg_per_s-obj.m_dot_f_kg_per_s)*...
                obj.x_CO2_amb + m_CO2_fire_kg_per_s * obj.Mw_air_g_per_mol/...
                obj.Mw_CO2_g_per_mol)./obj.m_e_kg_per_s;
            % obj.x_CO2_measure = m_CO2_fire_kg_per_s./obj.m_e_kg_per_s*obj.Mw_air_g_per_mol/...
            %     obj.Mw_CO2_g_per_mol;
            m_CO_fire_kg_per_s = obj.yield_CO_g_per_g_fuel.*...
                obj.m_dot_f_kg_per_s;
            obj.x_CO_measure = ((obj.m_e_kg_per_s-obj.m_dot_f_kg_per_s)*...
                obj.x_CO_amb + m_CO_fire_kg_per_s * obj.Mw_air_g_per_mol/...
                obj.Mw_CO_g_per_mol)./obj.m_e_kg_per_s;
            % obj.x_CO_measure = m_CO_fire_kg_per_s./obj.m_e_kg_per_s*obj.Mw_air_g_per_mol/...
            %     obj.Mw_CO_g_per_mol;

            % preallocate arrays
            obj.Te_K = NaN(length_HRR,1);
            obj.x_O2_measure = NaN(length_HRR,1);
            % Solve for the concentrations of O2 and temperature in the
            % duct
            % Assuming no heat losses from the convective component of the
            % HRR
            for i = 1:length_HRR
                % pull out specfic test we are looking out for simplicity
                Q_HRR = HRR(i);
                m_dot_e_kg_per_s = obj.m_e_kg_per_s(i);
                x_CO2 = obj.x_CO2_measure(i);
                x_CO = obj.x_CO_measure(i);
                % Remove Nan data
                if isnan(m_dot_e_kg_per_s)
                else
                    % Solve for Te
                    obj.Te_K(i) = X_conv*Q_HRR/(obj.Cp_kJ_per_kg_K*...
                        m_dot_e_kg_per_s) + obj.T_amb_K;
                    syms x_O2
                    phi = (obj.x_O2_amb*(1-x_CO2-x_CO)-x_O2*...
                        (1-obj.x_CO2_amb))/((1-x_O2-x_CO2-x_CO)*obj.x_O2_amb);
                    eq1 = Q_HRR == (obj.E_kJ_per_kg*phi-(obj.E_CO_kJ_per_kg...
                        -obj.E_kJ_per_kg)*(1-phi)/2*x_CO/x_O2)*...
                        m_dot_e_kg_per_s/(1+phi*(obj.alpha-1))*...
                        obj.Mw_O2_g_per_mol/obj.Mw_air_g_per_mol*obj.x_O2_amb;
        % % O2 only
        % % % % % 
        % % % % % phi = (obj.x_O2_amb-x_O2)/((1-x_O2)*obj.x_O2_amb);
        % % % % % eq1 = Q_HRR == obj.E_kJ_per_kg*phi/(1+phi*(obj.alpha-1))*...
        % % % % %     obj.Mw_O2_g_per_mol/obj.Mw_air_g_per_mol*m_dot_e_kg_per_s*obj.x_O2_amb;
                    mole_frac_O2 = vpasolve(eq1, x_O2, .19);
                    % find solution with x_O2 closest to ambinent, this
                    % removes extraneous soultions
                    [~, index_sol] = max(mole_frac_O2);
                  
                    % This corrects measurement for the effect of CO and 
                    % CO2 concentration on output O2 fraction
                    obj.x_O2_measure(i) = mole_frac_O2(index_sol)...
                        +x_CO*obj.CO_zero_offset_percent/100 ...
                        +x_CO2*obj.CO2_zero_offset_percent/100;
                end
            end
            % Get the denisty of the gas in the duct
            rho_e_kg_per_m3 = obj.P_amb_Pa*obj.Mw_air_g_per_mol*10^(-3)./...
                (obj.R_m3_Pa_per_K_mol*obj.Te_K);
            % Calculate the velocity in the duct
            A_duct_m2 = pi/4*obj.D_duct_m.^2;
            obj.v_e_m_per_s = obj.m_e_kg_per_s./rho_e_kg_per_m3./A_duct_m2;

            % From the velocity calculate the differential pressure
            % measurements
            obj.Delta_P_Pa = obj.f_Re^2*obj.v_e_m_per_s.^2.*rho_e_kg_per_m3/2;

            % Calculate the Re to verify the flow regime
            % get viscoisty 
            % preallocate
            visco_Pa_s = NaN(length_HRR,1);
            for i = 1:length_HRR
                visco_Pa_s(i) = air_viscoisty(rho_e_kg_per_m3(i), obj.Te_K(i));
            end
            obj.Re = rho_e_kg_per_m3.*obj.v_e_m_per_s*2.54/100 ...
                ./visco_Pa_s;
        end
    end


end
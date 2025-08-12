function [inlet_Mach_array, inlet_temps, inlet_pressures, temp_ratio, pressure_ratio, beta] = inlet_working(theta_array, M_infty, gamma_infty, beta_array, T_infty, P_infty)
    % Function to calculate the flow properties downstream of a four oblique shock inlet.
    % Inputs:
    %   theta_array - array of deflection angles for the ramps
    %   M_infty - freestream Mach number
    %   gamma_infty - freestream heat capacity ratio
    %   beta_array - array of shock angles
    %   T_infty - freestream static temperature
    %   P_infty - freestream static pressure
    % Outputs:
    %   inlet_Mach_array - array of Mach numbers at each shock
    %   inlet_temps - array of temperatures at each shock
    %   inlet_pressures - array of pressures at each shock
    %   temp_ratio - temperature ratio across each shock
    %   pressure_ratio - pressure ratio across each shock
    %   beta - shock angles for each ramp

    % Initialize freestream conditions for the first shock
    M = M_infty;
    T = T_infty;
    P = P_infty;

    % Loop through each deflection angle to calculate flow properties
    for i = 1:length(theta_array)
        % Calculate possible shock angles using tan function
        numerator = (M^2 * (sind(beta_array).^2)) - 1;
        denominator = (M^2 * (gamma_infty + cosd(2 .* beta_array))) + 2;
        tantheta = 2 .* cotd(beta_array) .* (numerator ./ denominator);

        possible_theta = atand(tantheta);
        possible_theta(possible_theta < 0) = NaN;

        % Check if deflection angle exceeds maximum possible value
        if theta_array(i) > max(possible_theta)
            % Set flow properties to NaN and display error message
            T(i) = NaN;
            p(i) = NaN;
            sprintf('Theta, %.f, is greater than the maximum deflection angle at this Mach number, %.f', theta_array(i), M)
        else
            %% Calculating exact value of beta for weak shock solutions
            delta = 1; % Weak shock solutions
            term1 = (M^2 - 1)^2;
            term2 = 3 * (1 + ((gamma_infty - 1) / 2) * M^2);
            term3 = (1 + ((gamma_infty + 1) / 2) * M^2) * (tand(theta_array(i)))^2;
            lambda = sqrt(term1 - (term2 * term3));

            term4 = (M^2 - 1)^3;
            term5 = 9 * (1 + ((gamma_infty - 1) / 2) * M^2);
            term6 = (1 + ((gamma_infty - 1) / 2) * M^2 + ((gamma_infty + 1) / 4) * M^4) * tand(theta_array(i))^2;
            chi = (term4 - (term5 * term6)) / lambda^3;

            term7 = (4 * pi * delta + acos(chi)) / 3;
            numerator = M^2 - 1 + (2 * lambda * cos(term7));
            denominator = 3 * (1 + ((gamma_infty - 1) / 2) * M^2) * tand(theta_array(i));
            
            beta(i) = atand(numerator / denominator);

            %% Ratios and Flow Properties

            % Mach Number
            Mn = M * sind(beta(i));
            M2n = sqrt((1 + ((gamma_infty - 1) / 2) * Mn^2) / (gamma_infty * Mn^2 - ((gamma_infty - 1) / 2)));
            inlet_Mach_array(i) = M2n / sind(beta(i) - theta_array(i));
            M = inlet_Mach_array(i); % Update M for the next loop

            % Temperature Ratio
            part1 = 1 + ((2 * gamma_infty) / (gamma_infty + 1)) * (Mn^2 - 1);
            part2 = (2 + (gamma_infty - 1) * Mn^2) / ((gamma_infty + 1) * Mn^2);
            temp_ratio(i) = part1 * part2;
            inlet_temps(i) = temp_ratio(i) * T;
            T = inlet_temps(i); % Update T for the next loop

            % Pressure Ratio
            pressure_ratio(i) = 1 + ((2 * gamma_infty) / (gamma_infty + 1)) * (Mn^2 - 1);
            inlet_pressures(i) = pressure_ratio(i) * P;
            P = inlet_pressures(i); % Update P for the next loop
        end
    end
end

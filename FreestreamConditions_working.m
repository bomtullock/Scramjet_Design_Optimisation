function [T_infty, P_infty] = FreestreamConditions_working(h)
    % Function to calculate freestream static temperature (T_infty) and
    % static pressure (P_infty) based on flight altitude (h).
    % Inputs:
    %   h - flight altitude in meters
    % Outputs:
    %   T_infty - freestream static temperature in Kelvin
    %   P_infty - freestream static pressure in Pascals

    % Check if the altitude is within the valid range
    if (h < 0) || (h > 3e4)
        % Print error message and set outputs to NaN if altitude is out of range
        fprintf('The flight altitude is outside of the range\n')
        T_infty = NaN;
        P_infty = NaN;
    
    % Calculate T_infty and P_infty for altitude below 11,000 meters
    elseif h < 11000
        T_infty = 288.19 - (0.00649 * h);
        P_infty = 1.01e5 * (T_infty / 288.19)^5.256;

    % Calculate T_infty and P_infty for altitude between 11,000 and 25,000 meters
    elseif (11000 <= h) && (h < 25000)
        T_infty = 216.69;
        P_infty = 2.27e4 * exp(1)^(1.73 - (0.000157 * h));

    % Calculate T_infty and P_infty for altitude above 25,000 meters
    else
        T_infty = 141.94 + (0.00299 * h);
        P_infty = 2.49e3 * (T_infty / 216.6)^-11.388;
    end
end

function [Tt_infty, Pt_infty] = StagnationConditions_working(T_infty, P_infty, M_infty, gamma_infty)
    % Function to calculate freestream stagnation temperature (Tt_infty) and
    % stagnation pressure (Pt_infty) based on freestream conditions and Mach number.
    % Inputs:
    %   T_infty - freestream static temperature in Kelvin
    %   P_infty - freestream static pressure in Pascals
    %   M_infty - freestream Mach number
    %   gamma_infty - freestream heat capacity ratio
    % Outputs:
    %   Tt_infty - freestream stagnation temperature in Kelvin
    %   Pt_infty - freestream stagnation pressure in Pascals

    % Calculate freestream stagnation temperature (Tt_infty)
    Tt_infty = (1 + ((gamma_infty - 1) / 2 * M_infty^2)) * T_infty;

    % Calculate freestream stagnation pressure (Pt_infty)
    Pt_infty = (1 + ((gamma_infty - 1) / 2 * M_infty^2))^(gamma_infty / (gamma_infty - 1)) * P_infty;
end

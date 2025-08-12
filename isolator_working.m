function [M2, T2, P2, m2_dot] = isolator_working(iso_pt_ratio, iso_p_ratio, M1, T1, P1, gamma_isolator, iso_area, R)
    % Function to calculate the flow properties at the isolator exit
    % Inputs:
    %   iso_pt_ratio - Total pressure recovery ratio (p_t2 / p_t1)
    %   iso_p_ratio - Static pressure ratio (p2 / p1)
    %   M1 - Mach number at isolator entrance
    %   T1 - Static temperature at isolator entrance (K)
    %   P1 - Static pressure at isolator entrance (Pa)
    %   gamma_isolator - Heat capacity ratio at isolator
    %   iso_area - Cross-sectional area of isolator (m^2)
    %   R - Specific gas constant for air (J/(kgK))
    % Outputs:
    %   M2 - Mach number at isolator exit
    %   T2 - Static temperature at isolator exit (K)
    %   P2 - Static pressure at isolator exit (Pa)
    %   m2_dot - Mass flow rate at isolator exit (kg/s)

    % Calculate total pressure at isolator entrance
    Pt1 = P1 * (1 + (gamma_isolator - 1) / 2 * M1^2)^(gamma_isolator / (gamma_isolator - 1));

    % Calculate total pressure at isolator exit
    Pt2 = Pt1 * iso_pt_ratio;

    % Calculate static pressure at isolator exit
    P2 = P1 * iso_p_ratio;

    % Calculate Mach number at isolator exit
    M2 = sqrt((2 / (gamma_isolator - 1)) * ((Pt2 / P2)^((gamma_isolator - 1) / gamma_isolator) - 1));

    % Calculate static temperature at isolator exit
    T2 = T1 * (1 + (gamma_isolator - 1) / 2 * M1^2) / (1 + (gamma_isolator - 1) / 2 * M2^2);

    % Calculate density at isolator exit
    rho2 = P2 / (R * T2);

    % Calculate velocity at isolator exit
    V2 = M2 * sqrt(gamma_isolator * R * T2);

    % Calculate mass flow rate at isolator exit
    m2_dot = rho2 * iso_area * V2;
end

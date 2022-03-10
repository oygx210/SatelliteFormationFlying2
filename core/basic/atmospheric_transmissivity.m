function tau = atmospheric_transmissivity(satellite_elevation_angle)

% input:
% satellite_elevation_angle [rad]
% vector import is applicable

tau = 0.1283+0.7559*exp(-0.3878*sec(pi/2 - satellite_elevation_angle)); % atmosphere transmissivity


end
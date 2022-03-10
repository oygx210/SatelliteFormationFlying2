clear all;
consts = startup_formation_control();

% orbital elements
oe = zeros(6,1);
orbit_radius = 7266463;
orbit_altitude = orbit_radius - consts.rEarth;
gamma_beam_Sun = 32/60*pi/180 / 2; % half anglular size of the Sun, rad
S_spot_Sun = pi * (gamma_beam_Sun * orbit_altitude/1000)^2;

S_spot_req = [S_spot_Sun; 100*1e6;200*1e6; 500*1e6];
gamma_beam_req = sqrt(S_spot_req/pi) / orbit_altitude;

plot(S_spot_req / 1e6, rad2deg(gamma_beam_req), S_spot_Sun / 1e6, rad2deg(gamma_beam_Sun), 'or');
xlabel('Footprint area, km^2');
ylabel('\gamma_{beam}, deg');
grid on;
legend('\gamma_{beam} depending on required footprint area', '\gamma_{beam} for flat reflector');

% reflector selection
I0 = 1360;                  % The average intensity of solar energy at the earth distance (W/m^2)
Iref = (2.54e-6)/683;       % Intensity of the light from the brightness star in W/m^2 
rho = 0.92;                 % The reflectivity coefficient for the thin Mylar film coated with aluminum 
m_required = -3; % required magnitude for demonstration
I_req = Iref*10^(-m_required/2.5);

theta_sat_min = deg2rad(10);
gamma_refl = deg2rad((90 + 65)/2);
beta_max = 0.2668;
d_max = sqrt(orbit_radius^2 + consts.rEarth^2 - 2*orbit_radius*consts.rEarth*cos(beta_max));

tau_max = atmospheric_transmissivity(theta_sat_min);
I_sp_min = tau_max * cos(gamma_refl) * sin(theta_sat_min) / d_max^2;

Ar = I_req*4*tan(gamma_beam_req).^2./(I0*rho*I_sp_min);
Ar_Sun = I_req*4*tan(gamma_beam_Sun).^2./(I0*rho*I_sp_min);
Ar_100 = I_req*4*tan(gamma_beam_req).^2./(I0*rho*I_sp_min);
Ar_500 = I_req*4*tan(gamma_beam_req(end)).^2./(I0*rho*I_sp_min);

plot(S_spot_req / 1e6, Ar, '-o');
hold on;
plot(S_spot_Sun, Ar_Sun, 'or');
xlabel('Required footprint area, km^2');
ylabel('A_{r}, m^2');
grid on;


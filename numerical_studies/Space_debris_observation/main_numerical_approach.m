clear all;

%% Initialization
global environment 
environment = 'J2';

consts = startup_formation_control();

spacecraft.dr_observation = 200e3; 

orbit_epoch = datetime(2021, 9, 7, 0, 0, 0);
misson_duration = seconds(days(1));

% Target orbit
target_orbit = [consts.rEarth + 700e3; 0; NaN; NaN; 0; 0];
target_orbit(3) = get_SSO_inclination(target_orbit(1), target_orbit(2), consts);
target_orbit(4) = get_RAAN_for_terminator_orbit(orbit_epoch);
rv_target_orbit = oe2rv(target_orbit, consts);

% Studied orbit
studied_orbit = target_orbit;
studied_orbit(1) = target_orbit(1) + 100e3;
studied_orbit(3) = get_SSO_inclination(studied_orbit(1), studied_orbit(2), consts);

target_orbit_rv = oe2rv(target_orbit, consts);
studied_orbit_rv = oe2rv(studied_orbit, consts);

n_nodes = 100;

%% Show orbits

% plot_orbits([target_orbit_rv; studied_orbit_rv], consts);

%% orbit covarage
tic;
tau_max = orbit_coverage(target_orbit, studied_orbit, orbit_epoch, n_nodes, spacecraft, consts);
toc;
clear all;
consts = startup_formation_control();

%% Initialization

% Spacecraft parameters
spacecraft.DragArea = 0.0481;           % [m^2] cross-section area perpendicular to the incoming airflow
spacecraft.Cdrag = 2.2;                 % atmospheric drag coefficient
spacecraft.mass = 4;                    % [kg]

% Launch vehicle (LV) orbit [sma[m], ecc[-], inc[rad], RAAN[rad], AOP[rad], MA[rad]]
LV.coe(1) = 400e3 + consts.rEarth;      % [m] sma
LV.coe(2) = 4e-4;                       % [-] ecc
LV.coe(3) = 52*consts.deg2rad;          %[rad] inc
LV.coe(4) = 0;                          % [rad] RAAN
LV.coe(5) = 0;                          % [rad] AOP
LV.coe(6) = 0;                          % [rad] Mean anomaly
LV.rv = oe2rv(LV.coe, consts);

formation.N_sats = 2;
formation.ISD = 100e3; % ISD - intersatellite distance, [m]
formation.ISD_acceptable_error = 10e3; % [m]
formation.deployment_interval = 30; % [s]
formation.target_sma = LV.coe(1); % [m]
formation.max_sma_error = 100e3; % [m] - the value is that big because orbit maintenance is not ready yet
formation.commissioning_period = consts.day2sec*14; % [s]
formation.experiment_duration = consts.day2sec*10; % [s]
formation.orbital_configuration = [];

% Generating initial conditions for deployment
rv_orb_deployment = [0 0;
                     0 0;
                     0 0;
                     -1.02 -0.98;
                     0 0;
                     0 0];

%% Stage 1 - undocking & commissioning period

% Undocking - consecutive deployment of satellites with 30s delay;
rv_ECI_deployment(1:6,1) = orb2ECI(LV.rv, rv_orb_deployment(:,1), consts);
rv_ECI_deployment(7:12,1) = LV.rv;

disp('Undocking...');
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t_vec_deployment, rv_ECI_deployment] = ode45(@(t, rv) J2_atmo_Formation(t, rv, consts, spacecraft), 0:1:formation.deployment_interval, rv_ECI_deployment, options);
disp('Sat1 - released');
rv_ECI_deployment = rv_ECI_deployment';

rv_ECI_deployment(7:12,end) = orb2ECI(rv_ECI_deployment(7:12,end), rv_orb_deployment(:,2), consts);
disp('Sat2 - released');

% disp('Commissioning time has started');
% [t_vec_commissioning, rv_ECI_commissioning] = ode45(@(t, rv) J2_atmo_Formation(t, rv, consts, spacecraft), 0:formation.commissioning_period	, rv_ECI_deployment(:,end), options);
% disp('Commissioning time has been completed');
% 
% rv_ECI_commissioning = rv_ECI_commissioning';
% 
% % plot_coe_diff(t_vec_commissioning,  rv_ECI_commissioning(1:6,:), rv_ECI_commissioning(7:12,:), consts, 'd');
% % plot_ISD(t_vec_commissioning, rv_ECI_commissioning(1:6,:), rv_ECI_commissioning(7:12,:), consts, 'd');
% 
% %% Stage 2 - formation flying 
% sat1.coe = rv2coe(rv_ECI_commissioning(1:6,end), consts);
% sat2.coe = rv2coe(rv_ECI_commissioning(7:12,end), consts);
% delta_coe = sat1.coe - sat2.coe;
% 
% % defining leader-follower configuration
% if sat1.coe(8,end) < sat2.coe(8,end)
%    disp('Sat2 - leader, Sat1 - follower');
%    leader = rv_ECI_commissioning(7:12,end);
%    follower = rv_ECI_commissioning(1:6,end);
%    rv_ECI = [leader ; follower];
% else
%    rv_ECI = rv_ECI_commissioning(:,end);
%   disp('Sat1 - leader, Sat2 - follower');
% end

% These state vectors for leader and follower appeares when commissioning period is finished
% It is used to speed up debugging
rv_ECI = [-1356459.24692174;
        -5318393.28754868;
        -3936774.60174553;
        4919.74802489928;
        -4274.76364247790;
        4067.27013798862;
        -1449197.55959714;
        -5236248.86573836;
        -4012908.07969179;
        4885.55233097234;
        -4403.74493988956;
        3969.76127089436]; % 14 days commissioning -1.02, -0.98 deployment
 
formation.ISD = vecnorm(rv_ECI(1:3,end) - rv_ECI(7:9,end));

%% Formation flying

[t_vec_FF,rv_ECI_FF, maneuvers_dV, maneuvers_t] = train_formati on_flying_ISD_keeping(rv_ECI, formation.experiment_duration, consts, spacecraft, formation);

% rv_ECI = [rv_ECI_commissioning, rv_ECI_FF];
% t_vec = [t_vec_commissioning; t_vec_FF + t_vec_commissioning(end)];

plot_ISD(t_vec_FF, rv_ECI_FF(1:6,:),  rv_ECI_FF(7:12,:), consts, formation, 'd');

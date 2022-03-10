% The script is used to calculate the semi-major orbit for a Sun-synchronous repeating
% ground track orbit

clear all;
consts = startup_formation_control();

k_rev2rep = 15;
kday2rep = 1;
h_min = 100e3; 
h_max = 1000e3;
sma = get_SSO_RGT_orbit_sma(k_rev2rep, kday2rep, consts);
h = sma - consts.rEarth;
incl = get_SSO_inclination(sma(1), 0, consts);

function sma = get_SSO_RGT_orbit_sma(k_rev2rep, k_day2rep, consts)
    % circular orbits are considered
    sma = [];
    T_greenwich_nodal = 2*pi / (consts.omegaEarth - consts.EarthMeanMotion);
    k_revpday = k_rev2rep / k_day2rep;
    n_assumption = k_revpday * consts.omegaEarth;
    a_assumption = (consts.muEarth*(1/n_assumption)^2)^(1/3);
    
    delta_a = [-50e3:0.1:50e3];
    
    for i = 1:length(delta_a)    
        eps = 1e-3;
    
        a_local = a_assumption + delta_a(i);
        incl_SSO = acos(-2*a_local^(7/2)*consts.EarthMeanMotion/(3*consts.rEarth_equatorial^2*consts.J2*sqrt(consts.muEarth)));
        T_kep = 2*pi*sqrt(a_local^3/consts.muEarth);
        T_sat_nodal_local = T_kep*(1 - 3/2*consts.J2*(consts.rEarth_equatorial/a_local)^2*(3 - 4*sin(incl_SSO)^2));
        delta(i) = abs((T_sat_nodal_local*k_rev2rep - T_greenwich_nodal*k_day2rep));
        if delta(i) < eps
            sma = [sma; a_local];
        end
    
    end
end


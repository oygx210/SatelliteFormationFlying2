function formation_fuel_level_out = write_off_fuel(formation_fuel_level_history, maneuvering_sats, sats_dV, spacecraft, consts)

formation_fuel_level_current = formation_fuel_level_history(:,end);
formation_fuel_level_updated = formation_fuel_level_current;

for i = 1:length(maneuvering_sats)    
    spacecraft_wet_mass_updated = (spacecraft.dry_mass + formation_fuel_level_current(maneuvering_sats(i)))*exp(-sats_dV(i)/spacecraft.thruster_Isp/consts.g);
    fuel_spent = (spacecraft.dry_mass + formation_fuel_level_current(maneuvering_sats(i))) - spacecraft_wet_mass_updated;
    formation_fuel_level_updated(maneuvering_sats(i)) = formation_fuel_level_current(maneuvering_sats(i)) - fuel_spent;
end

formation_fuel_level_out = [formation_fuel_level_history, formation_fuel_level_updated];
end


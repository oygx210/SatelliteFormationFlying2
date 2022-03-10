
cost_m_post = cost_matrix_dV.continuous_post_correction_maneuvers_dV;
cost_m_imp = cost_matrix_dV.impulsive_maneuvers_dV;

for i = 1:formation.N_active_sats
    predicted_imp_reconf(i,1) = cost_m_imp(match_matrix(i,1,2), match_matrix(i,2,2));
    predicted_post(i,1) = cost_m_post(match_matrix(i,1,2), match_matrix(i,2,2));
end

simulated_post = post_coorection_maneuvers(:,2);
simulated_imp = maneuvers(:,5) - simulated_post;

mean_error_imp = mean(simulated_imp - predicted_imp_reconf);
mean_error_post = mean(simulated_post - predicted_post);
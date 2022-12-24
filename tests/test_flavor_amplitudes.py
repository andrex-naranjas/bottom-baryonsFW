from bottomfw.common.flavor_amplitudes import FlavorAmplitudes

test_flavor_amplitudes = FlavorAmplitudes(baryons="omegas")

test_transition_lambda_0_sigma_0_1 = test_flavor_amplitudes.flavor_matrix_elements(n=1, x_i = test_flavor_amplitudes.Lambda_b_0(), y_j = test_flavor_amplitudes.Sigma_b_0())
print(test_transition_lambda_0_sigma_0_1)

test_transition_lambda_0_sigma_0_2 = test_flavor_amplitudes.flavor_matrix_elements(n=2, x_i = test_flavor_amplitudes.Lambda_b_0(), y_j = test_flavor_amplitudes.Sigma_b_0())
print(test_transition_lambda_0_sigma_0_2)

test_transition_lambda_0_sigma_0_3 = test_flavor_amplitudes.flavor_matrix_elements(n=3, x_i = test_flavor_amplitudes.Lambda_b_0(), y_j = test_flavor_amplitudes.Sigma_b_0())
print(test_transition_lambda_0_sigma_0_3)

test_transition_Xi_0_Xi_prime_0_1 = test_flavor_amplitudes.flavor_matrix_elements(n=1, x_i = test_flavor_amplitudes.Xi_b_0(), y_j = test_flavor_amplitudes.Xi_b_pr_0())
print(test_transition_Xi_0_Xi_prime_0_1)

test_transition_Xi_0_Xi_prime_0_2 = test_flavor_amplitudes.flavor_matrix_elements(n=2, x_i = test_flavor_amplitudes.Xi_b_0(), y_j = test_flavor_amplitudes.Xi_b_pr_0())
print(test_transition_Xi_0_Xi_prime_0_2)

test_transition_Xi_0_Xi_prime_0_3 = test_flavor_amplitudes.flavor_matrix_elements(n=3, x_i = test_flavor_amplitudes.Xi_b_0(), y_j = test_flavor_amplitudes.Xi_b_pr_0())
print(test_transition_Xi_0_Xi_prime_0_3)

test_transition_Xi_m_Xi_prime_m_1 = test_flavor_amplitudes.flavor_matrix_elements(n = 1, x_i = test_flavor_amplitudes.Xi_b_m(), y_j = test_flavor_amplitudes.Xi_b_pr_m()) 
print(test_transition_Xi_m_Xi_prime_m_1)

test_transition_Xi_m_Xi_prime_m_2 = test_flavor_amplitudes.flavor_matrix_elements(n=2, x_i = test_flavor_amplitudes.Xi_b_m(), y_j = test_flavor_amplitudes.Xi_b_pr_m()) 
print(test_transition_Xi_m_Xi_prime_m_2)

test_transition_Xi_m_Xi_prime_m_3 = test_flavor_amplitudes.flavor_matrix_elements(n = 3, x_i = test_flavor_amplitudes.Xi_b_m(), y_j = test_flavor_amplitudes.Xi_b_pr_m())
print(test_transition_Xi_m_Xi_prime_m_3)

test_transition_Lambda_0_Lambda_0_1 = test_flavor_amplitudes.flavor_matrix_elements(n = 1, x_i = test_flavor_amplitudes.Lambda_b_0(), y_j = test_flavor_amplitudes.Lambda_b_0())
print(test_transition_Lambda_0_Lambda_0_1)

test_transition_Lambda_0_Lambda_0_2 = test_flavor_amplitudes.flavor_matrix_elements(n = 2, x_i = test_flavor_amplitudes.Lambda_b_0(), y_j = test_flavor_amplitudes.Lambda_b_0())
print(test_transition_Lambda_0_Lambda_0_2)

test_transition_Lambda_0_Lambda_0_3 = test_flavor_amplitudes.flavor_matrix_elements(n = 3, x_i = test_flavor_amplitudes.Lambda_b_0(), y_j = test_flavor_amplitudes.Lambda_b_0())
print(test_transition_Lambda_0_Lambda_0_3)

test_transition_Xi_b_0_Xi_b_0_1 = test_flavor_amplitudes.flavor_matrix_elements(n = 1, x_i = test_flavor_amplitudes.Xi_b_0(), y_j = test_flavor_amplitudes.Xi_b_0())
print(test_transition_Xi_b_0_Xi_b_0_1)

test_transition_Xi_b_0_Xi_b_0_2 = test_flavor_amplitudes.flavor_matrix_elements(n = 2, x_i = test_flavor_amplitudes.Xi_b_0(), y_j = test_flavor_amplitudes.Xi_b_0())
print(test_transition_Xi_b_0_Xi_b_0_2)

test_transition_Xi_b_0_Xi_b_0_3 = test_flavor_amplitudes.flavor_matrix_elements(n = 3, x_i = test_flavor_amplitudes.Xi_b_0(), y_j = test_flavor_amplitudes.Xi_b_0())
print(test_transition_Xi_b_0_Xi_b_0_3)

test_transition_Xi_b_m_Xi_b_m_1 = test_flavor_amplitudes.flavor_matrix_elements(n = 1, x_i = test_flavor_amplitudes.Xi_b_m(), y_j = test_flavor_amplitudes.Xi_b_m())
print(test_transition_Xi_b_m_Xi_b_m_1)

test_transition_Xi_b_m_Xi_b_m_2 = test_flavor_amplitudes.flavor_matrix_elements(n = 2, x_i = test_flavor_amplitudes.Xi_b_m(), y_j = test_flavor_amplitudes.Xi_b_m())
print(test_transition_Xi_b_m_Xi_b_m_2)

test_transition_Xi_b_m_Xi_b_m_3 = test_flavor_amplitudes.flavor_matrix_elements(n = 3, x_i = test_flavor_amplitudes.Xi_b_m(), y_j = test_flavor_amplitudes.Xi_b_m())
print(test_transition_Xi_b_m_Xi_b_m_3)

test_transition_Sigma_b_p_Sigma_b_p_1 = test_flavor_amplitudes.flavor_matrix_elements(n = 1, x_i = test_flavor_amplitudes.Sigma_b_p(), y_j = test_flavor_amplitudes.Sigma_b_p())
print(test_transition_Sigma_b_p_Sigma_b_p_1)

test_transition_Sigma_b_p_Sigma_b_p_2 = test_flavor_amplitudes.flavor_matrix_elements(n = 2, x_i = test_flavor_amplitudes.Sigma_b_p(), y_j = test_flavor_amplitudes.Sigma_b_p())
print(test_transition_Sigma_b_p_Sigma_b_p_2)

test_transition_Sigma_b_p_Sigma_b_p_3 = test_flavor_amplitudes.flavor_matrix_elements(n = 3, x_i = test_flavor_amplitudes.Sigma_b_p(), y_j = test_flavor_amplitudes.Sigma_b_p())
print(test_transition_Sigma_b_p_Sigma_b_p_3)

test_transition_Sigma_b_0_Sigma_b_0_1 = test_flavor_amplitudes.flavor_matrix_elements(n = 1, x_i = test_flavor_amplitudes.Sigma_b_0(), y_j = test_flavor_amplitudes.Sigma_b_0())
print(test_transition_Sigma_b_0_Sigma_b_0_1)

test_transition_Sigma_b_0_Sigma_b_0_2 = test_flavor_amplitudes.flavor_matrix_elements(n = 2, x_i = test_flavor_amplitudes.Sigma_b_0(), y_j = test_flavor_amplitudes.Sigma_b_0())
print(test_transition_Sigma_b_0_Sigma_b_0_2)

test_transition_Sigma_b_0_Sigma_b_0_3 = test_flavor_amplitudes.flavor_matrix_elements(n = 3, x_i = test_flavor_amplitudes.Sigma_b_0(), y_j = test_flavor_amplitudes.Sigma_b_0())
print(test_transition_Sigma_b_0_Sigma_b_0_3)

test_transition_Sigma_b_m_Sigma_b_m_1 = test_flavor_amplitudes.flavor_matrix_elements(n = 1, x_i = test_flavor_amplitudes.Sigma_b_m(), y_j = test_flavor_amplitudes.Sigma_b_m())
print(test_transition_Sigma_b_m_Sigma_b_m_1)

test_transition_Sigma_b_m_Sigma_b_m_2 = test_flavor_amplitudes.flavor_matrix_elements(n = 2, x_i = test_flavor_amplitudes.Sigma_b_m(), y_j = test_flavor_amplitudes.Sigma_b_m())
print(test_transition_Sigma_b_m_Sigma_b_m_2)

test_transition_Sigma_b_m_Sigma_b_m_3 = test_flavor_amplitudes.flavor_matrix_elements(n = 3, x_i = test_flavor_amplitudes.Sigma_b_m(), y_j = test_flavor_amplitudes.Sigma_b_m())
print(test_transition_Sigma_b_m_Sigma_b_m_3)

test_transition_Xi_b_pr_0_Xi_b_pr_0_1 = test_flavor_amplitudes.flavor_matrix_elements(n = 1, x_i = test_flavor_amplitudes.Xi_b_pr_0(), y_j = test_flavor_amplitudes.Xi_b_pr_0())
print(test_transition_Xi_b_pr_0_Xi_b_pr_0_1)

test_transition_Xi_b_pr_0_Xi_b_pr_0_2 = test_flavor_amplitudes.flavor_matrix_elements(n = 2, x_i = test_flavor_amplitudes.Xi_b_pr_0(), y_j = test_flavor_amplitudes.Xi_b_pr_0())
print(test_transition_Xi_b_pr_0_Xi_b_pr_0_2)

test_transition_Xi_b_pr_0_Xi_b_pr_0_3 = test_flavor_amplitudes.flavor_matrix_elements(n = 3, x_i = test_flavor_amplitudes.Xi_b_pr_0(), y_j = test_flavor_amplitudes.Xi_b_pr_0())
print(test_transition_Xi_b_pr_0_Xi_b_pr_0_3)

test_transition_Xi_b_pr_m_Xi_b_pr_m_1 = test_flavor_amplitudes.flavor_matrix_elements(n = 1, x_i = test_flavor_amplitudes.Xi_b_pr_m(), y_j = test_flavor_amplitudes.Xi_b_pr_m())
print(test_transition_Xi_b_pr_m_Xi_b_pr_m_1)

test_transition_Xi_b_pr_m_Xi_b_pr_m_2 = test_flavor_amplitudes.flavor_matrix_elements(n = 2, x_i = test_flavor_amplitudes.Xi_b_pr_m(), y_j = test_flavor_amplitudes.Xi_b_pr_m())
print(test_transition_Xi_b_pr_m_Xi_b_pr_m_2)

test_transition_Xi_b_pr_m_Xi_b_pr_m_3 = test_flavor_amplitudes.flavor_matrix_elements(n = 3, x_i = test_flavor_amplitudes.Xi_b_pr_m(), y_j = test_flavor_amplitudes.Xi_b_pr_m())
print(test_transition_Xi_b_pr_m_Xi_b_pr_m_3)

test_transition_Omega_b_m_Omega_b_m_1 = test_flavor_amplitudes.flavor_matrix_elements(n = 1, x_i = test_flavor_amplitudes.Omega_b_m(), y_j = test_flavor_amplitudes.Omega_b_m())
print(test_transition_Omega_b_m_Omega_b_m_1)

test_transition_Omega_b_m_Omega_b_m_2 = test_flavor_amplitudes.flavor_matrix_elements(n = 2, x_i = test_flavor_amplitudes.Omega_b_m(), y_j = test_flavor_amplitudes.Omega_b_m())
print(test_transition_Omega_b_m_Omega_b_m_2)

test_transition_Omega_b_m_Omega_b_m_3 = test_flavor_amplitudes.flavor_matrix_elements(n = 3, x_i = test_flavor_amplitudes.Omega_b_m(), y_j = test_flavor_amplitudes.Omega_b_m())
print(test_transition_Omega_b_m_Omega_b_m_3)

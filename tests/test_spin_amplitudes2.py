from bottomfw.common.spin_amplitudes import SpinAmplitudes

test_amplitudes = SpinAmplitudes(baryons="omegas")

test_amp_slo1_s1_2_s3_2 = test_amplitudes.matrix_elements(n=1, x_i=test_amplitudes.symmetric_states(m_proj=1/2), y_j=test_amplitudes.symmetric_states(m_proj=3/2))
print(test_amp_slo1_s1_2_s3_2)

test_amp_slo1_s1_2m_s1_2 = test_amplitudes.matrix_elements(n=1, x_i=test_amplitudes.symmetric_states(m_proj=-1/2), y_j=test_amplitudes.symmetric_states(m_proj=1/2))
print(test_amp_slo1_s1_2m_s1_2)

test_amp_slo1_s3_2m_s1_2m = test_amplitudes.matrix_elements(n=1, x_i=test_amplitudes.symmetric_states(m_proj=-3/2), y_j=test_amplitudes.symmetric_states(m_proj=-1/2))
print(test_amp_slo1_s3_2m_s1_2m)

test_amp_slo1_s1_2m_r1_2 = test_amplitudes.matrix_elements(n=1, x_i=test_amplitudes.symmetric_states(m_proj=-1/2), y_j=test_amplitudes.rho_states(m_proj=1/2))
print(test_amp_slo1_s1_2m_r1_2)

test_amp_slo1_s3_2m_r1_2m = test_amplitudes.matrix_elements(n=1, x_i=test_amplitudes.symmetric_states(m_proj=-3/2), y_j=test_amplitudes.rho_states(m_proj=-1/2))
print(test_amp_slo1_s3_2m_r1_2m)

test_amp_slo1_s1_2m_l1_2 = test_amplitudes.matrix_elements(n=1, x_i=test_amplitudes.symmetric_states(m_proj=-1/2), y_j=test_amplitudes.lambda_states(m_proj=1/2))
print(test_amp_slo1_s1_2m_l1_2)

test_amp_slo1_s3_2m_l1_2m = test_amplitudes.matrix_elements(n=1, x_i=test_amplitudes.symmetric_states(m_proj=-3/2), y_j=test_amplitudes.lambda_states(m_proj=-1/2))
print(test_amp_slo1_s3_2m_l1_2m)

test_amp_slo1_r1_2_s3_2 = test_amplitudes.matrix_elements(n=1, x_i=test_amplitudes.rho_states(m_proj=1/2), y_j=test_amplitudes.symmetric_states(m_proj=3/2))
print(test_amp_slo1_r1_2_s3_2)

test_amp_slo1_r1_2m_s1_2 = test_amplitudes.matrix_elements(n=1, x_i=test_amplitudes.rho_states(m_proj=-1/2), y_j=test_amplitudes.symmetric_states(m_proj=1/2))
print(test_amp_slo1_r1_2m_s1_2)

test_amp_slo1_r1_2m_l1_2 = test_amplitudes.matrix_elements(n=1, x_i=test_amplitudes.rho_states(m_proj=-1/2), y_j=test_amplitudes.lambda_states(m_proj=1/2))
print(test_amp_slo1_r1_2m_l1_2)

test_amp_slo1_l1_2_s3_2 = test_amplitudes.matrix_elements(n=1, x_i=test_amplitudes.lambda_states(m_proj=1/2), y_j=test_amplitudes.symmetric_states(m_proj=3/2))
print(test_amp_slo1_l1_2_s3_2)

test_amp_slo1_l1_2m_s1_2 = test_amplitudes.matrix_elements(n=1, x_i=test_amplitudes.lambda_states(m_proj=-1/2), y_j=test_amplitudes.symmetric_states(m_proj=1/2))
print(test_amp_slo1_l1_2m_s1_2)

test_amp_slo1_l1_2m_r1_2 = test_amplitudes.matrix_elements(n=1, x_i=test_amplitudes.lambda_states(m_proj=-1/2), y_j=test_amplitudes.rho_states(m_proj=1/2))
print(test_amp_slo1_l1_2m_r1_2)

test_amp_slo1_l1_2m_l1_2 = test_amplitudes.matrix_elements(n=1, x_i=test_amplitudes.lambda_states(m_proj=-1/2), y_j=test_amplitudes.lambda_states(m_proj=1/2))
print(test_amp_slo1_l1_2m_l1_2)


test_amp_slo2_s1_2_s3_2 = test_amplitudes.matrix_elements(n=2, x_i=test_amplitudes.symmetric_states(m_proj=1/2), y_j=test_amplitudes.symmetric_states(m_proj=3/2))
print(test_amp_slo2_s1_2_s3_2)

test_amp_slo2_s1_2m_s1_2 = test_amplitudes.matrix_elements(n=2, x_i=test_amplitudes.symmetric_states(m_proj=-1/2), y_j=test_amplitudes.symmetric_states(m_proj=1/2))
print(test_amp_slo2_s1_2m_s1_2)

test_amp_slo2_s3_2m_s1_2m = test_amplitudes.matrix_elements(n=2, x_i=test_amplitudes.symmetric_states(m_proj=-3/2), y_j=test_amplitudes.symmetric_states(m_proj=-1/2))
print(test_amp_slo2_s3_2m_s1_2m)

test_amp_slo2_s1_2m_r1_2 = test_amplitudes.matrix_elements(n=2, x_i=test_amplitudes.symmetric_states(m_proj=-1/2), y_j=test_amplitudes.rho_states(m_proj=1/2))
print(test_amp_slo2_s1_2m_r1_2)

test_amp_slo2_s3_2m_r1_2m = test_amplitudes.matrix_elements(n=2, x_i=test_amplitudes.symmetric_states(m_proj=-3/2), y_j=test_amplitudes.rho_states(m_proj=-1/2))
print(test_amp_slo2_s3_2m_r1_2m)

test_amp_slo2_s1_2m_l1_2 = test_amplitudes.matrix_elements(n=2, x_i=test_amplitudes.symmetric_states(m_proj=-1/2), y_j=test_amplitudes.lambda_states(m_proj=1/2))
print(test_amp_slo2_s1_2m_l1_2)

test_amp_slo2_s3_2m_l1_2m = test_amplitudes.matrix_elements(n=2, x_i=test_amplitudes.symmetric_states(m_proj=-3/2), y_j=test_amplitudes.lambda_states(m_proj=-1/2))
print(test_amp_slo2_s3_2m_l1_2m)

test_amp_slo2_r1_2_s3_2 = test_amplitudes.matrix_elements(n=2, x_i=test_amplitudes.rho_states(m_proj=1/2), y_j=test_amplitudes.symmetric_states(m_proj=3/2))
print(test_amp_slo2_r1_2_s3_2)

test_amp_slo2_r1_2m_s1_2 = test_amplitudes.matrix_elements(n=2, x_i=test_amplitudes.rho_states(m_proj=-1/2), y_j=test_amplitudes.symmetric_states(m_proj=1/2))
print(test_amp_slo2_r1_2m_s1_2)

test_amp_slo2_r1_2m_l1_2 = test_amplitudes.matrix_elements(n=2, x_i=test_amplitudes.rho_states(m_proj=-1/2), y_j=test_amplitudes.lambda_states(m_proj=1/2))
print(test_amp_slo2_r1_2m_l1_2)

test_amp_slo2_l1_2_s3_2 = test_amplitudes.matrix_elements(n=2, x_i=test_amplitudes.lambda_states(m_proj=1/2), y_j=test_amplitudes.symmetric_states(m_proj=3/2))
print(test_amp_slo2_l1_2_s3_2)

test_amp_slo2_l1_2m_s1_2 = test_amplitudes.matrix_elements(n=2, x_i=test_amplitudes.lambda_states(m_proj=-1/2), y_j=test_amplitudes.symmetric_states(m_proj=1/2))
print(test_amp_slo2_l1_2m_s1_2)

test_amp_slo2_l1_2m_r1_2 = test_amplitudes.matrix_elements(n=2, x_i=test_amplitudes.lambda_states(m_proj=-1/2), y_j=test_amplitudes.rho_states(m_proj=1/2))
print(test_amp_slo2_l1_2m_r1_2)

test_amp_slo2_l1_2m_l1_2 = test_amplitudes.matrix_elements(n=2, x_i=test_amplitudes.lambda_states(m_proj=-1/2), y_j=test_amplitudes.lambda_states(m_proj=1/2))
print(test_amp_slo2_l1_2m_l1_2)


test_amp_slo3_s1_2_s3_2 = test_amplitudes.matrix_elements(n=3, x_i=test_amplitudes.symmetric_states(m_proj=1/2), y_j=test_amplitudes.symmetric_states(m_proj=3/2))
print(test_amp_slo3_s1_2_s3_2)

test_amp_slo3_s1_2m_s1_2 = test_amplitudes.matrix_elements(n=3, x_i=test_amplitudes.symmetric_states(m_proj=-1/2), y_j=test_amplitudes.symmetric_states(m_proj=1/2))
print(test_amp_slo3_s1_2m_s1_2)

test_amp_slo3_s3_2m_s1_2m = test_amplitudes.matrix_elements(n=3, x_i=test_amplitudes.symmetric_states(m_proj=-3/2), y_j=test_amplitudes.symmetric_states(m_proj=-1/2))
print(test_amp_slo3_s3_2m_s1_2m)

test_amp_slo3_s1_2m_l1_2 = test_amplitudes.matrix_elements(n=3, x_i=test_amplitudes.symmetric_states(m_proj=-1/2), y_j=test_amplitudes.lambda_states(m_proj=1/2))
print(test_amp_slo3_s1_2m_l1_2)

test_amp_slo3_s3_2m_l1_2m = test_amplitudes.matrix_elements(n=3, x_i=test_amplitudes.symmetric_states(m_proj=-3/2), y_j=test_amplitudes.lambda_states(m_proj=-1/2))
print(test_amp_slo3_s3_2m_l1_2m)

test_amp_slo3_r1_2m_r1_2 = test_amplitudes.matrix_elements(n=3, x_i=test_amplitudes.rho_states(m_proj=-1/2), y_j=test_amplitudes.rho_states(m_proj=1/2))
print(test_amp_slo3_r1_2m_r1_2)

test_amp_slo3_l1_2_s3_2 = test_amplitudes.matrix_elements(n=3, x_i=test_amplitudes.lambda_states(m_proj=1/2), y_j=test_amplitudes.symmetric_states(m_proj=3/2))
print(test_amp_slo3_l1_2_s3_2)

test_amp_slo3_l1_2m_s1_2 = test_amplitudes.matrix_elements(n=3, x_i=test_amplitudes.lambda_states(m_proj=-1/2), y_j=test_amplitudes.symmetric_states(m_proj=1/2))
print(test_amp_slo3_l1_2m_s1_2)

test_amp_slo3_l1_2m_l1_2 = test_amplitudes.matrix_elements(n=3, x_i=test_amplitudes.lambda_states(m_proj=-1/2), y_j=test_amplitudes.lambda_states(m_proj=1/2))
print(test_amp_slo3_l1_2m_l1_2)

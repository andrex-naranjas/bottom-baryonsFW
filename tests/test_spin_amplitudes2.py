from bottomfw.common.spin_amplitudes import SpinAmplitudes

test_amplitudes = SpinAmplitudes(baryons="omegas")

test_su2_1 = test_amplitudes.su2_generator(index=1)
print(test_su2_1)

test_su2_2 = test_amplitudes.su2_generator(index=2)
print(test_su2_2)

test_su2_3 = test_amplitudes.su2_generator(index=3)
print(test_su2_3)

test_identity = test_amplitudes.identity_matrix(dim=2)
print(test_identity)

test_ladder_operatorp = test_amplitudes.ladder_operator(sign=1, dim1=1, dim2=2)
print(test_ladder_operatorp)

test_ladder_operatorm = test_amplitudes.ladder_operator(sign=-1, dim1=1, dim2=2)
print(test_ladder_operatorm)

test_tensor_product2 = test_amplitudes.tensor_product(index=2, dim=2)
print(test_tensor_product2)

test_ladder_operator_tensorp = test_amplitudes.ladder_operator_tensor(sign=1, dim1=1, dim2=2)
print(test_ladder_operator_tensorp)

test_symmetric_state_3_2 = test_amplitudes.symmetric_states(m_proj=3/2)
print(test_symmetric_state_3_2)

test_symmetric_state_1_2 = test_amplitudes.symmetric_states(m_proj=1/2)
print(test_symmetric_state_1_2)

test_symmetric_state_1_2m = test_amplitudes.symmetric_states(m_proj=-1/2)
print(test_symmetric_state_1_2m)

test_symmetric_state_3_2m = test_amplitudes.symmetric_states(m_proj=-3/2)
print(test_symmetric_state_3_2m)

test_rho_state_1_2 = test_amplitudes.rho_states(m_proj=1/2)
print(test_rho_state_1_2)

test_rho_state_1_2m = test_amplitudes.rho_states(m_proj=-1/2)
print(test_rho_state_1_2m)

test_lambda_state_1_2 = test_amplitudes.lambda_states(m_proj=1/2)
print(test_lambda_state_1_2)

test_lambda_state_1_2m = test_amplitudes.lambda_states(m_proj=-1/2)
print(test_lambda_state_1_2m)

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

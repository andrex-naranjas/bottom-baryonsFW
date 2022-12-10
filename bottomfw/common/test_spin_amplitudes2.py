from spin_amplitudes import SpinAmplitudes

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


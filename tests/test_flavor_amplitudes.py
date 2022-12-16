from bottomfw.common.flavor_amplitudes import FlavorAmplitudes

test_flavor_amplitudes = FlavorAmplitudes(baryons="omegas")

test_magnetic_operators = test_flavor_amplitudes.magnetic_operators(index=1)
print(test_magnetic_operators)

test_magnetic_operators = test_flavor_amplitudes.magnetic_operators(index=1)
print(test_magnetic_operators)

test_tensor_product_operator_1 = test_flavor_amplitudes.tensor_product_operator_1
print(test_tensor_product_operator_1())

test_tensor_product_operator_2 = test_flavor_amplitudes.tensor_product_operator_2
print(test_tensor_product_operator_2())

test_tensor_product_operator_3 = test_flavor_amplitudes.tensor_product_operator_3
print(test_tensor_product_operator_3())

test_Sigma_0 = test_flavor_amplitudes.Sigma_0
print(test_Sigma_0())

test_Lambda_0 = test_flavor_amplitudes.Lambda_0
print(test_Lambda_0())

test_Xi_prime_0 = test_flavor_amplitudes.Xi_prime_0
print(test_Xi_prime_0())

test_Xi_0 = test_flavor_amplitudes.Xi_0
print(test_Xi_0())

test_Xi_prime_m = test_flavor_amplitudes.Xi_prime_m
print(test_Xi_prime_m())

test_Xi_m = test_flavor_amplitudes.Xi_m
print(test_Xi_m())




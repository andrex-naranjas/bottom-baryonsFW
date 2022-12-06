from bottomfw.common.spin_amplitudes import SpinAmplitudes

test_amplitudes = SpinAmplitudes(baryons="omegas")

test_su2_generator = test_amplitudes.su2_generator(index=3)
print(test_su2_generator)

test_identity = test_amplitudes.identity_matrix(dim=3)
print(test_identity)

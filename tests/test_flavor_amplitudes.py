import unittest
import pytest


class TestFlavorAmplitudes(unittest.TestCase):
    """
    Tests methods in flavor_amplitudes
    """
    @pytest.fixture(autouse=True)
    def initdir(self, tmpdir):
        print(tmpdir)
        tmpdir.chdir()

    def test_flavor_matrix_elements(self):
        """
        Test the matrix element method for flavor
        """
        from bottomfw.common.flavor_amplitudes import FlavorAmplitudes
        import numpy as np
        import sympy as sp
        test_flavor_amplitudes = FlavorAmplitudes(baryons="omegas")
        test_transition_lambda_0_sigma_0_1 = test_flavor_amplitudes.flavor_matrix_elements(n=1, x_i = test_flavor_amplitudes.Lambda_b_0(), y_j = test_flavor_amplitudes.Sigma_b_0())
        test_transition_lambda_0_sigma_0_2 = test_flavor_amplitudes.flavor_matrix_elements(n=2, x_i = test_flavor_amplitudes.Lambda_b_0(), y_j = test_flavor_amplitudes.Sigma_b_0())
        test_transition_lambda_0_sigma_0_3 = test_flavor_amplitudes.flavor_matrix_elements(n=3, x_i = test_flavor_amplitudes.Lambda_b_0(), y_j = test_flavor_amplitudes.Sigma_b_0())
        test_transition_Xi_0_Xi_prime_0_1 = test_flavor_amplitudes.flavor_matrix_elements(n=1, x_i = test_flavor_amplitudes.Xi_b_0(), y_j = test_flavor_amplitudes.Xi_b_pr_0())
        test_transition_Xi_0_Xi_prime_0_2 = test_flavor_amplitudes.flavor_matrix_elements(n=2, x_i = test_flavor_amplitudes.Xi_b_0(), y_j = test_flavor_amplitudes.Xi_b_pr_0())
        test_transition_Xi_0_Xi_prime_0_3 = test_flavor_amplitudes.flavor_matrix_elements(n=3, x_i = test_flavor_amplitudes.Xi_b_0(), y_j = test_flavor_amplitudes.Xi_b_pr_0())
        test_transition_Xi_m_Xi_prime_m_1 = test_flavor_amplitudes.flavor_matrix_elements(n = 1, x_i = test_flavor_amplitudes.Xi_b_m(), y_j = test_flavor_amplitudes.Xi_b_pr_m()) 
        test_transition_Xi_m_Xi_prime_m_2 = test_flavor_amplitudes.flavor_matrix_elements(n=2, x_i = test_flavor_amplitudes.Xi_b_m(), y_j = test_flavor_amplitudes.Xi_b_pr_m()) 
        test_transition_Xi_m_Xi_prime_m_3 = test_flavor_amplitudes.flavor_matrix_elements(n = 3, x_i = test_flavor_amplitudes.Xi_b_m(), y_j = test_flavor_amplitudes.Xi_b_pr_m())
        test_transition_Lambda_0_Lambda_0_1 = test_flavor_amplitudes.flavor_matrix_elements(n = 1, x_i = test_flavor_amplitudes.Lambda_b_0(), y_j = test_flavor_amplitudes.Lambda_b_0())
        test_transition_Lambda_0_Lambda_0_2 = test_flavor_amplitudes.flavor_matrix_elements(n = 2, x_i = test_flavor_amplitudes.Lambda_b_0(), y_j = test_flavor_amplitudes.Lambda_b_0())
        test_transition_Lambda_0_Lambda_0_3 = test_flavor_amplitudes.flavor_matrix_elements(n = 3, x_i = test_flavor_amplitudes.Lambda_b_0(), y_j = test_flavor_amplitudes.Lambda_b_0())
        
        self.assertAlmostEqual(test_transition_lambda_0_sigma_0_1, 0.0008361204013377926)
        self.assertAlmostEqual(test_transition_lambda_0_sigma_0_2, -0.0008361204013377926)
        self.assertAlmostEqual(test_transition_lambda_0_sigma_0_3, 0)
        self.assertAlmostEqual(test_transition_Xi_0_Xi_prime_0_1, 0.0007366250704259119)
        self.assertAlmostEqual(test_transition_Xi_0_Xi_prime_0_2, -0.0007366250704259119)
        self.assertAlmostEqual(test_transition_Xi_0_Xi_prime_0_3, 0)
        self.assertAlmostEqual(test_transition_Xi_m_Xi_prime_m_1, -9.94953309118807 * 10 **-5)
        self.assertAlmostEqual(test_transition_Xi_m_Xi_prime_m_2, 9.94953309118807 * 10 ** -5)
        self.assertAlmostEqual(test_transition_Xi_m_Xi_prime_m_3, 0)
        self.assertAlmostEqual(test_transition_Lambda_0_Lambda_0_1, 0.0002787068004459309)
        self.assertAlmostEqual(test_transition_Lambda_0_Lambda_0_2, 0.0002787068004459309)
        self.assertAlmostEqual(test_transition_Lambda_0_Lambda_0_3, -3.382034632034632 * 10 **-5)

        
if __name__ == '__main__':
    unittest.main()
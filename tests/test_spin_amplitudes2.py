import unittest
import pytest


class TestSpinAmplitudes(unittest.TestCase):
    """
    Tests methods in spin_amplitudes
    """
    @pytest.fixture(autouse=True)
    def initdir(self, tmpdir):
        print(tmpdir)
        tmpdir.chdir()

    def test_matrix_elements(self):
        """
        Test the matrix element method
        """
        from bottomfw.common.spin_amplitudes import SpinAmplitudes
        import numpy as np
        import sympy as sp
        test_amplitudes = SpinAmplitudes(baryons="omegas")
        test_ladder_operator_tensorp = test_amplitudes.ladder_operator_tensor(sign=1, dim1=1, dim2=2)
        test_amp_slo1_s1_2_s3_2 = test_amplitudes.matrix_elements(n=1, x_i=test_amplitudes.symmetric_states(m_proj=1/2), y_j=test_amplitudes.symmetric_states(m_proj=3/2))
        test_amp_slo1_s1_2m_s1_2 = test_amplitudes.matrix_elements(n=1, x_i=test_amplitudes.symmetric_states(m_proj=-1/2), y_j=test_amplitudes.symmetric_states(m_proj=1/2))
        test_amp_slo1_s3_2m_s1_2m = test_amplitudes.matrix_elements(n=1, x_i=test_amplitudes.symmetric_states(m_proj=-3/2), y_j=test_amplitudes.symmetric_states(m_proj=-1/2))
        test_amp_slo1_s1_2m_r1_2 = test_amplitudes.matrix_elements(n=1, x_i=test_amplitudes.symmetric_states(m_proj=-1/2), y_j=test_amplitudes.rho_states(m_proj=1/2))
        test_amp_slo1_s3_2m_r1_2m = test_amplitudes.matrix_elements(n=1, x_i=test_amplitudes.symmetric_states(m_proj=-3/2), y_j=test_amplitudes.rho_states(m_proj=-1/2))
        test_amp_slo1_s1_2m_l1_2 = test_amplitudes.matrix_elements(n=1, x_i=test_amplitudes.symmetric_states(m_proj=-1/2), y_j=test_amplitudes.lambda_states(m_proj=1/2))
        test_amp_slo1_s3_2m_l1_2m = test_amplitudes.matrix_elements(n=1, x_i=test_amplitudes.symmetric_states(m_proj=-3/2), y_j=test_amplitudes.lambda_states(m_proj=-1/2))
        test_amp_slo1_r1_2_s3_2 = test_amplitudes.matrix_elements(n=1, x_i=test_amplitudes.rho_states(m_proj=1/2), y_j=test_amplitudes.symmetric_states(m_proj=3/2))

        self.assertAlmostEqual(test_ladder_operator_tensorp, sp.Matrix([[0, 1, 1, 0, 1, 0, 0, 0], [0, 0, 0, 1, 0, 1, 0, 0], [0, 0, 0, 1, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 1, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 0]]))
        self.assertAlmostEqual(test_amp_slo1_s1_2_s3_2, np.sqrt(3) / 3)
        self.assertAlmostEqual(test_amp_slo1_s1_2m_s1_2, 2 / 3)
        self.assertAlmostEqual(test_amp_slo1_s3_2m_s1_2m, np.sqrt(3) / 3)
        self.assertAlmostEqual(test_amp_slo1_s1_2m_r1_2, np.sqrt(6) / 6)
        self.assertAlmostEqual(test_amp_slo1_s3_2m_r1_2m, np.sqrt(2) / 2)

        self.assertAlmostEqual(test_amp_slo1_s1_2m_l1_2, np.sqrt(2) / 6)
        self.assertAlmostEqual(test_amp_slo1_s3_2m_l1_2m, np.sqrt(6) / 6)
        self.assertAlmostEqual(test_amp_slo1_r1_2_s3_2, -np.sqrt(2) / 2)

        
if __name__ == '__main__':
    unittest.main()


        

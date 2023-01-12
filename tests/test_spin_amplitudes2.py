import unittest
import pytest


class TestBoostedSVM(unittest.TestCase):
    """Tests classes in common_methods"""
    @pytest.fixture(autouse=True)
    def initdir(self, tmpdir):
        print(tmpdir)
        tmpdir.chdir()

    def test_matrix_elements(self):
        """
        Test the matrix element method
        """
        from bottomfw.common.spin_amplitudes import SpinAmplitudes
        test_amplitudes = SpinAmplitudes(baryons="omegas")
        test_amp_slo1_s1_2_s3_2 = test_amplitudes.matrix_elements(n=1, x_i=test_amplitudes.symmetric_states(m_proj=1/2), y_j=test_amplitudes.symmetric_states(m_proj=3/2))
        test_amp_slo1_s1_2m_s1_2 = test_amplitudes.matrix_elements(n=1, x_i=test_amplitudes.symmetric_states(m_proj=-1/2), y_j=test_amplitudes.symmetric_states(m_proj=1/2))
        test_amp_slo1_s3_2m_s1_2m = test_amplitudes.matrix_elements(n=1, x_i=test_amplitudes.symmetric_states(m_proj=-3/2), y_j=test_amplitudes.symmetric_states(m_proj=-1/2))
        test_amp_slo1_s1_2m_r1_2 = test_amplitudes.matrix_elements(n=1, x_i=test_amplitudes.symmetric_states(m_proj=-1/2), y_j=test_amplitudes.rho_states(m_proj=1/2))
        test_amp_slo1_s3_2m_r1_2m = test_amplitudes.matrix_elements(n=1, x_i=test_amplitudes.symmetric_states(m_proj=-3/2), y_j=test_amplitudes.rho_states(m_proj=-1/2))

        self.assertAlmostEqual(test_amp_slo1_s1_2_s3_2, 0.5773502691896257)
        self.assertAlmostEqual(test_amp_slo1_s1_2m_s1_2, 0.6666666666666666)
        self.assertAlmostEqual(test_amp_slo1_s3_2m_s1_2m, 0.5773502691896257)
        self.assertAlmostEqual(test_amp_slo1_s1_2m_r1_2, 0.408248290463863)
        self.assertAlmostEqual(test_amp_slo1_s3_2m_r1_2m, 0.7071067811865476)

        
if __name__ == '__main__':
    unittest.main()


        

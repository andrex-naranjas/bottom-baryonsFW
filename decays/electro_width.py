"""
---------------------------------------------------------------
  Authors: A. Ramirez-Morales (andres.ramirez.morales@cern.ch)
           H. Garcia-Tecocoatzi
---------------------------------------------------------------
"""
from decays.decay_wrapper import decay
import decays.decay_utils as du
import numpy as np


class ElectroWidths:
    """
    Class that administrates the decay width calculations of the hevay baryon widths done by the C++ class
    The class calls the python wrapper and feeds the functions with the needed quatumn numbers
    and masses
    """
    def __init__(self, bootstrap=False, baryons='', workpath="."):
        self.m_width = decay(workpath)
        # self.fetch_decay_masses(bootstrap)
        self.channel_widths_vector = []

    def total_decay_width(self, baryons, k_prim, massA, SA_val, LA_val, JA_val, SL_val, ModEx_val, bootstrap=False, m1=0, m2=0, m3=0):
        """
        Method that calls the wrapper and sums the individual decay widths
        """
        MassA = massA #/1000.0
        SA_qm = SA_val
        LA_qm = LA_val
        JA_qm = JA_val
        SL_qm = SL_val
        baryon= 0 # self.baryon_flag(baryons)
        ModEx = 0 # self.ModEx_flag(ModEx_val)
        nChannels = 0 # self.n_channels(baryons)
        m_lam, m_rho = self.reduced_masses(baryons, m1, m2, m3)
        channel_widths = ([])
        
        alpha_lam = self.alphas(k_prim, m_lam)
        alpha_rho = self.alphas(k_prim, m_rho)

        alpha_lam = 1
        alpha_rho = 1

        MassB = 5.935
        excMode = 1
        prodDecay = 1

        decay_value = self.m_width.electro_width(MassA, MassB,
                                                 SA_qm, LA_qm, JA_qm, SL_qm, alpha_lam, alpha_rho,
                                                 m1, m2, baryon, excMode, prodDecay)


        return decay_value

    def reduced_masses(self, baryons, m1_input, m2_input, m3_input):
        """
        Method to calculate reduced masses of the harmonic oscillator
        """
        m_lam,m_rho=0,0
        if(baryons=='omegas'):
            m_rho = m2_input
            m_lam = (3*m2_input*m1_input)/(2*m2_input+m1_input)
        elif(baryons=='cascades' or baryons =='cascades_anti3'):
            m_rho = (m2_input+m3_input)/2
            m_lam = (1.5*(m2_input+m3_input)*m1_input)/(m1_input+m2_input+m3_input)
        elif(baryons=='sigmas' or baryons=='lambdas'):
             m_rho = m3_input
             m_lam = (3*m3_input*m1_input)/(2*m3_input+m1_input)
        return m_lam, m_rho


    def alphas(self, k_prim, m_lam_rho):
        """
        Method to calculate the decay alphas
        """
        value1 = (np.sqrt(3./m_lam_rho)) * k_prim
        value2 = value1*m_lam_rho
        return np.sqrt(value2)/1000. # transform from GeV -> MeV


test_electro = ElectroWidths()
m1 = 4.928
m2 = 0.299
m3 = 0.465
massA = 6.190
k_prim = 1
SA_val = 0.5
LA_val = 1
JA_val = 0.5
SL_val = 1
ModEx_val = 1
baryons = "cascades"


value = test_electro.total_decay_width(baryons, k_prim, massA, SA_val, LA_val, JA_val, SL_val, ModEx_val, bootstrap=False, m1=m1, m2=m2, m3=m3)

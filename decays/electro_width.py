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

    baryon FLAG: 1 -> omega, 2->cascade_6, 3->sigma,# 4 -> lambda, 5-> cascade_3
    ModEx  FLAG: 0-> ground(grd), 1 -> lambda(lam), 2->rho, 3->rop_lam(rpl), 4->rop_rho(rpr), 5->mix  Excitation
    decPr  FLAG: 0->.....  decayProduct Flag
    """
    def __init__(self, bootstrap=False, baryons='', workpath="."):
        self.m_width = decay(workpath)
        # self.fetch_decay_masses(bootstrap)
        self.channel_widths_vector = []


    def total_decay_width(self, baryons, k_prim, massA, SA_val, JA_val, LA_val, SlA_val, LlA_val, LrA_val,
                          ModEx_val, bootstrap=False, m1=0, m2=0, m3=0):
        """
        Method that calls the wrapper and sums the individual decay widths
        """
        MassA = massA/1000.0
        SA_qm = SA_val
        JA_qm = JA_val
        LA_qm = LA_val
        SlA_qm = SlA_val
        LlA_qm = LlA_val
        LrA_qm = LrA_val

        baryon= 5 # self.baryon_flag(baryons)
        ModEx = 0 # self.ModEx_flag(ModEx_val)
        nChannels = 0 # self.n_channels(baryons)
        m_lam, m_rho = self.reduced_masses(baryons, m1, m2, m3)
        channel_widths = ([])
        
        alpha_lam = self.alphas(k_prim, m_lam)
        alpha_rho = self.alphas(k_prim, m_rho)
        # print(alpha_lam, alpha_rho)

        mbottom  = m1/1000
        mupdown  = m2/1000
        mstrange = m3/1000

        # these depend on each decay
        MassB = 5.935
        SB_qm = 0
        JB_qm = 0
        LB_qm = 0
        SlB_qm = 0
        LlB_qm = 0
        LrB_qm = 0
        
        # first AMPS
        SB_qm = 0.5
        JB_qm = 0.5
        SlB_qm = 0.0

        # Orbital lambda amps
        SB_qm = 0.5
        SlB_qm = 1.0
        LB_qm = 0.0
        LlB_qm = 0.0
        LrB_qm = 0.0
        JB_qm = 0.5

        # Orbital lambda amps second test
        # SB_qm = 0.5 
        # SlB_qm = 1.0
        # LB_qm = 0.0
        # LlB_qm = 0.0
        # LrB_qm = 0.0
        # JB_qm = 0.5

        # hugo test
        # MassB = 5.792
        # SB_qm = 0.5 
        # JB_qm = 0.5
        # LB_qm = 0.0
        # LlB_qm = 0.0
        # LrB_qm = 0.0
        # SlB_qm = 0.0


        excMode = 1
        prodDecay = -1
        decay_value = self.m_width.electro_width(MassA, SA_qm, JA_qm, LA_qm, SlA_qm, LlA_qm, LrA_qm,
                                                 MassB, SB_qm, JB_qm, LB_qm, SlB_qm, LlB_qm, LrB_qm,
                                                 alpha_lam, alpha_rho,
                                                 mbottom, mupdown, mstrange,
                                                 baryon, excMode, prodDecay)
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
        return np.sqrt(value2)/1000. # transform from MeV -> GeV


test_electro = ElectroWidths()
m1 = 4.928 * 1000
m2 = 0.299 * 1000
m3 = 0.465 * 1000
massA = 6.190 * 1000
k_prim = 5044.8
ModEx_val = 1
baryons = "cascades"

# first AMPS
SA_val = 1.5
JA_val = 1.5
LA_val = 1
SlA_val = 1.
LlA_val = 0
LrA_val = 0

# OrbitalLAMBDA amps
SA_val = 0.5
SlA_val = 1.0
LA_val = 1.0
LlA_val = 1.0
LrA_val = 0.0
JA_val = 0.5

# OrbitalLAMBDA amps second test
massA = 6.198 * 1000
SA_val = 0.5
SlA_val = 1.0
LA_val = 1.0
LlA_val = 1.0
LrA_val = 0.0
JA_val = 1.5


# SA_val = 0.5
# SlA_val = 1.0
# LA_val = 1
# LlA_val = 1
# LrA_val = 0
# JA_val = 0.5



value = test_electro.total_decay_width(baryons, k_prim, massA, SA_val, JA_val, LA_val, SlA_val, LlA_val, LrA_val,
                                       ModEx_val, bootstrap=False, m1=m1, m2=m2, m3=m3)
print("EM decay width:  ", value)

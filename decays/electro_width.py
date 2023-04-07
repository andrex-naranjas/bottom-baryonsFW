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
    ModEx  FLAG: 0 -> ground(grd), 1 -> lambda(lam), 2->rho, 3->rop_lam(rpl), 4->rop_rho(rpr), 5->mix  Excitation
    decPr  FLAG: 0 ->.....  decayProduct Flag
    """
    def __init__(self, bootstrap=False, baryons='', workpath="."):
        self.m_width = decay(workpath)
        self.fetch_decay_masses(bootstrap)
        self.channel_widths_vector = []

    def total_decay_width(self, baryons, k_prim, massA, SA_val, JA_val, LA_val, SlA_val, LlA_val, LrA_val,
                          ModEx_val, bootstrap=False, m1=0, m2=0, m3=0):
        """
        Method that calls the wrapper and sums the individual decay widths
        """
        MassA = massA/1000.0
        mbottom  = m1/1000.0
        mupdown  = m2/1000.0
        mstrange = m3/1000.0

        SA_qm = SA_val
        JA_qm = JA_val
        LA_qm = LA_val
        SlA_qm = SlA_val
        LlA_qm = LlA_val
        LrA_qm = LrA_val

        baryon = self.baryon_flag(baryons)
        ModEx = self.ModEx_flag(ModEx_val)
        nChannels = self.n_channels(baryons)
        m_lam, m_rho = self.reduced_masses(baryons, m1, m2, m3)
        channel_widths = ([])
        
        alpha_lam = self.alphas(k_prim, m_lam)
        alpha_rho = self.alphas(k_prim, m_rho)

        for i in range(nChannels):
            decPr = i+1
            decPr = 3
            MassB = self.decay_mass(bootstrap, baryons, decPr)
            single_decay_value = self.m_width.electro_width(MassA, SA_qm, JA_qm, LA_qm, SlA_qm, LlA_qm, LrA_qm,
                                                            MassB,
                                                            alpha_lam, alpha_rho,
                                                            mbottom, mupdown, mstrange,
                                                            baryon, ModEx, decPr)

            channel_widths = np.append(channel_widths, single_decay_value)
            baryon_name, ModEx_name, decPr_name =  "test", "test", "test" # du.state_labels(baryon, ModEx, decPr, LA_qm)
            if not bootstrap:
                print('%6s |  %10s | %12s |  %5.3f |   %5.3f |  %5.1f |  %5.1f |  %5.1f | %5.6f '
                      %(baryon_name, ModEx_name, decPr_name, MassA, MassB, JA_qm, LA_qm, SA_qm, single_decay_value))
                    
        # sum the individual width to obtain total width
        total_decay_width = np.sum(channel_widths)
        # print(alpha_lam,alpha_rho)
        if not bootstrap:
            print('          ******************   TOTAL WIDTH FOR', baryons, ModEx_name, round(total_decay_width,4), '   ******************')
            print('-------------------------------------------------------------------------------------------------------------')
            
        self.channel_widths_vector.append(channel_widths) # for individual decay tables, this is a list of arrays!
        return total_decay_width


    def baryon_flag(self, baryons):
        """
        Method to parse the baryons names to integers
        """
        if(baryons=='omegas'):           return 1
        elif(baryons=='cascades'):       return 2
        elif(baryons=='sigmas'):         return 3
        elif(baryons=='lambdas'):        return 4
        elif(baryons=='cascades_anti3'): return 5

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

    def ModEx_flag(self, ModEx_val):
        """
        Method to parse the h.o mode to integers
        grd=0, lam =1 , rho=2, rpl=3, rpr=4, mix=5
        """
        if(ModEx_val=='grd'):   return 0
        elif(ModEx_val=='lam'): return 1
        elif(ModEx_val=='rho'): return 2
        elif(ModEx_val=='rpl'): return 3
        elif(ModEx_val=='rpr'): return 4
        elif(ModEx_val=='mix'): return 5

    def n_channels(self, baryons):
        """
        Method to set number of decay channels has each baryon
        """
        if(baryons=='omegas'):           return 2
        elif(baryons=='cascades'):       return 6
        elif(baryons=='sigmas'):         return 7
        elif(baryons=='lambdas'):        return 3
        elif(baryons=='cascades_anti3'): return 6

    def decay_mass(self, bootstrap, baryons, decPr):
        """
        Method to fetch mass of the decay products
        """
        if(baryons=='omegas'):
            if(decPr==1):
                if not bootstrap:  return self.omega_mass
                else: return np.random.choice(self.gauss_omega, size=None)
            if(decPr==2):
                if not bootstrap:  return self.omega_s_mass
                else: return np.random.choice(self.gauss_omega_s, size=None)         
        elif(baryons=='cascades'):
            if(decPr==1):
                if not bootstrap: return self.xi_mass
                else: return np.random.choice(self.gauss_xi)
            elif(decPr==2):
                if not bootstrap: return self.xi_mass
                else: return np.random.choice(self.gauss_xi, size=None)
            elif(decPr==3):
                if not bootstrap: return self.xi_p_mass
                else: return np.random.choice(self.gauss_xi_p, size=None)
            elif(decPr==4):
                if not bootstrap: return self.xi_p_mass
                else: return np.random.choice(self.gauss_xi_p, size=None)
            elif(decPr==5):
                if not bootstrap: return self.xi_p_s_mass
                else: return np.random.choice(self.gauss_xi_p_s, size=None)
            elif(decPr==6):
                if not bootstrap: return self.xi_p_s_mass
                else: return np.random.choice(self.gauss_xi_p_s, size=None)
        elif(baryons=='sigmas'):
            if(decPr==1):
                if not bootstrap: return self.sigma_mass
                else: return np.random.choice(self.gauss_sigma, size=None)
            elif(decPr==2):
                if not bootstrap: return self.sigma_mass
                else: return np.random.choice(self.gauss_sigma, size=None)
            elif(decPr==3):
                if not bootstrap: return self.sigma_mass
                else: return np.random.choice(self.gauss_sigma, size=None)
            elif(decPr==4):
                if not bootstrap: return self.lambda_mass
                else: return np.random.choice(self.gauss_lambda, size=None)
            elif(decPr==5):
                if not bootstrap: return self.sigma_s_mass
                else: return np.random.choice(self.gauss_sigma_s, size=None)
            elif(decPr==6):
                if not bootstrap: return self.sigma_s_mass
                else: return np.random.choice(self.gauss_sigma_s, size=None)
            elif(decPr==7):
                if not bootstrap: return self.sigma_s_mass
                else: return np.random.choice(self.gauss_sigma_s, size=None)
        elif(baryons=='lambdas'):
            if(decPr==1):
                if not bootstrap: return self.lambda_mass
                else: return np.random.choice(self.gauss_lambda, size=None)
            elif(decPr==2):
                if not bootstrap: return self.sigma_mass
                else: return np.random.choice(self.gauss_sigma, size=None)
            elif(decPr==3):
                if not bootstrap: return self.sigma_s_mass
                else: return np.random.choice(self.gauss_sigma_s, size=None)
        elif(baryons=='cascades_anti3'):
            if(decPr==1):
                if not bootstrap: return self.xi_mass,     self.pion_mass
                else: return np.random.choice(self.gauss_xi, size=None)
            elif(decPr==2):
                if not bootstrap: return self.xi_mass,     self.pion_mass
                else: return np.random.choice(self.gauss_xi, size=None)
            elif(decPr==3):
                if not bootstrap: return self.xi_p_mass,   self.pion_mass
                else: return np.random.choice(self.gauss_xi_p, size=None)
            elif(decPr==4):
                if not bootstrap: return self.xi_p_mass,   self.pion_mass
                else: return np.random.choice(self.gauss_xi_p, size=None)
            elif(decPr==5):
                if not bootstrap: return self.xi_p_s_mass
                else: return np.random.choice(self.gauss_xi_p_s, size=None)
            elif(decPr==6):
                if not bootstrap: return self.xi_p_s_mass
                else: return np.random.choice(self.gauss_xi_p_s, size=None)

    def fetch_decay_masses(self, bootstrap):
        # Bottom hadrons
        self.lambda_mass      = 5.61960 # +- 0.0001
        self.xi_p_mass        = 5.9350#2 # +- 0.00005        
        self.xi_p_s_mass      = 5.9350#2 # +- 0.00005        CHECK!!
        self.xi_mass          = 5.79700 # +- 0.00060.... Difference with Xb0=5.9 +- 0.6 MeV
        self.xi_s_mass        = 6.07800 # +- 0.00006 (predicted mass)$6078^{+10}_{-10}$  CHECK!!
        self.sigma_mass       = 5.81056 # +- 0.00025.... Difference of + and - == 5.06+-0.18 MeV
        self.sigma_s_mass     = 5.83032 # +- 0.00030.... Difference of + and - == 4.37+-0.33 OK
        self.omega_mass       = 6.04520 # +- 0.00120
        self.omega_s_mass     = 6.09300 # +- 0.00060 (predicted mass) # $6093^{+10}_{-10}$ CHECK!!
        self.B0_mass          = 5.27966 # +- 0.00012
        self.Bs_mass          = 5.36692 # +- 0.00010
        self.B_star_mass      = 5.32471 # +- 0.00021

        if(bootstrap):
            # Bottom hadrons
            self.gauss_lambda      = np.random.normal(5.61960, 0.00017, 10000)
            self.gauss_xi_p        = np.random.normal(5.93502, 0.00005, 10000)
            self.gauss_xi          = np.random.normal(5.79700, 0.00060, 10000)
            self.gauss_xi_s        = np.random.normal(6.07800, 0.00100, 10000) # predicted massA
            self.gauss_sigma       = np.random.normal(5.81056, 0.00025, 10000)
            self.gauss_sigma_s     = np.random.normal(5.83032, 0.00030, 10000)
            self.gauss_omega       = np.random.normal(6.04520, 0.00120, 10000)
            self.gauss_omega_s     = np.random.normal(6.09300, 0.00060, 10000) # predicted massA    
            self.gauss_B0          = np.random.normal(5.27966, 0.00012, 10000)
            self.gauss_Bs          = np.random.normal(5.36692, 0.00010, 10000)
            self.gauss_B_star      = np.random.normal(5.32471, 0.00021, 10000)




test_electro = ElectroWidths()
m1 = 4.928 * 1000
m2 = 0.299 * 1000
m3 = 0.465 * 1000
massA = 6.190 * 1000
k_prim = 5044.8
ModEx_val = "lam"
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


        # SB_qm = 0
        # JB_qm = 0
        # LB_qm = 0
        # SlB_qm = 0
        # LlB_qm = 0
        # LrB_qm = 0
        
        # # first AMPS
        # SB_qm = 0.5
        # JB_qm = 0.5
        # SlB_qm = 0.0

        # # Orbital lambda amps
        # SB_qm = 0.5
        # SlB_qm = 1.0
        # LB_qm = 0.0
        # LlB_qm = 0.0
        # LrB_qm = 0.0
        # JB_qm = 0.5

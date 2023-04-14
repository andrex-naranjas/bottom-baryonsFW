"""
---------------------------------------------------------------
  Authors: A. Ramirez-Morales (andres.ramirez.morales@cern.ch)
           H. Garcia-Tecocoatzi
---------------------------------------------------------------
"""
from decays.decay_wrapper import decay
import decays.decay_utils_em as du
import numpy as np


class ElectroWidths:
    """
    Class that administrates the decay width calculations of the hevay baryon widths done by the C++ class
    The class calls the python wrapper and feeds the functions with the needed quatumn numbers
    and masses

    baryon FLAG: 1 -> omega, 2->cascade_6, 3->sigma,# 4 -> lambda, 5-> cascade_3
    ModEx  FLAG: 0 -> ground(grd), 1 -> lambda(lam), 2->rho
    decPr  FLAG: 0 -> ...  decayProduct Flag
    """
    def __init__(self, bootstrap=False, baryons='', workpath="."):
        self.m_width = decay(workpath)
        self.fetch_decay_masses(bootstrap)
        self.channel_widths_vector = []

    def total_decay_width(self, baryons, k_prim, massA, SA_val, JA_val, LA_val, SlA_val,
                          ModEx_val, bootstrap=False, m1=0, m2=0, m3=0):
        """
        Method that calls the wrapper and sums the individual decay widths
        """
        mb = m1
        ms = m2
        ml = m3
        
        MassA = massA/1000.0
        mbottom  = mb/1000.0
        mupdown  = ml/1000.0
        mstrange = ms/1000.0

        SA_qm = SA_val
        JA_qm = JA_val
        LA_qm = LA_val
        SlA_qm = SlA_val
        LlA_qm, LrA_qm = self.orbital_projections(ModEx_val, LA_val)
            
        baryon = self.baryon_flag(baryons)
        ModEx = self.ModEx_flag(ModEx_val)
        nChannels = self.n_channels(baryons)
        m_lam, m_rho = self.reduced_masses(baryons, mb, ml, ms)
        channel_widths = ([])
        
        alpha_lam = self.alphas(k_prim, m_lam)
        alpha_rho = self.alphas(k_prim, m_rho)

        for i in range(nChannels):
            decPr = i+1
            MassB = self.decay_mass(bootstrap, baryons, decPr)
            single_decay_value = self.m_width.electro_width(MassA, SA_qm, JA_qm, LA_qm, SlA_qm, LlA_qm, LrA_qm,
                                                            MassB,
                                                            alpha_lam, alpha_rho,
                                                            mbottom, mupdown, mstrange,
                                                            baryon, ModEx, decPr)
            channel_widths = np.append(channel_widths, single_decay_value)
            baryon_name, ModEx_name, decPr_name =  du.state_labels(baryon, ModEx, decPr, LA_qm)
            if not bootstrap:
                print('%6s |  %10s | %12s |  %5.3f |   %5.3f |  %5.1f |  %5.1f |  %5.1f | %5.1f | %5.6f '
                      %(baryon_name, ModEx_name, decPr_name, MassA, MassB, JA_qm, LA_qm, SA_qm, SlA_qm, single_decay_value))
                    
        # sum the individual width to obtain total width
        total_decay_width = np.sum(channel_widths)
        # print(alpha_lam,alpha_rho)
        if not bootstrap:
            print(' ******************   TOTAL ELECTROMAGNETIC WIDTH FOR', baryons, ModEx_name, round(total_decay_width,4), '   ******************')
            print('-------------------------------------------------------------------------------------------------------------')
            
        self.channel_widths_vector.append(channel_widths) # for individual decay tables, this is a list of arrays!
        return total_decay_width


    def orbital_projections(self, ModEx_val, LA_val):
        """
        Method to fecth the orbital projection (up to P-wave)
        """
        if(ModEx_val=="grd"):
            LlA = LA_val
            LrA = LA_val
        elif(ModEx_val=="lam"):
            LlA = LA_val
            LrA = 0
        elif(ModEx_val=="rho"):
            LlA = 0
            LrA = LA_val
        else:
            LlA = -1
            LrA = -1
        return LlA, LrA

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
        grd=0, lam =1 , rho=2
        """
        if(ModEx_val=='grd'):   return 0
        elif(ModEx_val=='lam'): return 1
        elif(ModEx_val=='rho'): return 2

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
        elif(baryons=='cascades' or baryons=='cascades_anti3'):
            if(decPr==1):
                if not bootstrap: return self.xi_mass
                else: return np.random.choice(self.gauss_xi, size=None)
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
        
    def fetch_decay_masses(self, bootstrap):
        '''
        Method to fetch the decay products coming from our fit (mA)
        '''
        self.omega_mass    = 6.06400
        self.omega_s_mass  = 6.09300
        self.sigma_mass    = 5.80500
        self.sigma_s_mass  = 5.83400
        self.xi_p_mass     = 5.92500
        self.xi_p_s_mass   = 5.95500
        self.xi_mass       = 5.80600
        self.lambda_mass   = 5.61400
       
        if(bootstrap):
            self.gauss_omega    = np.random.normal(6.06400, 0.00600, 10000)
            self.gauss_omega_s  = np.random.normal(6.09300, 0.00700, 10000)
            self.gauss_sigma    = np.random.normal(5.80500, 0.00600, 10000)
            self.gauss_sigma_s  = np.random.normal(5.83400, 0.00700, 10000)
            self.gauss_xi_p     = np.random.normal(5.92500, 0.00500, 10000)
            self.gauss_xi_p_s   = np.random.normal(5.95500, 0.00500, 10000)
            self.gauss_xi       = np.random.normal(5.80600, 0.00700, 10000)
            self.gauss_lambda   = np.random.normal(5.61400, 0.00700, 10000)
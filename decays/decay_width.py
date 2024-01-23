"""
---------------------------------------------------------------
  Authors: A. Ramirez-Morales (andres.ramirez.morales@cern.ch)
           H. Garcia-Tecocoatzi
---------------------------------------------------------------
"""
from decays.decay_wrapper import decay
import decays.decay_utils as du
import numpy as np


class DecayWidths:
    """
    Class that administrates the decay width calculations of the hevay baryon widths done by the C++ class
    The class calls the python wrapper and feeds the functions with the needed quatumn numbers
    and masses. 

    baryon FLAG: 1 -> omega, 2->cascade_6, 3->sigma,# 4 -> lambda, 5-> cascade_3
    ModEx  FLAG: 0-> ground(grd), 1 -> lambda(lam), 2->rho, 3->rop_lam(rpl), 4->rop_rho(rpr), 5->mix  Excitation
    decPr  FLAG: 3->Xi'+Pi, 5->Sigma+K  decayProduct Flag
    """
    def __init__(self, bootstrap=False, baryons='', workpath="."):
        self.m_width = decay(workpath)
        self.fetch_decay_masses(bootstrap)
        self.set_gamma_val(bootstrap)
        self.channel_widths_vector = []

    def load_average_mass(self, mass_avg=0):
        self.MassA_avg = mass_avg/1000.0 # central value of the bootstrap distribution
        
    def total_decay_width(self, baryons, k_prim, massA, SA_val, L_val, JA_val, SL_val, ModEx_val, bootstrap=False, gamma_val=None, m1=0, m2=0, m3=0):
        """
        Method that calls the wrapper and sums the individual decay widths
        """
        MassA = massA/1000.0
        MassA_avg = self.MassA_avg
        SA_qm = SA_val
        LA_qm = L_val
        JA_qm = JA_val
        SL_qm = SL_val
        baryon= self.baryon_flag(baryons)
        ModEx = self.ModEx_flag(ModEx_val)
        nChannels = self.n_channels(baryons)
        m_lam, m_rho = self.reduced_masses(baryons, m1, m2, m3)
        channel_widths = ([])
        
        alpha_lam = self.alphas(k_prim, m_lam)
        alpha_rho = self.alphas(k_prim, m_rho)

        if gamma_val is None:
            gamma = self.gamma_fitted(bootstrap)
        else:
            gamma = gamma_val

        for i in range(nChannels):
            decPr = i+1
            MassB,MassC = self.decay_masses(bootstrap, baryons, decPr)
            MassB_avg,MassC_avg = self.decay_masses(False, baryons, decPr)
            single_decay_value = self.m_width.decay_width(MassA_avg, MassB_avg, MassC_avg, MassA, MassB, MassC,
                                                          gamma, SA_qm,
                                                          LA_qm, JA_qm, SL_qm, alpha_lam, alpha_rho,
                                                          baryon, ModEx, decPr)
            if (MassA_avg<MassB_avg+MassC_avg): # energy conservation
                single_decay_value = 0.
            channel_widths = np.append(channel_widths, single_decay_value)
            baryon_name, ModEx_name, decPr_name = du.state_labels(baryon, ModEx, decPr, LA_qm)
            if not bootstrap:
                print('%6s |  %10s | %12s |  %5.3f |  %5.3f | %5.3f |  %5.1f |  %5.1f |  %5.1f |  %5.1f | %5.6f '
                      %(baryon_name, ModEx_name, decPr_name, MassA, MassB, MassC, JA_qm, LA_qm, SA_qm, SL_qm, single_decay_value))
                    
        # sum the individual width to obtain total width
        total_decay_width = np.sum(channel_widths)
        # print(alpha_lam,alpha_rho)
        if not bootstrap:
            print('          ******************   TOTAL WIDTH FOR', baryons, ModEx_name, round(total_decay_width,4), '   ******************')
            print('-------------------------------------------------------------------------------------------------------------')
            
        self.channel_widths_vector.append(channel_widths) # for individual decay tables, this is a list of arrays!
        return total_decay_width


    def gamma_fitted(self, bootstrap=False):
        """
        gamma=17.25/np.sqrt(7) # old value
        """
        if not bootstrap: return self.gamma_fit
        else:
            return np.random.choice(self.gauss_gamma, size=None)
        
    def baryon_flag(self, baryons):
        """
        Method to parse the baryons names to integers
        """
        if(baryons=='omegas'):           return 1
        elif(baryons=='cascades'):       return 2
        elif(baryons=='sigmas'):         return 3
        elif(baryons=='lambdas'):        return 4
        elif(baryons=='cascades_anti3'): return 5
        
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
        if(baryons=='omegas'):           return 12  +2
        elif(baryons=='cascades'):       return 7+17+4
        elif(baryons=='sigmas'):         return 5+13+8
        elif(baryons=='lambdas'):        return 3+10+1
        elif(baryons=='cascades_anti3'): return 7+17+4        
        
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
             
        return m_lam,m_rho
                    
    def alphas(self, k_prim, m_lam_rho):
        """
        Method to calculate the decay alphas
        """
        value1 = (np.sqrt(3./m_lam_rho)) * k_prim
        value2 = value1*m_lam_rho
        return np.sqrt(value2)/1000. # transform from GeV -> MeV
            
    def decay_masses(self, bootstrap, baryons, decPr):
        """
        Method to fetch mass of the decay products
        """
        if(baryons=='omegas'):
            if(decPr==1):
                if not bootstrap: return self.xi_mass,   self.kaon_mass
                else: return np.random.choice(self.gauss_xi, size=None), np.random.choice(self.gauss_kaon, size=None)
            elif(decPr==2):                
                if not bootstrap: return self.xi_p_mass, self.kaon_mass
                else: return np.random.choice(self.gauss_xi_p, size=None), np.random.choice(self.gauss_kaon, size=None)
            elif(decPr==3):
                if not bootstrap: return self.xi_s_mass, self.kaon_mass
                else: return np.random.choice(self.gauss_xi_s, size=None), np.random.choice(self.gauss_kaon, size=None)
            elif(decPr==4):
                if not bootstrap: return self.xi_mass, self.kaon_s_mass
                else: return np.random.choice(self.gauss_xi, size=None), np.random.choice(self.gauss_kaon_s, size=None)
            elif(decPr==5):
                if not bootstrap: return self.xi_p_mass, self.kaon_s_mass
                else: return np.random.choice(self.gauss_xi_p, size=None), np.random.choice(self.gauss_kaon_s, size=None)
            elif(decPr==6):
                if not bootstrap: return self.xi_s_mass, self.kaon_s_mass
                else: return np.random.choice(self.gauss_xi_s, size=None), np.random.choice(self.gauss_kaon_s, size=None)
            elif(decPr==7):
                if not bootstrap: return self.omega_mass, self.eta_mass
                else: return np.random.choice(self.gauss_omega, size=None), np.random.choice(self.gauss_eta, size=None)
            elif(decPr==8):
                if not bootstrap: return self.omega_s_mass, self.eta_mass
                else: return np.random.choice(self.gauss_omega_s, size=None), np.random.choice(self.gauss_eta, size=None)
            elif(decPr==9):
                if not bootstrap: return self.omega_mass, self.phi_mass
                else: return np.random.choice(self.gauss_omega, size=None), np.random.choice(self.gauss_phi, size=None)
            elif(decPr==10):
                if not bootstrap: return self.omega_s_mass, self.phi_mass
                else: return np.random.choice(self.gauss_omega_s, size=None), np.random.choice(self.gauss_phi, size=None)
            elif(decPr==11):
                if not bootstrap: return self.omega_mass, self.eta_p_mass
                else: return np.random.choice(self.gauss_omega, size=None), np.random.choice(self.gauss_eta_p, size=None)
            elif(decPr==12):
                if not bootstrap: return self.omega_s_mass, self.eta_p_mass
                else: return np.random.choice(self.gauss_omega_s, size=None), np.random.choice(self.gauss_eta_p, size=None)
            # START OF SECOND DIAGRAM DECAYS
            elif(decPr==13):
                if not bootstrap: return self.Xi_8_mass, self.B0_mass
                else: return np.random.choice(self.gauss_Xi_8, size=None), np.random.choice(self.gauss_B0, size=None)
            elif(decPr==14):
                if not bootstrap: return self.Xi_10_mass, self.B0_mass
                else: return np.random.choice(self.gauss_Xi_10, size=None), np.random.choice(self.gauss_B0, size=None)

        elif(baryons=='cascades'):
            if(decPr==1):
                if not bootstrap: return self.lambda_mass, self.kaon_mass
                else: return np.random.choice(self.gauss_lambda, size=None), np.random.choice(self.gauss_kaon, size=None)
            elif(decPr==2):
                if not bootstrap: return self.xi_mass,     self.pion_mass
                else: return np.random.choice(self.gauss_xi, size=None), np.random.choice(self.gauss_pion, size=None)
            elif(decPr==3):
                if not bootstrap: return self.xi_p_mass,   self.pion_mass
                else: return np.random.choice(self.gauss_xi_p, size=None), np.random.choice(self.gauss_pion, size=None)
            elif(decPr==4):
                if not bootstrap: return self.xi_s_mass,   self.pion_mass
                else: return np.random.choice(self.gauss_xi_s, size=None), np.random.choice(self.gauss_pion, size=None)
            elif(decPr==5):
                if not bootstrap: return self.sigma_mass,  self.kaon_mass
                else: return np.random.choice(self.gauss_sigma, size=None), np.random.choice(self.gauss_kaon, size=None)
            elif(decPr==6):
                if not bootstrap: return self.sigma_s_mass,self.kaon_mass
                else: return np.random.choice(self.gauss_sigma_s, size=None), np.random.choice(self.gauss_kaon, size=None)
            elif(decPr==7):
                if not bootstrap: return self.xi_mass,     self.eta_mass
                else: return np.random.choice(self.gauss_xi, size=None), np.random.choice(self.gauss_eta, size=None)
            elif(decPr==8):
                if not bootstrap: return self.lambda_mass,     self.kaon_s_mass
                else: return np.random.choice(self.gauss_lambda, size=None), np.random.choice(self.gauss_kaon_s, size=None)
            elif(decPr==9):
                if not bootstrap: return self.xi_mass,     self.meson_rho_mass
                else: return np.random.choice(self.gauss_xi, size=None), np.random.choice(self.gauss_meson_rho, size=None)
            elif(decPr==10):
                if not bootstrap: return self.xi_p_mass,     self.meson_rho_mass
                else: return np.random.choice(self.gauss_xi_p, size=None), np.random.choice(self.gauss_meson_rho, size=None)
            elif(decPr==11):
                if not bootstrap: return self.xi_s_mass,     self.meson_rho_mass
                else: return np.random.choice(self.gauss_xi_s, size=None), np.random.choice(self.gauss_meson_rho, size=None)
            elif(decPr==12):
                if not bootstrap: return self.sigma_mass,     self.kaon_s_mass
                else: return np.random.choice(self.gauss_sigma, size=None), np.random.choice(self.gauss_kaon_s, size=None)
            elif(decPr==13):
                if not bootstrap: return self.sigma_s_mass,     self.kaon_s_mass
                else: return np.random.choice(self.gauss_sigma_s, size=None), np.random.choice(self.gauss_kaon_s, size=None)
            elif(decPr==14):
                if not bootstrap: return self.xi_p_mass,     self.eta_mass
                else: return np.random.choice(self.gauss_xi_p, size=None), np.random.choice(self.gauss_eta, size=None)
            elif(decPr==15):
                if not bootstrap: return self.xi_s_mass,     self.eta_mass
                else: return np.random.choice(self.gauss_xi_s, size=None), np.random.choice(self.gauss_eta, size=None)
            elif(decPr==16):
                if not bootstrap: return self.xi_mass,     self.eta_p_mass
                else: return np.random.choice(self.gauss_xi, size=None), np.random.choice(self.gauss_eta_p, size=None)
            elif(decPr==17):
                if not bootstrap: return self.xi_p_mass,     self.eta_p_mass
                else: return np.random.choice(self.gauss_xi_p, size=None), np.random.choice(self.gauss_eta_p, size=None)
            elif(decPr==18):
                if not bootstrap: return self.xi_s_mass,     self.eta_p_mass
                else: return np.random.choice(self.gauss_xi_s, size=None), np.random.choice(self.gauss_eta_p, size=None)
            elif(decPr==19):
                if not bootstrap: return self.xi_mass,     self.small_omega_mass
                else: return np.random.choice(self.gauss_xi, size=None), np.random.choice(self.gauss_small_omega, size=None)
            elif(decPr==20):
                if not bootstrap: return self.xi_p_mass,     self.small_omega_mass
                else: return np.random.choice(self.gauss_xi_p, size=None), np.random.choice(self.gauss_small_omega, size=None)
            elif(decPr==21):
                if not bootstrap: return self.xi_s_mass,     self.small_omega_mass
                else: return np.random.choice(self.gauss_xi_s, size=None), np.random.choice(self.gauss_small_omega, size=None)
            elif(decPr==22):
                if not bootstrap: return self.xi_mass,     self.phi_mass
                else: return np.random.choice(self.gauss_xi, size=None), np.random.choice(self.gauss_phi, size=None)
            elif(decPr==23):
                if not bootstrap: return self.xi_p_mass,     self.phi_mass
                else: return np.random.choice(self.gauss_xi_p, size=None), np.random.choice(self.gauss_phi, size=None)
            elif(decPr==24):
                if not bootstrap: return self.xi_s_mass,     self.phi_mass
                else: return np.random.choice(self.gauss_xi_s, size=None), np.random.choice(self.gauss_phi, size=None)
                # START OF SECOND DIAGRAM DECAYS
            elif(decPr==25):
                if not bootstrap: return self.Sigma_8_mass, self.B0_mass #//Sigma_8_lam+B
                else: return np.random.choice(self.gauss_Sigma_8, size=None), np.random.choice(self.gauss_B0, size=None)
            elif(decPr==26):
                if not bootstrap: return self.Xi_8_mass, self.Bs_mass # //Xi_8_lam+Bs
                else: return np.random.choice(self.gauss_Xi_8, size=None), np.random.choice(self.gauss_Bs, size=None)
            elif(decPr==27):
                if not bootstrap: return self.Sigma_8_mass, self.B_star_mass #  //Sigma_8_lam+B*
                else: return np.random.choice(self.gauss_Sigma_8, size=None), np.random.choice(self.gauss_B_star, size=None)                
            elif(decPr==28):
                if not bootstrap: return self.Sigma_10_mass, self.B0_mass #  //Sigma_10_lam+B
                else: return np.random.choice(self.gauss_Sigma_10, size=None), np.random.choice(self.gauss_B0, size=None)
                
        elif(baryons=='sigmas'):
            if(decPr==1):
                if not bootstrap: return self.sigma_mass,  self.pion_mass
                else: return np.random.choice(self.gauss_sigma, size=None), np.random.choice(self.gauss_pion, size=None)
            elif(decPr==2):
                if not bootstrap: return self.sigma_s_mass, self.pion_mass
                else: return np.random.choice(self.gauss_sigma_s, size=None), np.random.choice(self.gauss_pion, size=None)
            elif(decPr==3):
                if not bootstrap: return self.lambda_mass,  self.pion_mass
                else: return np.random.choice(self.gauss_lambda, size=None), np.random.choice(self.gauss_pion, size=None)
            elif(decPr==4):
                if not bootstrap: return self.sigma_mass,   self.eta_mass
                else: return np.random.choice(self.gauss_sigma, size=None), np.random.choice(self.gauss_eta, size=None)
            elif(decPr==5):
                if not bootstrap: return self.xi_mass,      self.kaon_mass
                else: return np.random.choice(self.gauss_xi, size=None), np.random.choice(self.gauss_kaon, size=None)
            elif(decPr==6):
                if not bootstrap: return self.sigma_mass,      self.meson_rho_mass
                else: return np.random.choice(self.gauss_sigma, size=None), np.random.choice(self.gauss_meson_rho, size=None)
            elif(decPr==7):
                if not bootstrap: return self.sigma_s_mass,      self.meson_rho_mass
                else: return np.random.choice(self.gauss_sigma_s, size=None), np.random.choice(self.gauss_meson_rho, size=None)
            elif(decPr==8):
                if not bootstrap: return self.lambda_mass,      self.meson_rho_mass
                else: return np.random.choice(self.gauss_lambda, size=None), np.random.choice(self.gauss_meson_rho, size=None)
            elif(decPr==9):
                if not bootstrap: return self.sigma_s_mass,      self.eta_mass
                else: return np.random.choice(self.gauss_sigma_s, size=None), np.random.choice(self.gauss_eta, size=None)
            elif(decPr==10):
                if not bootstrap: return self.sigma_mass,      self.eta_p_mass
                else: return np.random.choice(self.gauss_sigma, size=None), np.random.choice(self.gauss_eta_p, size=None)
            elif(decPr==11):
                if not bootstrap: return self.sigma_s_mass,      self.eta_p_mass
                else: return np.random.choice(self.gauss_sigma_s, size=None), np.random.choice(self.gauss_eta_p, size=None)
            elif(decPr==12):
                if not bootstrap: return self.xi_p_mass,      self.kaon_mass
                else: return np.random.choice(self.gauss_xi_p, size=None), np.random.choice(self.gauss_kaon, size=None)
            elif(decPr==13):
                if not bootstrap: return self.xi_s_mass,      self.kaon_mass
                else: return np.random.choice(self.gauss_xi_s, size=None), np.random.choice(self.gauss_kaon, size=None)
            elif(decPr==14):
                if not bootstrap: return self.xi_mass,      self.kaon_s_mass
                else: return np.random.choice(self.gauss_xi, size=None), np.random.choice(self.gauss_kaon_s, size=None)
            elif(decPr==15):
                if not bootstrap: return self.xi_p_mass,      self.kaon_s_mass
                else: return np.random.choice(self.gauss_xi_p, size=None), np.random.choice(self.gauss_kaon_s, size=None)
            elif(decPr==16):
                if not bootstrap: return self.xi_s_mass,      self.kaon_s_mass
                else: return np.random.choice(self.gauss_xi_s, size=None), np.random.choice(self.gauss_kaon_s, size=None)
            elif(decPr==17):
                if not bootstrap: return self.sigma_mass,      self.small_omega_mass
                else: return np.random.choice(self.gauss_sigma, size=None), np.random.choice(self.gauss_small_omega, size=None)
            elif(decPr==18):
                if not bootstrap: return self.sigma_s_mass,      self.small_omega_mass
                else: return np.random.choice(self.gauss_sigma_s, size=None), np.random.choice(self.gauss_small_omega, size=None)
                # START OF SECOND DIAGRAM DECAYS
            elif(decPr==19):
                if not bootstrap: return self.proton_mass, self.B0_mass #  //N+D
                else: return np.random.choice(self.gauss_proton, size=None), np.random.choice(self.gauss_B0, size=None)
            elif(decPr==20):
                if not bootstrap: return self.Sigma_8_mass, self.Bs_mass # //Sigma_8+Ds
                else: return np.random.choice(self.gauss_Sigma_8, size=None), np.random.choice(self.gauss_Bs, size=None)
            elif(decPr==21):
                if not bootstrap: return self.proton_mass, self.B_star_mass # //N+D*
                else: return np.random.choice(self.gauss_proton, size=None), np.random.choice(self.gauss_B_star, size=None)
            elif(decPr==22):
                if not bootstrap: return self.Delta_mass, self.B0_mass #  //Delta+D    
                else: return np.random.choice(self.gauss_Delta, size=None), np.random.choice(self.gauss_B0, size=None)
            elif(decPr==23):
                if not bootstrap: return self.N_star_1520_mass, self.B0_mass # //N*(1520)+D
                else: return np.random.choice(self.gauss_N_star_1520, size=None), np.random.choice(self.gauss_B0, size=None)
            elif(decPr==24):
                if not bootstrap: return self.N_star_1535_mass, self.B0_mass #  //N*(1535)+D    
                else: return np.random.choice(self.gauss_N_star_1535, size=None), np.random.choice(self.gauss_B0, size=None)
            elif(decPr==25):
                if not bootstrap: return self.N_star_1680_mass, self.B0_mass # //N*(1680)+D
                else: return np.random.choice(self.gauss_N_star_1680, size=None), np.random.choice(self.gauss_B0, size=None)
            elif(decPr==26):
                if not bootstrap: return self.N_star_1720_mass, self.B0_mass # //N*(1720)+D
                else: return np.random.choice(self.gauss_N_star_1720, size=None), np.random.choice(self.gauss_B0, size=None)

        elif(baryons=='lambdas'):
            if(decPr==1):
                if not bootstrap: return self.sigma_mass,   self.pion_mass
                else: return np.random.choice(self.gauss_sigma, size=None), np.random.choice(self.gauss_pion, size=None)
            elif(decPr==2):
                if not bootstrap: return self.sigma_s_mass, self.pion_mass
                else: return np.random.choice(self.gauss_sigma_s, size=None), np.random.choice(self.gauss_pion, size=None)
            elif(decPr==3):
                if not bootstrap: return self.lambda_mass,  self.eta_mass
                else: return np.random.choice(self.gauss_lambda, size=None), np.random.choice(self.gauss_eta, size=None)
            elif(decPr==4):
                if not bootstrap: return self.sigma_mass,  self.meson_rho_mass
                else: return np.random.choice(self.gauss_sigma, size=None), np.random.choice(self.gauss_meson_rho, size=None)
            elif(decPr==5):
                if not bootstrap: return self.sigma_s_mass,  self.meson_rho_mass
                else: return np.random.choice(self.gauss_sigma_s, size=None), np.random.choice(self.gauss_meson_rho, size=None)
            elif(decPr==6):
                if not bootstrap: return self.lambda_mass,  self.eta_p_mass
                else: return np.random.choice(self.gauss_lambda, size=None), np.random.choice(self.gauss_eta_p, size=None)
            elif(decPr==7):
                if not bootstrap: return self.lambda_mass,  self.small_omega_mass
                else: return np.random.choice(self.gauss_lambda, size=None), np.random.choice(self.gauss_small_omega, size=None)
            elif(decPr==8):
                if not bootstrap: return self.xi_mass,  self.kaon_mass
                else: return np.random.choice(self.gauss_xi, size=None), np.random.choice(self.gauss_kaon, size=None)
            elif(decPr==9):
                if not bootstrap: return self.xi_p_mass,  self.kaon_mass
                else: return np.random.choice(self.gauss_xi_p, size=None), np.random.choice(self.gauss_kaon, size=None)
            elif(decPr==10):
                if not bootstrap: return self.xi_s_mass,  self.kaon_mass
                else: return np.random.choice(self.gauss_xi_s, size=None), np.random.choice(self.gauss_kaon, size=None)
            elif(decPr==11):
                if not bootstrap: return self.xi_mass,  self.kaon_s_mass
                else: return np.random.choice(self.gauss_xi, size=None), np.random.choice(self.gauss_kaon_s, size=None)
            elif(decPr==12):
                if not bootstrap: return self.xi_p_mass,  self.kaon_s_mass
                else: return np.random.choice(self.gauss_xi_p, size=None), np.random.choice(self.gauss_kaon_s, size=None)
            elif(decPr==13):
                if not bootstrap: return self.xi_s_mass,  self.kaon_s_mass
                else: return np.random.choice(self.gauss_xi_s, size=None), np.random.choice(self.gauss_kaon_s, size=None)
                # START OF SECOND DIAGRAM DECAYS
            elif(decPr==14):
                if not bootstrap: return self.proton_mass, self.B0_mass #  //N+D
                else: return np.random.choice(self.gauss_proton, size=None), np.random.choice(self.gauss_B0, size=None)
            # elif(decPr==15):
            #     if not bootstrap: return self.proton_mass, self.B_star_mass # //N+D*
            #     else: return np.random.choice(self.gauss_proton, size=None), np.random.choice(self.gauss_B_star, size=None)

        elif(baryons=='cascades_anti3'):
            if(decPr==1):
                if not bootstrap: return self.lambda_mass, self.kaon_mass
                else: return np.random.choice(self.gauss_lambda, size=None), np.random.choice(self.gauss_kaon, size=None)
            elif(decPr==2):
                if not bootstrap: return self.xi_mass,     self.pion_mass
                else: return np.random.choice(self.gauss_xi, size=None), np.random.choice(self.gauss_pion, size=None)
            elif(decPr==3):
                if not bootstrap: return self.xi_p_mass,   self.pion_mass
                else: return np.random.choice(self.gauss_xi_p, size=None), np.random.choice(self.gauss_pion, size=None)
            elif(decPr==4):
                if not bootstrap: return self.xi_s_mass,   self.pion_mass
                else: return np.random.choice(self.gauss_xi_s, size=None), np.random.choice(self.gauss_pion, size=None)
            elif(decPr==5):
                if not bootstrap: return self.sigma_mass,  self.kaon_mass
                else: return np.random.choice(self.gauss_sigma, size=None), np.random.choice(self.gauss_kaon, size=None)
            elif(decPr==6):
                if not bootstrap: return self.sigma_s_mass,self.kaon_mass
                else: return np.random.choice(self.gauss_sigma_s, size=None), np.random.choice(self.gauss_kaon, size=None)
            elif(decPr==7):
                if not bootstrap: return self.xi_mass,     self.eta_mass
                else: return np.random.choice(self.gauss_xi, size=None), np.random.choice(self.gauss_eta, size=None)
            elif(decPr==8):
                if not bootstrap: return self.lambda_mass,     self.kaon_s_mass
                else: return np.random.choice(self.gauss_lambda, size=None), np.random.choice(self.gauss_kaon_s, size=None)
            elif(decPr==9):
                if not bootstrap: return self.xi_mass,     self.meson_rho_mass
                else: return np.random.choice(self.gauss_xi, size=None), np.random.choice(self.gauss_meson_rho, size=None)
            elif(decPr==10):
                if not bootstrap: return self.xi_p_mass,     self.meson_rho_mass
                else: return np.random.choice(self.gauss_xi_p, size=None), np.random.choice(self.gauss_meson_rho, size=None)
            elif(decPr==11):
                if not bootstrap: return self.xi_s_mass,     self.meson_rho_mass
                else: return np.random.choice(self.gauss_xi_s, size=None), np.random.choice(self.gauss_meson_rho, size=None)
            elif(decPr==12):
                if not bootstrap: return self.sigma_mass,     self.kaon_s_mass
                else: return np.random.choice(self.gauss_sigma, size=None), np.random.choice(self.gauss_kaon_s, size=None)
            elif(decPr==13):
                if not bootstrap: return self.sigma_s_mass,     self.kaon_s_mass
                else: return np.random.choice(self.gauss_sigma_s, size=None), np.random.choice(self.gauss_kaon_s, size=None)
            elif(decPr==14):
                if not bootstrap: return self.xi_p_mass,     self.eta_mass
                else: return np.random.choice(self.gauss_xi_p, size=None), np.random.choice(self.gauss_eta, size=None)
            elif(decPr==15):
                if not bootstrap: return self.xi_s_mass,     self.eta_mass
                else: return np.random.choice(self.gauss_xi_s, size=None), np.random.choice(self.gauss_eta, size=None)
            elif(decPr==16):
                if not bootstrap: return self.xi_mass,     self.eta_p_mass
                else: return np.random.choice(self.gauss_xi, size=None), np.random.choice(self.gauss_eta_p, size=None)
            elif(decPr==17):
                if not bootstrap: return self.xi_p_mass,     self.eta_p_mass
                else: return np.random.choice(self.gauss_xi_p, size=None), np.random.choice(self.gauss_eta_p, size=None)
            elif(decPr==18):
                if not bootstrap: return self.xi_s_mass,     self.eta_p_mass
                else: return np.random.choice(self.gauss_xi_s, size=None), np.random.choice(self.gauss_eta_p, size=None)
            elif(decPr==19):
                if not bootstrap: return self.xi_mass,     self.small_omega_mass
                else: return np.random.choice(self.gauss_xi, size=None), np.random.choice(self.gauss_small_omega, size=None)
            elif(decPr==20):
                if not bootstrap: return self.xi_p_mass,     self.small_omega_mass
                else: return np.random.choice(self.gauss_xi_p, size=None), np.random.choice(self.gauss_small_omega, size=None)
            elif(decPr==21):
                if not bootstrap: return self.xi_s_mass,     self.small_omega_mass
                else: return np.random.choice(self.gauss_xi_s, size=None), np.random.choice(self.gauss_small_omega, size=None)
            elif(decPr==22):
                if not bootstrap: return self.xi_mass,     self.phi_mass
                else: return np.random.choice(self.gauss_xi, size=None), np.random.choice(self.gauss_phi, size=None)
            elif(decPr==23):
                if not bootstrap: return self.xi_p_mass,     self.phi_mass
                else: return np.random.choice(self.gauss_xi_p, size=None), np.random.choice(self.gauss_phi, size=None)
            elif(decPr==24):
                if not bootstrap: return self.xi_s_mass,     self.phi_mass
                else: return np.random.choice(self.gauss_xi_s, size=None), np.random.choice(self.gauss_phi, size=None)
                # START OF SECOND DIAGRAM DECAYS
            elif(decPr==25):
                if not bootstrap: return self.Lambda_8_mass, self.B0_mass # //Lambda_8_rho D
                else: return np.random.choice(self.gauss_Lambda_8, size=None), np.random.choice(self.gauss_B0, size=None)
            elif(decPr==26):
                if not bootstrap: return self.Lambda_8_mass, self.B_star_mass # //Lambda_8_rho D* 
                else: return np.random.choice(self.gauss_Lambda_8, size=None), np.random.choice(self.gauss_B_star, size=None)
            elif(decPr==27):
                if not bootstrap: return self.Sigma_8_mass, self.B0_mass #  //Sigma_8_rho D
                else: return np.random.choice(self.gauss_Sigma_8, size=None), np.random.choice(self.gauss_B0, size=None)
            elif(decPr==28):
                if not bootstrap: return self.Lambda_star_mass, self.B0_mass # //Lambda*_8_rho D
                else: return np.random.choice(self.gauss_Lambda_star, size=None), np.random.choice(self.gauss_B0, size=None)


    def fetch_decay_masses(self, bootstrap):
        self.pion_mass        = 0.13725
        self.kaon_mass        = 0.49564
        self.eta_mass         = 0.54786
        self.eta_p_mass       = 0.95778
        self.meson_rho_mass   = 0.77518
        self.kaon_s_mass      = 0.89555
        self.small_omega_mass = 0.78266
        self.phi_mass         = 1.01946        
        self.proton_mass      = 0.93891
        self.N_star_1520_mass = 1.51500
        self.N_star_1535_mass = 1.53000
        self.N_star_1680_mass = 1.68500
        self.N_star_1720_mass = 1.72000
        self.Delta_mass       = 1.23200
        self.Lambda_8_mass    = 1.11568
        self.Lambda_star_mass = 1.51900
        self.Xi_8_mass        = 1.31820
        self.Xi_10_mass       = 1.53370
        self.Sigma_8_mass     = 1.11932
        self.Sigma_10_mass    = 1.38460
        # Bottom hadrons
        self.lambda_mass      = 5.61960 # +- 0.0001
        self.xi_p_mass        = 5.93502 # +- 0.00005        
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
            self.gauss_pion        = np.random.normal(0.13725, 0.00295, 10000)
            self.gauss_kaon        = np.random.normal(0.49564, 0.00279, 10000)
            self.gauss_eta         = np.random.normal(0.54786, 0.00002, 10000)
            self.gauss_eta_p       = np.random.normal(0.95778, 0.00006, 10000)
            self.gauss_meson_rho   = np.random.normal(0.77518, 0.00045, 10000)
            self.gauss_kaon_s      = np.random.normal(0.89555, 0.00100, 10000)
            self.gauss_small_omega = np.random.normal(0.78266, 0.00002, 10000)
            self.gauss_phi         = np.random.normal(1.01946, 0.00002, 10000)
            self.gauss_proton      = np.random.normal(0.93891, 0.00091, 10000)
            self.gauss_N_star_1520 = np.random.normal(1.51500, 0.00500, 10000)
            self.gauss_N_star_1535 = np.random.normal(1.53000, 0.01500, 10000)
            self.gauss_N_star_1680 = np.random.normal(1.68500, 0.00500, 10000)
            self.gauss_N_star_1720 = np.random.normal(1.72000, 0.03500, 10000)
            self.gauss_Delta       = np.random.normal(1.23200, 0.00200, 10000)
            self.gauss_Lambda_8    = np.random.normal(1.11568, 0.00001, 10000)
            self.gauss_Lambda_star = np.random.normal(1.51900, 0.00010, 10000)
            self.gauss_Xi_8        = np.random.normal(1.31820, 0.00360, 10000)
            self.gauss_Xi_10       = np.random.normal(1.53370, 0.00250, 10000)
            self.gauss_Sigma_8     = np.random.normal(1.11932, 0.00340, 10000)
            self.gauss_Sigma_10    = np.random.normal(1.38460, 0.00460, 10000)
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

            
    def set_gamma_val(self, bootstrap):
        
        self.gamma_fit = 20.8 # this reproduces an specific decay width
        self.gauss_gamma = np.random.normal(20.8, np.power( 1.47**2 + 4.44**2, 0.5) , 10000)
        if(bootstrap):
            self.gauss_gamma = np.random.normal(20.8, np.power( 1.5**2 + 2.0**2, 0.5) , 10000) # error 2, is arbitrary

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
----------------------------------------------------------------
 Authors: A. Ramirez-Morales (andres.ramirez.morales@cern.ch)
          H. Garcia-Tecocoatzi
---------------------------------------------------------------
"""
import os
import numpy as np
import pandas as pd
# framework includes
from bottomfw.baryons import data_preparation as dp
from decays.decay_width import DecayWidths
from bottomfw.baryons import bottom_states as bs


class BottomDiquark:
    """
    Class that computes the mass and decay widths of three-quark heavy baryons.
    """    
    def __init__(self, baryons, params, sampled, corr_mat, asymmetric, batch_number=None, workpath='.'):
        self.params = params
        self.sampled = sampled
        self.corr_mat = corr_mat
        self.asymmetric = asymmetric
        self.m_batch_number = batch_number
        self.m_baryons = baryons
        self.m_workpath = workpath
    
    def model_mass(self, i, j, sampled=False):
        """
        Method that computes baryon masses according to model
        """
        if not sampled: # here the mass is calculated using the average value of parameters
            return  self.md1*self.is_omega_param[i] + self.md2*self.is_cascade_param[i] + self.md3*self.is_sigma_param[i] + self.mb + \
                self.omega_harmonic(self.Kp, self.md1, self.md2, self.md3, self.mb,
                                    self.is_omega_param[i], self.is_cascade_param[i], self.is_sigma_param[i], self.v_param[i] ) + \
                self.A*self.w_param[i] + self.B*self.x_param[i] + self.E*self.y_param[i] + self.G*self.z_param[i]
        else: # here the mass is calculated using every single bootstrap simulated parameter
            return  self.sampled_md1[j]*self.is_omega_param[i] + self.sampled_md2[j]*self.is_cascade_param[i] + self.sampled_md3[j]*self.is_sigma_param[i] + \
                self.sampled_mb[j] + \
                self.omega_harmonic(self.sampled_k[j], self.sampled_md1[j], self.sampled_md2[j], self.sampled_md3[j], self.sampled_mb[j],
                                    self.is_omega_param[i], self.is_cascade_param[i], self.is_sigma_param[i], self.v_param[i] ) + \
                self.sampled_a[j]*self.w_param[i] + self.sampled_b[j]*self.x_param[i] + \
                self.sampled_e[j]*self.y_param[i] + self.sampled_g[j]*self.z_param[i]

    
    def omega_harmonic(self, kp_input, md1_input, md2_input, md3_input, mb_input, is_omega_input, is_cascade_input, is_sigma_input, v_input):
        """
        Method that calculates the h.o. frequencies(omegas)
        """
        if(v_input==0.): return 0 # no contribution from harmonic oscillator, i.e. base states
        m_d = is_omega_input*md1_input + is_cascade_input*md2_input + is_sigma_input*md3_input
        m_u = m_d*mb_input/(m_d + mb_input)
        den_m_u = np.sqrt(1./m_u)
        omega = v_input * den_m_u  * kp_input
        return omega
                
    def paper_results_predictions(self, bootstrap=False, bootstrap_width=False, prev_params=False, decay_width=False):
        """
        Method to call the calculations and save them (paper latex tables, here we include all the states)
        The parameters are redifined. Will clean the input param
        """
        self.reload_quantum_param(self.m_baryons) #available: baryons=omegas,cascades,sigmas,lambdas,cascades_anti3
        
        if prev_params: # use the parameters of the previous paper
            self.previous_parameters()
            self.previous_parameters_uncertainty(N_boots=len(self.sampled_k)) # set params. with previous-paper gauss shape (arbitrary error)

        # compute masses/decays and save them in csv files
        self.compute_save_predictions(self.m_baryons, bootstrap=bootstrap, decay_width=decay_width, bootstrap_width=bootstrap_width)           
                    
        
    def compute_save_predictions(self, baryons='', bootstrap=False, decay_width=False, bootstrap_width=False):
        """
        Methods that calculates the masses and decay widths with the bootstraped fitted parameters
        """
        masses_csv = []
        decays_csv = []        
        for i in range(len(self.v_param)): # states loop
            dummy,dummy_decay = ([]),([])
            if bootstrap:
                for j in range(len(self.sampled_k)): # sampled data loop (e.g. 10^5)
                    mass = self.model_mass(i, j, sampled=True)
                    dummy = np.append(dummy, mass)
                    if decay_width and bootstrap_width: # and self.L_tot[i]==1  # decayWidth calculation, DecayWidths class
                        decay_value = self.baryon_decay.total_decay_width(baryons, self.sampled_k[j], mass,
                                                                          self.S_tot[i], self.L_tot[i], self.J_tot[i],
                                                                          self.SL[i], self.ModEx[i], bootstrap=bootstrap_width,
                                                                          m1=self.sampled_m1[j], m2=self.sampled_m2[j], m3=self.sampled_m3[j])
                        dummy_decay = np.append(dummy_decay, decay_value)
            else: # no boostrap, only one prediction using the average of the fitted parameters
                mass = self.model_mass(i, 0, sampled=False)
                dummy = np.append(dummy, round(mass))
                if decay_width:
                    decay = self.baryon_decay.total_decay_width(baryons, self.Kp, mass,
                                                                self.S_tot[i], self.L_tot[i], self.J_tot[i], self.SL[i],
                                                                self.ModEx[i], bootstrap=False, m1=self.m1, m2=self.m1, m3=self.m3)
                    dummy_decay = np.append(dummy_decay, decay)
                                    
            # prepare data frame to save as csv
            masses_csv.append(pd.Series(dummy))
            if decay_width:
                decays_csv.append(pd.Series(dummy_decay))
                
                # for individual decay channels table (this is really tricky!!!)
                # TODO: save to a csv file
                channels_indi = np.array(self.baryon_decay.channel_widths_vector)
                self.baryon_decay.channel_widths_vector = [] # clean,get ready for next mass[i]
                # last tine of states loop
                    
        df_masses = None
        df_decays = None        
        if len(masses_csv) != 0:
            keys_names = [str(i)+"_state" for i in range(len(masses_csv))]
            df_masses = pd.concat(masses_csv, axis=1, keys=keys_names)

        if len(decays_csv) != 0:
            keys_names = [str(i)+"_state" for i in range(len(decays_csv))]
            df_decays = pd.concat(decays_csv, axis=1, keys=keys_names)            
        # save results in cvs files
        if self.m_batch_number is None:
            if df_masses is not None:
                if not os.path.exists(self.m_workpath+"/tables/"):
                    os.mkdir(self.m_workpath+"/tables/")
                df_masses.to_csv(self.m_workpath+"/tables/masses_states_diquark_"+self.m_baryons+".csv", index=False)
            if df_decays is not None:
                if not os.path.exists(self.m_workpath+"/tables/"):
                    os.mkdir(self.m_workpath+"/tables/")
                df_decays.to_csv(self.m_workpath+"/tables/decays_states_diquark_"+self.m_baryons+".csv", index=False)                
        else:
            if df_masses is not None:
                if not os.path.exists(self.m_workpath+"/batch_results_diquark/"+self.m_baryons+"/mass_states/"):
                    os.mkdir(self.m_workpath+"/batch_results_diquark/"+self.m_baryons+"/mass_states/")
                df_masses.to_csv(self.m_workpath+"/batch_results_diquark/"+self.m_baryons+"/mass_states/"+str(self.m_batch_number)+".csv", index=False)
            if df_decays is not None:
                if not os.path.exists(self.m_workpath+"/batch_results_diquark/"+self.m_baryons+"/decay_states/"):
                    os.mkdir(self.m_workpath+"/batch_results_diquark/"+self.m_baryons+"/decay_states/")                
                df_decays.to_csv(self.m_workpath+"/batch_results_diquark/"+self.m_baryons+"/decay_states/"+str(self.m_batch_number)+".csv", index=False)
                
    def fetch_values(self):
        """
        Fetch the values from the input dictionaries    
        """
        # sampled        
        self.sampled_md1 = self.sampled['sampled_md1']
        self.sampled_md2 = self.sampled['sampled_md2']
        self.sampled_md3 = self.sampled['sampled_md3']
        self.sampled_mb  = self.sampled['sampled_mb']

        self.sampled_k = self.sampled['sampled_k']
        self.sampled_a = self.sampled['sampled_a']
        self.sampled_b = self.sampled['sampled_b']
        self.sampled_e = self.sampled['sampled_e']
        self.sampled_g = self.sampled['sampled_g']
        
        self.md1,self.delta_md1= np.mean(self.sampled['sampled_md1']), np.std(self.sampled['sampled_md1'],ddof=1)
        self.md2,self.delta_md2= np.mean(self.sampled['sampled_md2']), np.std(self.sampled['sampled_md2'],ddof=1)
        self.md3,self.delta_md3= np.mean(self.sampled['sampled_md3']), np.std(self.sampled['sampled_md3'],ddof=1)
        self.mb,self.delta_mb= np.mean(self.sampled['sampled_mb']), np.std(self.sampled['sampled_mb'],ddof=1)

        self.Kp,self.delta_Kp= np.mean(self.sampled['sampled_k']), np.std(self.sampled['sampled_k'],ddof=1)
        self.A, self.delta_A = np.mean(self.sampled['sampled_a']), np.std(self.sampled['sampled_a'],ddof=1)
        self.B, self.delta_B = np.mean(self.sampled['sampled_b']), np.std(self.sampled['sampled_b'],ddof=1)
        self.E, self.delta_E = np.mean(self.sampled['sampled_e']), np.std(self.sampled['sampled_e'],ddof=1)
        self.G, self.delta_G = np.mean(self.sampled['sampled_g']), np.std(self.sampled['sampled_g'],ddof=1)

        # parameters
        self.is_omega_param  = self.params['is_omega']
        self.is_cascade_param= self.params['is_cascade']
        self.is_sigma_param  = self.params['is_sigma']
                
        self.v_param = self.params['V']
        self.w_param = self.params['W']
        self.x_param = self.params['X']
        self.y_param = self.params['Y']
        self.z_param = self.params['Z']
        
        # correlation matrix
        self.rho_md2md1 = np.mean(self.corr_mat['rho_md2md1'])
        self.rho_md3md1 = np.mean(self.corr_mat['rho_md3md1'])
        self.rho_mbmd1  = np.mean(self.corr_mat['rho_mbmd1'])
        self.rho_kmd1   = np.mean(self.corr_mat['rho_kmd1'])
        self.rho_amd1   = np.mean(self.corr_mat['rho_amd1'])
        self.rho_bmd1   = np.mean(self.corr_mat['rho_bmd1'])
        self.rho_emd1   = np.mean(self.corr_mat['rho_emd1'])
        self.rho_gmd1   = np.mean(self.corr_mat['rho_gmd1'])
        self.rho_md3md2 = np.mean(self.corr_mat['rho_md3md2'])
        self.rho_mbmd2  = np.mean(self.corr_mat['rho_mbmd2'])
        self.rho_kmd2   = np.mean(self.corr_mat['rho_kmd2'])
        self.rho_amd2   = np.mean(self.corr_mat['rho_amd2'])
        self.rho_bmd2   = np.mean(self.corr_mat['rho_bmd2'])
        self.rho_emd2   = np.mean(self.corr_mat['rho_emd2'])
        self.rho_gmd2   = np.mean(self.corr_mat['rho_gmd2'])
        self.rho_mbmd3  = np.mean(self.corr_mat['rho_mbmd3'])
        self.rho_kmd3   = np.mean(self.corr_mat['rho_kmd3'])
        self.rho_amd3   = np.mean(self.corr_mat['rho_amd3'])
        self.rho_bmd3   = np.mean(self.corr_mat['rho_bmd3'])
        self.rho_emd3   = np.mean(self.corr_mat['rho_emd3'])
        self.rho_gmd3   = np.mean(self.corr_mat['rho_gmd3'])
        self.rho_kmb    = np.mean(self.corr_mat['rho_kmb'])
        self.rho_amb    = np.mean(self.corr_mat['rho_amb'])
        self.rho_bmb    = np.mean(self.corr_mat['rho_bmb'])
        self.rho_emb    = np.mean(self.corr_mat['rho_emb'])
        self.rho_gmb    = np.mean(self.corr_mat['rho_gmb'])
        self.rho_ak     = np.mean(self.corr_mat['rho_ak'])
        self.rho_bk     = np.mean(self.corr_mat['rho_bk'])
        self.rho_ek     = np.mean(self.corr_mat['rho_ek'])
        self.rho_gk     = np.mean(self.corr_mat['rho_gk'])
        self.rho_ba     = np.mean(self.corr_mat['rho_ba'])
        self.rho_ea     = np.mean(self.corr_mat['rho_ea'])
        self.rho_ga     = np.mean(self.corr_mat['rho_ga'])
        self.rho_eb     = np.mean(self.corr_mat['rho_eb'])
        self.rho_gb     = np.mean(self.corr_mat['rho_gb'])
        self.rho_ge     = np.mean(self.corr_mat['rho_ge'])

              
    def reload_quantum_param(self,baryons):
        """
        Load the quantum parameters that multiply the A,B,E,G in Hamiltonian
        """
        # get the original quantum numbers and experimental data
        state,sum_mass,J_tot,S_tot,L_tot,I_tot,SU_tot,HO_n,SL,ModEx = bs.states_mass_diquark(baryons)
        # get the hamiltonian factors based on the quatum numbers
        param_v, param_q1, param_q2, param_q3, param_is_rho, param_is_lam,\
            param_is_omega, param_is_cascade, param_is_sigma,\
            param_w,param_x,param_y,param_z  = \
                dp.hamiltonian_quantum_factors_extended(state,J_tot,
                                                        S_tot,L_tot,I_tot,SU_tot,HO_n,ModEx)
        # parameters re-assigment
        self.v_param = param_v
        self.w_param = param_w
        self.x_param = param_x
        self.y_param = param_y
        self.z_param = param_z
        # self.sum_mass 
        self.q1_param= param_q1
        self.q2_param= param_q2
        self.q3_param= param_q3
        self.is_rho_param    = param_is_rho
        self.is_lam_param    = param_is_lam
        self.is_omega_param  = param_is_omega
        self.is_cascade_param= param_is_cascade
        self.is_sigma_param  = param_is_sigma        
        # quantum number assigment
        self.state  = state
        self.J_tot  = J_tot
        self.S_tot  = S_tot
        self.L_tot  = L_tot
        self.I_tot  = I_tot
        self.SU_tot = SU_tot
        self.SL     = SL
        self.HO_n   = HO_n 
        self.ModEx  = ModEx

    def previous_parameters(self):
        """
        Set parameters used in previous publication
        """
        self.Kp =  4847.67
        self.A  =  6.74
        self.B  =  5.15
        self.E  =  26
        self.G  =  70.91
        # quark masses
        self.m1 = 4920
        self.m2 = 450
        self.m3 = 295
        
    def previous_parameters_uncertainty(self, N_boots):
        """
        set the values of the previous paper as gaussian (this is for testing decay widths)
        """
        self.sampled_k = np.random.normal(4847.67, self.delta_Kp, N_boots)
        self.sampled_a = np.random.normal(6.74,   0.37, N_boots)
        self.sampled_b = np.random.normal(5.15,   0.31, N_boots)
        self.sampled_e = np.random.normal(26,   0.23, N_boots)
        self.sampled_g = np.random.normal(70.91,   0.49, N_boots)

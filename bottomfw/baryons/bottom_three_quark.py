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
from decays.electro_width import ElectroWidths
from bottomfw.baryons import bottom_states as bs


class BottomThreeQuark:
    """
    Class that computes the mass and decay widths of three-quark heavy baryons.
    """
    def __init__(self, baryons, params, sampled, corr_mat, asymmetric, decay_width=False, bootstrap_width=False, decay_width_em=False, bootstrap_width_em=False, batch_number=None, workpath='.'):
        self.params = params
        self.sampled = sampled
        self.corr_mat = corr_mat
        self.asymmetric = asymmetric
        self.m_decay_width = decay_width
        self.m_decay_width_em = decay_width_em # check
        self.m_bootstrap_width_em = bootstrap_width_em # check
        self.m_batch_number = batch_number
        self.m_baryons = baryons
        self.m_workpath = workpath
        if decay_width:
            self.baryon_decay = DecayWidths(bootstrap_width, baryons, workpath=workpath)
        if decay_width_em:
            self.electro_decay = ElectroWidths(bootstrap_width_em, baryons, workpath=workpath)
    
    def model_mass(self, i, j, sampled=False):
        """
        Method that computes baryon masses according to model
        """
        if not sampled: # here the mass is calculated using the average value of parameters
            return  self.m1*self.q1_param[i] + self.m2*self.q2_param[i] + self.m3*self.q3_param[i] + \
                self.omega_harmonic(self.Kp, self.m1, self.m2, self.m3,
                                    self.is_rho_param[i], self.is_lam_param[i],
                                    self.is_omega_param[i], self.is_cascade_param[i], self.is_sigma_param[i], self.v_param[i] ) + \
                self.A*self.w_param[i]  + self.B*self.x_param[i] + self.E*self.y_param[i] + self.G*self.z_param[i]
        else: # here the mass is calculated using every single bootstrap simulated parameter
            return  self.sampled_m1[j]*self.q1_param[i] + self.sampled_m2[j]*self.q2_param[i] + self.sampled_m3[j]*self.q3_param[i] + \
                self.omega_harmonic(self.sampled_k[j], self.sampled_m1[j], self.sampled_m2[j], self.sampled_m3[j],
                                    self.is_rho_param[i], self.is_lam_param[i],
                                    self.is_omega_param[i], self.is_cascade_param[i], self.is_sigma_param[i], self.v_param[i] ) + \
                self.sampled_a[j]*self.w_param[i] + self.sampled_b[j]*self.x_param[i] + \
                self.sampled_e[j]*self.y_param[i] + self.sampled_g[j]*self.z_param[i]
        
    def get_omega_harmonic(self, i, j, sampled=False):
        """
        Method that fetches the h.o. frequencies(omegas)
        """
        if not sampled:
            return self.omega_harmonic(self.Kp, self.m1, self.m2, self.m3,
                                       self.is_rho_param[i], self.is_lam_param[i],
                                       self.is_omega_param[i], self.is_cascade_param[i], self.is_sigma_param[i], self.v_param[i])
        else:
            return self.omega_harmonic(self.sampled_k[j], self.sampled_m1[j], self.sampled_m2[j], self.sampled_m3[j],
                                       self.is_rho_param[i], self.is_lam_param[i],
                                       self.is_omega_param[i], self.is_cascade_param[i], self.is_sigma_param[i], self.v_param[i])

    @staticmethod
    def omega_harmonic(kp_input, m1_input, m2_input, m3_input, is_rho_input, is_lam_input, is_omega_input, is_cascade_input, is_sigma_input, v_input):
        """
        Method that calculates the h.o. frequencies(omegas)
        """
        if(v_input==0.): return 0 # no contribution from harmonic oscillator, i.e. base states

        m_rho = is_rho_input*(is_omega_input*m2_input + is_cascade_input*((m2_input+m3_input)/2) + is_sigma_input*m3_input )
        m_lam = is_lam_input*(is_omega_input*((3*m2_input*m1_input)/(2*m2_input+m1_input)) + \
                           is_cascade_input*((1.5*(m2_input+m3_input)*m1_input)/(m1_input+m2_input+m3_input)) + is_sigma_input*((3*m3_input*m1_input)/(2*m3_input+m1_input)))
        if(m_lam == 0):
            den_lam = 0            
        else:
            den_lam = np.sqrt(1./m_lam)

        if(m_rho == 0):
            den_rho = 0
        else:            
            den_rho = np.sqrt(1./m_rho)
        
        omega = v_input * (den_lam + den_rho) * kp_input
        return omega
                
    def paper_results_predictions(self, bootstrap=False, bootstrap_width=False, prev_params=False):
        """
        Method to call the calculations and save them (paper latex tables, here we include all the states)
        The parameters are redifined. Will clean the input param
        """
        self.reload_quantum_param(self.m_baryons) #available: baryons=omegas,cascades,sigmas,lambdas,cascades_anti3
        if prev_params: # use the parameters of the previous paper
            self.previous_parameters()
            self.previous_parameters_uncertainty(N_boots=len(self.sampled_k)) # set params. with previous-paper gauss shape (arbitrary error)
        # compute masses/decays and save them in csv files
        self.compute_save_predictions(self.m_baryons, bootstrap=bootstrap, decay_width=self.m_decay_width, bootstrap_width=bootstrap_width, decay_width_em=self.m_decay_width_em, bootstrap_width_em=self.m_bootstrap_width_em)

    def compute_save_predictions(self, baryons='', bootstrap=False, decay_width=False, bootstrap_width=False, decay_width_em=False, bootstrap_width_em=False):
        """
        Methods that calculates the masses and decay widths with the bootstraped fitted parameters
        """
        masses_csv = []
        omegas_csv = []
        decays_csv = []
        decays_csv_em = []
        for i in range(len(self.v_param)): # states loop
            dummy,dummy_decay,dummy_omega,dummy_decay_em = ([]),([]),([]),([])
            decays_indi_csv = []
            decays_electro_indi_csv = []
            if bootstrap:
                for j in range(len(self.sampled_k)): # sampled data loop (e.g. 10^4)
                    mass = self.model_mass(i, j, sampled=True)
                    mass_avg = self.model_mass(i, 0, sampled=False)
                    omega_ho = self.get_omega_harmonic(i, j, sampled=False)
                    dummy = np.append(dummy, mass)
                    dummy_omega = np.append(dummy_omega, omega_ho)
                    if decay_width and bootstrap_width: # and self.L_tot[i]==1  # decayWidth calculation, DecayWidths class
                        self.baryon_decay.load_average_mass(mass_avg) # load central value of massA to enforce energy conservation
                        decay_value = self.baryon_decay.total_decay_width(baryons, self.sampled_k[j], mass,
                                                                          self.S_tot[i], self.L_tot[i], self.J_tot[i], self.SL[i],
                                                                          self.ModEx[i], bootstrap=bootstrap_width, m1=self.sampled_m1[j],
                                                                          m2=self.sampled_m2[j], m3=self.sampled_m3[j])            
                        dummy_decay = np.append(dummy_decay, decay_value) # total decay
                        decays_indi_csv.append(self.baryon_decay.channel_widths_vector[0]) # individual channels decays
                        self.baryon_decay.channel_widths_vector=[]

                    if decay_width_em and bootstrap_width_em and self.L_tot[i]<=1 and (self.ModEx[i]=="grd" or self.ModEx[i]=="lam" or self.ModEx[i]=="rho"): # bootstrap mass but not bootstrap widhts
                        electro_decay = self.electro_decay.total_decay_width(baryons, self.sampled_k[j], mass,
                                                                             self.S_tot[i], self.J_tot[i], self.L_tot[i], self.SL[i],
                                                                             self.ModEx[i], bootstrap=bootstrap_width_em, m1=self.sampled_m1[j], m2=self.sampled_m2[j], m3=self.sampled_m3[j])
                        dummy_decay_em = np.append(dummy_decay_em, electro_decay) # total electro decay
                        decays_electro_indi_csv.append(self.electro_decay.channel_widths_vector[0]) # individual channel electro decays
                        self.electro_decay.channel_widths_vector=[] # clean decay object for next iteration

                if decay_width and not bootstrap_width: # bootstrap mass but not bootstrap widhts
                    mass_single = self.model_mass(i, 0, sampled=False)
                    decay = self.baryon_decay.total_decay_width(baryons, self.Kp, mass_single,
                                                                self.S_tot[i], self.L_tot[i], self.J_tot[i], self.SL[i],
                                                                self.ModEx[i], bootstrap=False, m1=self.m1, m2=self.m2, m3=self.m3)
                    dummy_decay = np.append(dummy_decay, decay) # total decay
                    decays_indi_csv.append(self.baryon_decay.channel_widths_vector[0]) # individual channel decays
                    self.baryon_decay.channel_widths_vector=[] # clean decay object for next iteration

                if decay_width_em and not bootstrap_width_em and self.L_tot[i]<=1 and (self.ModEx[i]=="grd" or self.ModEx[i]=="lam" or self.ModEx[i]=="rho"): # bootstrap mass but not bootstrap widhts
                    mass_single = self.model_mass(i, 0, sampled=False)
                    electro_decay = self.electro_decay.total_decay_width(baryons, self.Kp, mass_single,
                                                                         self.S_tot[i], self.J_tot[i], self.L_tot[i], self.SL[i],
                                                                         self.ModEx[i], bootstrap=False, m1=self.m1, m2=self.m2, m3=self.m3)
                    dummy_decay_em = np.append(dummy_decay_em, electro_decay) # total electro decay
                    decays_electro_indi_csv.append(self.electro_decay.channel_widths_vector[0]) # individual channel electro decays
                    self.electro_decay.channel_widths_vector=[] # clean decay object for next iteration

            else: # no bootstrap at all, only one prediction using the average of the fitted parameters
                mass = self.model_mass(i, 0, sampled=False)
                dummy = np.append(dummy, mass)
                if decay_width: # strong decays
                    strong_decay = self.baryon_decay.total_decay_width(baryons, self.Kp, mass,
                                                                       self.S_tot[i], self.L_tot[i], self.J_tot[i], self.SL[i],
                                                                       self.ModEx[i], bootstrap=False, m1=self.m1, m2=self.m2, m3=self.m3)
                    dummy_decay    = np.append(dummy_decay, strong_decay) # total decay
                    decays_indi_csv.append(pd.DataFrame(self.baryon_decay.channel_widths_vector)) # individual channel strong decays
                    self.baryon_decay.channel_widths_vector=[]  # clean decay object for next iteration
                if decay_width_em  and self.L_tot[i]<=1  and (self.ModEx[i]=="grd" or self.ModEx[i]=="lam" or self.ModEx[i]=="rho"): # electromagnetic decays only up to P-wave
                    electro_decay = self.electro_decay.total_decay_width(baryons, self.Kp, mass,
                                                                         self.S_tot[i], self.J_tot[i], self.L_tot[i], self.SL[i],
                                                                         self.ModEx[i], bootstrap=False, m1=self.m1, m2=self.m2, m3=self.m3)
                    dummy_decay_em = np.append(dummy_decay_em, electro_decay) # total electro decay
                    decays_electro_indi_csv.append(pd.DataFrame(self.electro_decay.channel_widths_vector)) # individual channel electro decays
                    self.electro_decay.channel_widths_vector=[] # clean decay object for next iteration

            # prepare data to a csv file
            masses_csv.append(pd.Series(dummy))
            omegas_csv.append(pd.Series(dummy_omega))
            if decay_width: # strong decays
                decays_csv.append(pd.Series(dummy_decay))
                df_decays_indi = None
                if len(decays_indi_csv) != 0:
                    decay_columns = []
                    for k in range(len(decays_indi_csv[0])):
                        dummy_column = []
                        for j in range(len(decays_indi_csv)):
                            dummy_column.append(decays_indi_csv[j][k])
                        decay_columns.append(dummy_column)
                    
                    channel_names = [str(k)+"_channel" for k in range(len(decays_indi_csv[0]))]
                    df_decays_indi = pd.DataFrame({channel_names[0]: decay_columns[0]})
                    for k in range(len(decays_indi_csv[0])-1):            
                        df_decays_indi[channel_names[k+1]]=decay_columns[k+1]                    
                    if self.m_batch_number is None:
                        if df_decays_indi is not None:
                            if not os.path.exists(self.m_workpath+"/tables/decays_indi/"):
                                os.makedirs(self.m_workpath+"/tables/decays_indi/")                            
                            df_decays_indi.to_csv(self.m_workpath+"/tables/decays_indi/decays_state_"+str(i)+"_"+self.m_baryons+".csv", index=False)
                    else:
                        if df_decays_indi is not None:  # save results for batch a given batch job
                            dec_dir = self.m_workpath+"/batch_results/"+self.m_baryons+"/decays_indi/state_"+str(i)
                            if not os.path.exists(dec_dir):
                                os.makedirs(dec_dir)
                            df_decays_indi.to_csv(dec_dir+"/"+str(self.m_batch_number)+".csv", index=False)                         
                # last line of states loop
            if decay_width_em and self.L_tot[i]<=1 and (self.ModEx[i]=="grd" or self.ModEx[i]=="lam" or self.ModEx[i]=="rho"): # electromagnetic decays
                decays_csv_em.append(pd.Series(dummy_decay_em))
                df_decays_indi_em = None
                if len(decays_electro_indi_csv) != 0:
                    decay_columns_em = []
                    for k in range(len(decays_electro_indi_csv[0])):
                        dummy_column_em = []
                        for j in range(len(decays_electro_indi_csv)):
                            dummy_column_em.append(decays_electro_indi_csv[j][k])
                        decay_columns_em.append(dummy_column_em)
                    
                    channel_names_em = [str(k)+"_channel" for k in range(len(decays_electro_indi_csv[0]))]
                    df_decays_indi_em = pd.DataFrame({channel_names_em[0]: decay_columns_em[0]})
                    for k in range(len(decays_electro_indi_csv[0])-1):
                        df_decays_indi_em[channel_names_em[k+1]]=decay_columns_em[k+1]
                    if self.m_batch_number is None:
                        if df_decays_indi_em is not None:
                            if not os.path.exists(self.m_workpath+"/tables/decays_indi_em/"):
                                os.makedirs(self.m_workpath+"/tables/decays_indi_em/")                            
                            df_decays_indi_em.to_csv(self.m_workpath+"/tables/decays_indi_em/decays_state_"+str(i)+"_"+self.m_baryons+".csv", index=False)
                    else:
                        if df_decays_indi_em is not None:  # save results for batch a given batch job
                            dec_dir_em = self.m_workpath+"/batch_results/"+self.m_baryons+"/decays_indi_em/state_"+str(i)
                            if not os.path.exists(dec_dir_em):
                                os.makedirs(dec_dir_em)
                            df_decays_indi_em.to_csv(dec_dir_em+"/"+str(self.m_batch_number)+".csv", index=False)                                                   
                # last line of states loop
        df_masses = None
        df_omegas = None
        df_decays = None
        df_electro= None
        if len(masses_csv) != 0:
            keys_names = [str(i)+"_state" for i in range(len(masses_csv))]
            df_masses = pd.concat(masses_csv, axis=1, keys=keys_names)

        if len(omegas_csv) != 0:
            keys_names = [str(i)+"_state" for i in range(len(omegas_csv))]
            df_omegas = pd.concat(omegas_csv, axis=1, keys=keys_names)

        if len(decays_csv) != 0:
            keys_names = [str(i)+"_state" for i in range(len(decays_csv))]
            df_decays = pd.concat(decays_csv, axis=1, keys=keys_names)

        if len(decays_csv_em) != 0:
            keys_names = [str(i)+"_state" for i in range(len(decays_csv_em))]
            df_electro = pd.concat(decays_csv_em, axis=1, keys=keys_names)

        # save in csv files
        if self.m_batch_number is None:
            if df_masses is not None:
                if not os.path.exists(self.m_workpath+"/tables/"):
                    os.makedirs(self.m_workpath+"/tables/")
                df_masses.to_csv(self.m_workpath+"/tables/masses_states_"+self.m_baryons+".csv", index=False)
            if df_omegas is not None:
                if not os.path.exists(self.m_workpath+"/tables/"):
                    os.makedirs(self.m_workpath+"/tables/")
                df_omegas.to_csv(self.m_workpath+"/tables/harmonic_states_"+self.m_baryons+".csv", index=False)
            if df_decays is not None:
                if not os.path.exists(self.m_workpath+"/tables/"):
                    os.makedirs(self.m_workpath+"/tables/")
                df_decays.to_csv(self.m_workpath+"/tables/decays_states_"+self.m_baryons+".csv", index=False)
            if df_electro is not None:
                if not os.path.exists(self.m_workpath+"/tables/"):
                    os.makedirs(self.m_workpath+"/tables/")
                df_electro.to_csv(self.m_workpath+"/tables/decays_electro_states_"+self.m_baryons+".csv", index=False)
        else:
            if df_masses is not None:  # save results for batch a given batch job
                if not os.path.exists(self.m_workpath+"/batch_results/"+self.m_baryons+"/mass_states/"):
                    os.makedirs(self.m_workpath+"/batch_results/"+self.m_baryons+"/mass_states/")
                df_masses.to_csv(self.m_workpath+"/batch_results/"+self.m_baryons+"/mass_states/"+str(self.m_batch_number)+".csv", index=False)
            if df_decays is not None:
                if not os.path.exists(self.m_workpath+"/batch_results/"+self.m_baryons+"/decay_states/"):
                    os.makedirs(self.m_workpath+"/batch_results/"+self.m_baryons+"/decay_states/")
                df_decays.to_csv(self.m_workpath+"/batch_results/"+self.m_baryons+"/decay_states/"+str(self.m_batch_number)+".csv", index=False)
            if df_electro is not None:
                if not os.path.exists(self.m_workpath+"/batch_results/"+self.m_baryons+"/decay_states_electro/"):
                    os.makedirs(self.m_workpath+"/batch_results/"+self.m_baryons+"/decay_states_electro/")
                df_electro.to_csv(self.m_workpath+"/batch_results/"+self.m_baryons+"/decay_states_electro/"+str(self.m_batch_number)+".csv", index=False)
    
    def fetch_values(self):
        """
        Fetch the values from the input dictionaries    
        """
        # sampled
        self.sampled_m1 = self.sampled['sampled_m1']        
        self.sampled_m2 = self.sampled['sampled_m2']
        self.sampled_m3 = self.sampled['sampled_m3']
        
        self.sampled_k = self.sampled['sampled_k']
        self.sampled_a = self.sampled['sampled_a']
        self.sampled_b = self.sampled['sampled_b']
        self.sampled_e = self.sampled['sampled_e']
        self.sampled_g = self.sampled['sampled_g']
        
        self.m1,self.delta_m1= np.mean(self.sampled['sampled_m1']), np.std(self.sampled['sampled_m1'],ddof=1)
        self.m2,self.delta_m2= np.mean(self.sampled['sampled_m2']), np.std(self.sampled['sampled_m2'],ddof=1)
        self.m3,self.delta_m3= np.mean(self.sampled['sampled_m3']), np.std(self.sampled['sampled_m3'],ddof=1)

        self.Kp,self.delta_Kp= np.mean(self.sampled['sampled_k']), np.std(self.sampled['sampled_k'],ddof=1)
        self.A, self.delta_A = np.mean(self.sampled['sampled_a']), np.std(self.sampled['sampled_a'],ddof=1)
        self.B, self.delta_B = np.mean(self.sampled['sampled_b']), np.std(self.sampled['sampled_b'],ddof=1)
        self.E, self.delta_E = np.mean(self.sampled['sampled_e']), np.std(self.sampled['sampled_e'],ddof=1)
        self.G, self.delta_G = np.mean(self.sampled['sampled_g']), np.std(self.sampled['sampled_g'],ddof=1)

        # parameters
        self.q1_param= self.params['q1']
        self.q2_param= self.params['q2']
        self.q3_param= self.params['q3']

        self.is_rho_param    = self.params['is_rho']
        self.is_lam_param    = self.params['is_lam']
        self.is_omega_param  = self.params['is_omega']
        self.is_cascade_param= self.params['is_cascade']
        self.is_sigma_param  = self.params['is_sigma']
                
        self.v_param = self.params['V']
        self.w_param = self.params['W']
        self.x_param = self.params['X']
        self.y_param = self.params['Y']
        self.z_param = self.params['Z']
        # correlation matrix
        self.rho_ak = np.mean(self.corr_mat['rho_ak'])
        self.rho_bk = np.mean(self.corr_mat['rho_bk'])
        self.rho_ba = np.mean(self.corr_mat['rho_ba'])
        self.rho_ek = np.mean(self.corr_mat['rho_ek'])
        self.rho_ea = np.mean(self.corr_mat['rho_ea'])
        self.rho_eb = np.mean(self.corr_mat['rho_eb'])
        self.rho_gk = np.mean(self.corr_mat['rho_gk'])
        self.rho_ga = np.mean(self.corr_mat['rho_ga'])
        self.rho_gb = np.mean(self.corr_mat['rho_gb'])
        self.rho_ge = np.mean(self.corr_mat['rho_ge'])
        # extended
        self.rho_m2m1 = np.mean(self.corr_mat['rho_m2m1'])
        self.rho_m3m1 = np.mean(self.corr_mat['rho_m3m1'])
        self.rho_km1  = np.mean(self.corr_mat['rho_km1'])
        self.rho_am1  = np.mean(self.corr_mat['rho_am1'])
        self.rho_bm1  = np.mean(self.corr_mat['rho_bm1'])
        self.rho_em1  = np.mean(self.corr_mat['rho_em1'])
        self.rho_gm1  = np.mean(self.corr_mat['rho_gm1'])
        self.rho_m3m2 = np.mean(self.corr_mat['rho_m3m2'])
        self.rho_km2  = np.mean(self.corr_mat['rho_km2'])
        self.rho_am2  = np.mean(self.corr_mat['rho_am2'])
        self.rho_bm2  = np.mean(self.corr_mat['rho_bm2'])
        self.rho_em2  = np.mean(self.corr_mat['rho_em2'])
        self.rho_gm2  = np.mean(self.corr_mat['rho_gm2'])
        self.rho_km3  = np.mean(self.corr_mat['rho_km3'])
        self.rho_am3  = np.mean(self.corr_mat['rho_am3'])
        self.rho_bm3  = np.mean(self.corr_mat['rho_bm3'])
        self.rho_em3  = np.mean(self.corr_mat['rho_em3'])
        self.rho_gm3  = np.mean(self.corr_mat['rho_gm3'])        
        self.rho_ak   = np.mean(self.corr_mat['rho_ak'])
        self.rho_bk   = np.mean(self.corr_mat['rho_bk'])
        self.rho_ek   = np.mean(self.corr_mat['rho_ek'])
        self.rho_gk   = np.mean(self.corr_mat['rho_gk'])
        self.rho_ba   = np.mean(self.corr_mat['rho_ba'])
        self.rho_ea   = np.mean(self.corr_mat['rho_ea'])
        self.rho_ga   = np.mean(self.corr_mat['rho_ga'])
        self.rho_eb   = np.mean(self.corr_mat['rho_eb'])
        self.rho_gb   = np.mean(self.corr_mat['rho_gb'])
        self.rho_ge   = np.mean(self.corr_mat['rho_ge'])

        
    def reload_quantum_param(self,baryons):
        """
        Load the quantum parameters that multiply the A,B,E,G in Hamiltonian
        """
        # get the original quantum numbers and experimental data
        state,sum_mass,J_tot,S_tot,L_tot,I_tot,SU_tot,HO_n,SL,ModEx = bs.states_mass(baryons)
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

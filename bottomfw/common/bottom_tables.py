#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
---------------------------------------------------------------
 Authors: A. Ramirez-Morales (andres.ramirez.morales@cern.ch)
          H. Garcia-Tecocoatzi
 ---------------------------------------------------------------
"""
import sys
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import pandas as pd
# framework includes
from bottomfw.common import data_utils as du


class BottomTables:
    """
    Class to produce latex tables containing mass and decay widths
    """
    def __init__(self, baryons, workpath='.', batch_results=False):
        self.m_workpath = workpath
        self.m_baryons = baryons
        self.m_batch = batch_results
        self.m_load_data(baryons)        

    def single_model_table(self):
        """
        Diquark or three quark to table to appear in the paper
        """
        if not os.path.exists(self.m_workpath+"/tables/"):
            os.mkdir(self.m_workpath+"/tables/")
        f_paper = open(self.m_workpath+'/tables/masses_'+self.m_baryons+'_paper.tex', "w")
        print("\\begin{tabular}{c | c  c c c }\hline \hline ", file=f_paper)
        print(" State     & Predicted Mass   & Experimental Mass & Predicted Width & Experimental Width   \\\ ", file=f_paper)
        print("           &      (MeV)       &    (MeV)          &      (MeV)      & $\Gamma_{tot}$ (MeV) \\\ \hline", file=f_paper)

        for i in range(len(self.m_mass)):

            if self.m_SU_tot[i] > 3 and self.m_SU_tot[i] < 3.5 : SU_tot_val = 10/3 # horrible bug fix
            else: SU_tot_val = 4/3
            
            quantum_state = du.name_quantum_state(self.m_baryons, self.m_J_tot[i], self.m_S_tot[i], self.m_L_tot[i], self.m_ModEx[i], SU_tot_val)
            exp_mass,exp_width= du.exp_mass_width(self.m_baryons, self.m_J_tot[i], self.m_S_tot[i], self.m_L_tot[i], self.m_ModEx[i], SU_tot_val)
            exp_mass_val,exp_mass_err_val= du.exp_mass_val(self.m_baryons,  self.m_J_tot[i], self.m_S_tot[i], self.m_L_tot[i], self.m_ModEx[i], SU_tot_val)

            print(quantum_state,'$',round(self.m_mass[i]),'^{+',round(self.m_error_up[i]),'}_{-',round(self.m_error_dn[i]),'}$',  '& ', exp_mass, ' & $',
                  round(self.m_decay[i],1),'^{+', round(self.m_decay_up[i],1),'}_{-', round(self.m_decay_dn[i],1),'}$', ' &', exp_width, '\\\ ', file=f_paper)

        name = self.m_baryons
        label = 'three_quark'
        print('\hline \hline', file=f_paper)
        print('\end{tabular}', file=f_paper)
        print("\caption{Every quantity is in MeV, except for percentage differences. States: $", self.m_baryons, "$}",file=f_paper)
        print("\label{tab:"+name+"_mass_"+label+"}", file=f_paper)

        
    def combined_model_table(self):
        """
        Combined three-quark and diquark
        """
        if not os.path.exists(self.m_workpath+"/tables/"):
            os.mkdir(self.m_workpath+"/tables/")
        f_paper = open(self.m_workpath+'/tables/masses_'+self.m_baryons+'_paper.tex', "w")
        baryon_name = du.baryon_name(self.m_baryons)
        flavor_name = du.flavor_label(self.m_baryons)
        
        print("\\begin{tabular}{c  c| c c c c c }\hline \hline", file=f_paper)
        print("            &  & Three-quark &  Quark-diquark    &               &              &  \\\ ", file=f_paper)
        print(baryon_name+"&  & Predicted   &    Predicted   &  Experimental &  Predicted            & Experimental \\\ ", file=f_paper)
        print(flavor_name+"  & $^{2S+1}L_{J}$ & Mass (MeV)  &   Mass (MeV)   &  Mass (MeV)   &  $\Gamma_{tot}$ (MeV) & $\Gamma$ (MeV) \\\ \hline", file=f_paper)

        s_wave_count,p_wave_count,d_wave_count=0,0,0
        for i in range(len(self.m_mass)):

            if self.m_HO_n[i] == 0:                
                if s_wave_count==0:
                    s_wave_count+=1
                    print('\hline', file=f_paper)
                    print(" $N=0$  &  &  &  &  &  \\\ ", file=f_paper)
            elif self.m_HO_n[i] == 1:
                if p_wave_count==0:
                    p_wave_count+=1
                    print('\hline', file=f_paper)
                    print(" $N=1$  &  &  &  &  &  \\\ ", file=f_paper)
            elif self.m_HO_n[i] == 2:
                if d_wave_count==0:
                    d_wave_count+=1
                    print('\hline', file=f_paper)
                    print(" $N=2$  &  &  &  &  &  \\\ ", file=f_paper)


            if self.m_SU_tot[i] > 3 and self.m_SU_tot[i] < 3.5 : SU_tot_val = 10/3 # horrible bug fix
            else: SU_tot_val = 4/3

            quantum_state = du.name_quantum_state(self.m_baryons, self.m_J_tot[i], self.m_S_tot[i], self.m_L_tot[i], self.m_ModEx[i], SU_tot_val)
            exp_mass,exp_width= du.exp_mass_width(self.m_baryons, self.m_J_tot[i], self.m_S_tot[i], self.m_L_tot[i], self.m_ModEx[i], SU_tot_val)
            exp_mass_val,exp_mass_err_val= du.exp_mass_val(self.m_baryons, self.m_J_tot[i], self.m_S_tot[i], self.m_L_tot[i], self.m_ModEx[i], SU_tot_val)

            wave_label= du.wave_label(self.m_S_tot[i], self.m_J_tot[i], self.m_L_tot[i])

            mass_lat = '$'+str(abs(round(self.m_mass[i])))+'^{+'+str(abs(round(self.m_error_up[i])))+'}_{-'+str(abs(round(self.m_error_dn[i])))+'}$'
            decay_lat = '$'+str(round(self.m_decay[i],1))+'^{+'+ str(abs(round(self.m_decay_up[i],1)))+'}_{-'+str(abs(round(self.m_decay_dn[i],1)))+'}$'
    
            mass_di_lat = "$\\dagger\\dagger$"
            for j in range(len(self.m_mass_di)): # match lambda states from the di-quark
                if self.m_SU_tot_di[j] > 3 and self.m_SU_tot_di[j] < 3.5 : SU_tot_di_val = 10/3 # horrible bug fix
                else: SU_tot_di_val = 4/3
                if self.m_J_tot[i] == self.m_J_tot_di[j] and self.m_S_tot[i]==self.m_S_tot_di[j] and self.m_L_tot[i]==self.m_L_tot_di[j] and self.m_ModEx[i]==self.m_ModEx_di[j] and SU_tot_val==SU_tot_di_val:
                    mass_di_lat = '$'+str(abs(round(self.m_mass_di[j])))+'^{+'+str(abs(round(self.m_error_up_di[j])))+'}_{-'+str(abs(round(self.m_error_dn_di[j])))+'}$'
                    break
        
            print(quantum_state, wave_label,'&', mass_lat, '&', mass_di_lat,'&', exp_mass, '&', decay_lat,'&', exp_width, '\\\ ', file=f_paper)
        
    
        name = self.m_baryons
        label = 'three_di_quark'
        print('\hline \hline', file=f_paper)
        print('\end{tabular}', file=f_paper)
        #print("\caption{Every quantity is in MeV, except for percentage differences. States: $", baryons, "$}",file=f_paper)
        #print("\label{tab:"+name+"_mass_"+label+"}", file=f_paper)
        f_paper.close()

        
    def decay_indi_table(self):
        
        df = pd.read_csv(self.m_workpath+'/tables/decays_indi_'+self.m_baryons+'_summary.csv')
        if not os.path.exists(self.m_workpath+"/tables/"):
            os.mkdir(self.m_workpath+"/tables/")
        f_decay_indi = open(self.m_workpath+'/tables/decay_indi_'+self.m_baryons+'_paper.tex', "w")
        n_decay_channels = int((len(df.columns)-8)/3)
        baryons = self.m_baryons

        import decay_utils as dec
        name_decays=[]
        name_decays.append('State')
        for k in range(n_decay_channels):
            name_decays.append(dec.latex_decay_label(baryons,k+1))
        name_decays.append('Tot $\\Gamma$')
        dec.print_header_latex(name_decays, f_decay_indi)

        for i in range(len(self.m_SU_tot)):
            channel_widths = []
            errors_up = []
            errors_dn = []
            for k in range(n_decay_channels):
                channel_widths.append(df['decay_'+str(k)][i])
                errors_up.append(df['dec_up_'+str(k)][i])
                errors_dn.append(df['dec_dn_'+str(k)][i])

            if self.m_SU_tot[i] > 3 and self.m_SU_tot[i] < 3.5 : SU_tot_val = 10/3
            else: SU_tot_val = 4/3
            quantum_state = du.name_quantum_state(self.m_baryons, self.m_J_tot[i],
                                                  self.m_S_tot[i], self.m_L_tot[i],
                                                  self.m_ModEx[i], SU_tot_val)
            wave_label= du.wave_label(self.m_S_tot[i], self.m_J_tot[i], self.m_L_tot[i])+'&'
            dec.print_row_latex(wave_label, channel_widths, errors_up, errors_dn, f_decay_indi)

        dec.print_bottom_latex(baryons, f_decay_indi)


    def latex_string_value_error(self, value, decimals=2, units='Mev'):
        """
        Method to calculate the mean value, the asymmetric error and return a latex string
        """
        N = len(value)
        qntl_up = int(N*0.975) # 95% C.L.
        qntl_dn = int(N*0.025)
        qntl_up = int(N*0.1587) # 68% C.L.
        qntl_dn = int(N*0.8413)
        
        return '$'+str(round(np.mean(value),2)) +'^{+'+ str(abs(round(np.sort(value)[qntl_dn-1] - np.mean(value),2))) +'}_{-'+\
            str(abs(round(np.sort(value)[qntl_up-1] - np.mean(value), 2)))+'}$ '+units

        
    def parameter_combined(self):
        """
        Method to produce the combined parameter table
        """     
        M1 = self.latex_string_value_error(self.m_sampled_m1, decimals=2, units='MeV')
        M2 = self.latex_string_value_error(self.m_sampled_m2, decimals=2, units='MeV')
        M3 = self.latex_string_value_error(self.m_sampled_m3, decimals=2, units='MeV')                                   
        K  = self.latex_string_value_error(self.m_sampled_k,  decimals=4, units='GeV$^{3}$')                                  
        A  = self.latex_string_value_error(self.m_sampled_a, decimals=2, units='MeV')
        B  = self.latex_string_value_error(self.m_sampled_b, decimals=2, units='MeV')
        E  = self.latex_string_value_error(self.m_sampled_e, decimals=2, units='MeV')
        G  = self.latex_string_value_error(self.m_sampled_g, decimals=2, units='MeV')

        Md1 = self.latex_string_value_error(self.m_sampled_md1,  decimals=2, units='MeV')
        Md2 = self.latex_string_value_error(self.m_sampled_md2,  decimals=2, units='MeV')
        Md3 = self.latex_string_value_error(self.m_sampled_md3,  decimals=2, units='MeV')
        MB  = self.latex_string_value_error(self.m_sampled_mb ,  decimals=2, units='MeV')
        K_di= self.latex_string_value_error(self.m_sampled_k_di, decimals=4, units='GeV$^{3}$')
        A_di= self.latex_string_value_error(self.m_sampled_a_di, decimals=2, units='MeV')
        B_di= self.latex_string_value_error(self.m_sampled_b_di, decimals=2, units='MeV')
        E_di= self.latex_string_value_error(self.m_sampled_e_di, decimals=2, units='MeV')
        G_di= self.latex_string_value_error(self.m_sampled_g_di, decimals=2, units='MeV')

        dd = '$\\dagger$'
        if not os.path.exists(self.m_workpath+"/tables/"):
            os.mkdir(self.m_workpath+"/tables/")
        f = open(self.m_workpath+'/tables/fit_parameters_combined.tex', "w")
        print("\\begin{tabular}{c | c c}\hline \hline", file=f)
        print("            &  three-quark & diquark \\\ ", file=f)
        print(" Parameter  &  Value       & Value    \\\ \hline", file=f)
        print(" $m_{b}$ &",                   M1, '&', MB,  "\\\ ", file=f)
        print(" $m_{s}$ &",                   M2, '&', dd,  "\\\ ", file=f)
        print(" $m_{n}$ &",                   M3, '&', dd,  "\\\ ", file=f)
        print(" $m_{D_{\\Omega}}$          &",dd, '&', Md1, "\\\ ", file=f)
        print(" $m_{D_{\\Xi}}$             &",dd, '&', Md2, "\\\ ", file=f)
        print(" $m_{D_{\\Sigma,\\Lambda}}$ &",dd, '&', Md3, "\\\ ", file=f)
        print(" $K_b$   &"                   ,K , '&', K_di,"\\\ ", file=f)
        print(" $P_S$     &"                   ,A , '&', A_di,"\\\ ", file=f)
        print(" $P_{SL}$     &"                   ,B , '&', B_di,"\\\ ", file=f)
        print(" $P_{I}$     &"                   ,E , '&', E_di,"\\\ ", file=f)
        print(" $P_{f}$     &"                   ,G , '&', G_di,"\\\ ", file=f)
        print("\hline\hline", file=f)
        print("\end{tabular}", file=f)
        print("\caption{Model fitted paremeters parameters.}",file=f)
        print("\label{tab:comb_fit}", file=f)
        f.close()

    def correlation_table_three(self):
        """
        Method to write correlation matrix table for the three quark system
        """
        if not os.path.exists(self.m_workpath+"/tables/"):
            os.mkdir(self.m_workpath+"/tables/")
        f = open(self.m_workpath+'/tables/correlation_3quark.tex', "w")
        print("\\begin{tabular}{c  c  c  c  c  c  c  c  c}\hline \hline", file=f)
        print("         &  $m_{b}$       &     $m_{s}$    &    $m_{n}$  &      $K_b$    & $P_S$ & $P_{SL}$ & $P_I$ & $P_f$ \\\ \hline", file=f)
        print(" $m_{b}$ &     1   &   &   &   &    &   &   &  \\\ ", file=f)
        print(" $m_{s}$ &",self.m_rho_m2m1, "&  1   &   &   &   &   &   &  \\\ ", file=f)
        print(" $m_{n}$ &",self.m_rho_m3m1, "&", self.m_rho_m3m2,"&  1   &   &   &   &   & \\\ ", file=f)
        print(" $K_b$   &",self.m_rho_km1,  "&", self.m_rho_km2 ,"&", self.m_rho_km3, "&  1   &   &   &   &   \\\ ", file=f)
        print(" $P_S$     &",self.m_rho_am1,  "&", self.m_rho_am2 ,"&", self.m_rho_am3, "&", self.m_rho_ak,"& 1 &   &   & \\\ ", file=f)
        print(" $P_{SL}$     &",self.m_rho_bm1,  "&", self.m_rho_bm2 ,"&", self.m_rho_bm3, "&", self.m_rho_bk,"&",self.m_rho_ba,"& 1  &   & \\\ ", file=f)
        print(" $P_I$     &",self.m_rho_em1,  "&", self.m_rho_em2 ,"&", self.m_rho_em3, "&", self.m_rho_ek,"&",self.m_rho_ea,"&",self.m_rho_eb,"& 1  &  \\\ ", file=f)
        print(" $P_f$     &",self.m_rho_gm1,  "&", self.m_rho_gm2 ,"&", self.m_rho_gm3, "&", self.m_rho_gk,"&",self.m_rho_ga,"&",self.m_rho_gb,"&",self.m_rho_ge,"& 1 \\\ \hline \hline", file=f)
        print('\end{tabular}', file=f)
        print("\caption{Correlation between fitted parameters, three-quark stystem }",file=f)
        print("\label{tab:3quark_corr}", file=f)
        f.close()

    def correlation_table_di(self):
        """
        Method to write correlation matrix table for the diquark system    
        """                
        md1="$m_{D_{\\Omega}}$"
        md2="$m_{D_{\\Xi}}$"
        md3="$m_{D_{\\Sigma,\\Lambda}}$"
        if not os.path.exists(self.m_workpath+"/tables/"):
            os.mkdir(self.m_workpath+"/tables/")
        f = open(self.m_workpath+'/tables/correlation_diquark.tex', "w")
        print("\\begin{tabular}{c | c c c c c c c c c}\hline \hline", file=f)
        print("         &  $m_{b}$",  "&",  md1, "&",  md2,  "&", md3  ,"&  $K_b$   & $P_S$ & $P_{SL}$ & $P_{I}$ & $P_f$ \\\ \hline", file=f)
        print(" $m_{b}$ &     1   &   &   &   &    &   &   &  &  \\\ ", file=f)
        print(md1,     "&",self.m_rho_md2md1, "&  1   &   &   &   &    &   &   &  \\\ ", file=f)
        print(md2,     "&",self.m_rho_md3md1, "&",self.m_rho_md3md2,"&  1   &   &   &   &    &   & \\\ ", file=f)
        print(md3,     "&",self.m_rho_mbmd1,  "&",self.m_rho_mbmd2 ,"&",self.m_rho_mbmd3, "&  1   &   &   &   &    &   \\\ ", file=f)
        print(" $K_b$   &",self.m_rho_kmd1 ,  "&",self.m_rho_kmd2  ,"&",self.m_rho_kmd3 , "&", self.m_rho_kmb,"& 1   &   &   &   & \\\ ", file=f)
        print(" $P_S$     &",self.m_rho_amd1 ,  "&",self.m_rho_amd2  ,"&",self.m_rho_amd3 , "&", self.m_rho_amb,"&",self.m_rho_ak,"& 1   &   &   & \\\ ", file=f)
        print(" $P_{SL}$     &",self.m_rho_bmd1 ,  "&",self.m_rho_bmd2  ,"&",self.m_rho_bmd3 , "&", self.m_rho_bmb,"&",self.m_rho_bk,"&",self.m_rho_ba,"& 1   &   &  \\\ ", file=f)
        print(" $P_I$     &",self.m_rho_emd1 ,  "&",self.m_rho_emd2  ,"&",self.m_rho_emd3 , "&", self.m_rho_emb,"&",self.m_rho_ek,"&",self.m_rho_ea,"&",self.m_rho_eb,"& 1   & \\\ ", file=f)
        print(" $P_dddf$     &",self.m_rho_gmd1 ,  "&",self.m_rho_gmd2  ,"&",self.m_rho_gmd3 , "&", self.m_rho_gmb,"&",self.m_rho_gk,"&",self.m_rho_ga,"&",self.m_rho_gb,"&", self.m_rho_ge, "&","1 \\\ \hline \hline", file=f) 
        print('\end{tabular}', file=f)
        print("\caption{Correlation between fitted parameters, diquark system.}",file=f)
        print("\label{tab:diquark_corr}", file=f)
        f.close()


    def m_load_data(self, baryons):
        """
        Method to load the data for all the tables
        -- three and diquark both already computed
        """
        data_frame = pd.read_csv(self.m_workpath+"/tables/masses_" + baryons + "_summary.csv")
        self.m_mass=        data_frame['mass']
        self.m_error_up=    data_frame['error_up']
        self.m_error_dn=    data_frame['error_dn']
        self.m_exp_mass=    data_frame['exp_mass']
        self.m_exp_mass_err=data_frame['exp_mass_err']
        self.m_decay=       data_frame['decay']
        self.m_decay_up=    data_frame['decay_up']
        self.m_decay_dn=    data_frame['decay_dn']
        self.m_J_tot=       data_frame['J_tot']
        self.m_S_tot=       data_frame['S_tot']
        self.m_L_tot=       data_frame['L_tot']
        self.m_ModEx=       data_frame['ModEx']
        self.m_HO_n=        data_frame['HO_n']
        self.m_SU_tot =     data_frame['SU_tot']
        
        data_frame_di = pd.read_csv(self.m_workpath+"/tables/masses_diquark_" + baryons + "_summary.csv")
        self.m_mass_di=        data_frame_di['mass']
        self.m_error_up_di=    data_frame_di['error_up']
        self.m_error_dn_di=    data_frame_di['error_dn']
        self.m_exp_mass_di=    data_frame_di['exp_mass']
        self.m_exp_mass_err_di=data_frame_di['exp_mass_err']
        self.m_decay_di=       data_frame_di['decay']
        self.m_decay_up_di=    data_frame_di['decay_up']
        self.m_decay_dn_di=    data_frame_di['decay_dn']
        self.m_J_tot_di=       data_frame_di['J_tot']
        self.m_S_tot_di=       data_frame_di['S_tot']
        self.m_L_tot_di=       data_frame_di['L_tot']
        self.m_ModEx_di=       data_frame_di['ModEx']
        self.m_HO_n_di=        data_frame_di['HO_n']
        self.m_SU_tot_di =     data_frame_di['SU_tot']

        if self.m_batch:
            all_files = glob.glob(os.path.join(self.m_workpath+"/batch_results/"+baryons+"/parameters/", "*.csv"))
            df_from_each_file = (pd.read_csv(f) for f in all_files)
            data_frame = pd.concat(df_from_each_file, ignore_index=True)            
        else:
            data_frame = pd.read_csv(self.m_workpath+"/tables/bootstrap_param_"+baryons+".csv")
        
        self.m_sampled_m1 = data_frame["M1"]
        self.m_sampled_m2 = data_frame["M2"]
        self.m_sampled_m3 = data_frame["M3"]
        self.m_sampled_k   = data_frame["K"].pow(2).div(pow(1000,3))
        self.m_sampled_a   = data_frame["A"]
        self.m_sampled_b   = data_frame["B"]
        self.m_sampled_e   = data_frame["E"]
        self.m_sampled_g   = data_frame["G"]

        data_frame_di = pd.read_csv(self.m_workpath+"/tables/bootstrap_param_diquark_"+baryons+".csv")
        self.m_sampled_md1  = data_frame_di["Md1"]
        self.m_sampled_md2  = data_frame_di["Md2"]
        self.m_sampled_md3  = data_frame_di["Md3"]
        self.m_sampled_mb   = data_frame_di["MB"]
        self.m_sampled_k_di = data_frame_di["K"].pow(2).div(pow(1000,3))
        self.m_sampled_a_di = data_frame_di["A"]
        self.m_sampled_b_di = data_frame_di["B"]
        self.m_sampled_e_di = data_frame_di["E"]
        self.m_sampled_g_di = data_frame_di["G"]

        if self.m_batch:
            all_files = glob.glob(os.path.join(self.m_workpath+"/batch_results/"+baryons+"/correlation/", "*.csv"))
            df_from_each_file = (pd.read_csv(f) for f in all_files)
            data_frame = pd.concat(df_from_each_file, ignore_index=True)            
        else:
            data_frame = pd.read_csv(self.m_workpath+"/tables/bootstrap_correlation_"+baryons+".csv")
            
        self.m_rho_m2m1 = round(np.mean(data_frame['rho_m2m1']), 2)
        self.m_rho_m3m1 = round(np.mean(data_frame['rho_m3m1']), 2)
        self.m_rho_km1  = round(np.mean(data_frame['rho_km1']), 2)
        self.m_rho_am1  = round(np.mean(data_frame['rho_am1']), 2)
        self.m_rho_bm1  = round(np.mean(data_frame['rho_bm1']), 2)
        self.m_rho_em1  = round(np.mean(data_frame['rho_em1']), 2)
        self.m_rho_gm1  = round(np.mean(data_frame['rho_gm1']), 2)
        self.m_rho_m3m2 = round(np.mean(data_frame['rho_m3m2']), 2)
        self.m_rho_km2  = round(np.mean(data_frame['rho_km2']), 2)
        self.m_rho_am2  = round(np.mean(data_frame['rho_am2']), 2)
        self.m_rho_bm2  = round(np.mean(data_frame['rho_bm2']), 2)
        self.m_rho_em2  = round(np.mean(data_frame['rho_em2']), 2)
        self.m_rho_gm2  = round(np.mean(data_frame['rho_gm2']), 2)
        self.m_rho_km3  = round(np.mean(data_frame['rho_km3']), 2)
        self.m_rho_am3  = round(np.mean(data_frame['rho_am3']), 2)
        self.m_rho_bm3  = round(np.mean(data_frame['rho_bm3']), 2)
        self.m_rho_em3  = round(np.mean(data_frame['rho_em3']), 2)
        self.m_rho_gm3  = round(np.mean(data_frame['rho_gm3']), 2)        
        self.m_rho_ak   = round(np.mean(data_frame['rho_ak']), 2)
        self.m_rho_bk   = round(np.mean(data_frame['rho_bk']), 2)
        self.m_rho_ek   = round(np.mean(data_frame['rho_ek']), 2)
        self.m_rho_gk   = round(np.mean(data_frame['rho_gk']), 2)
        self.m_rho_ba   = round(np.mean(data_frame['rho_ba']), 2)
        self.m_rho_ea   = round(np.mean(data_frame['rho_ea']), 2)
        self.m_rho_ga   = round(np.mean(data_frame['rho_ga']), 2)
        self.m_rho_eb   = round(np.mean(data_frame['rho_eb']), 2)
        self.m_rho_gb   = round(np.mean(data_frame['rho_gb']), 2)
        self.m_rho_ge   = round(np.mean(data_frame['rho_ge']), 2)

        if self.m_batch:
            all_files = glob.glob(os.path.join(self.m_workpath+"/batch_results_diquark/"+baryons+"/correlation/", "*.csv"))
            df_from_each_file = (pd.read_csv(f) for f in all_files)
            data_frame = pd.concat(df_from_each_file, ignore_index=True)            
        else:
            data_frame = pd.read_csv(self.m_workpath+"/tables/bootstrap_correlation_diquark_"+baryons+".csv")

        self.m_rho_md2md1 = round(np.mean(data_frame['rho_md2md1']), 2)
        self.m_rho_md3md1 = round(np.mean(data_frame['rho_md3md1']), 2)
        self.m_rho_mbmd1  = round(np.mean(data_frame['rho_mbmd1']), 2)
        self.m_rho_kmd1   = round(np.mean(data_frame['rho_kmd1']), 2)
        self.m_rho_amd1   = round(np.mean(data_frame['rho_amd1']), 2)
        self.m_rho_bmd1   = round(np.mean(data_frame['rho_bmd1']), 2)
        self.m_rho_emd1   = round(np.mean(data_frame['rho_emd1']), 2)
        self.m_rho_gmd1   = round(np.mean(data_frame['rho_gmd1']), 2)
        self.m_rho_md3md2 = round(np.mean(data_frame['rho_md3md2']), 2)
        self.m_rho_mbmd2  = round(np.mean(data_frame['rho_mbmd2']), 2)
        self.m_rho_kmd2   = round(np.mean(data_frame['rho_kmd2']), 2)
        self.m_rho_amd2   = round(np.mean(data_frame['rho_amd2']), 2)
        self.m_rho_bmd2   = round(np.mean(data_frame['rho_bmd2']), 2)
        self.m_rho_emd2   = round(np.mean(data_frame['rho_emd2']), 2)
        self.m_rho_gmd2   = round(np.mean(data_frame['rho_gmd2']), 2)
        self.m_rho_mbmd3  = round(np.mean(data_frame['rho_kmd3']), 2)
        self.m_rho_kmd3   = round(np.mean(data_frame['rho_kmd3']), 2)
        self.m_rho_amd3   = round(np.mean(data_frame['rho_amd3']), 2)
        self.m_rho_bmd3   = round(np.mean(data_frame['rho_bmd3']), 2)
        self.m_rho_emd3   = round(np.mean(data_frame['rho_emd3']), 2)
        self.m_rho_gmd3   = round(np.mean(data_frame['rho_gmd3']), 2)
        self.m_rho_kmb    = round(np.mean(data_frame['rho_kmb']), 2)
        self.m_rho_amb    = round(np.mean(data_frame['rho_amb']), 2)
        self.m_rho_bmb    = round(np.mean(data_frame['rho_bmb']), 2)
        self.m_rho_emb    = round(np.mean(data_frame['rho_emb']), 2)
        self.m_rho_gmb    = round(np.mean(data_frame['rho_gmb']), 2)
        self.m_rho_ak     = round(np.mean(data_frame['rho_ak']), 2)
        self.m_rho_bk     = round(np.mean(data_frame['rho_bk']), 2)
        self.m_rho_ek     = round(np.mean(data_frame['rho_ek']), 2)
        self.m_rho_gk     = round(np.mean(data_frame['rho_gk']), 2)
        self.m_rho_ba     = round(np.mean(data_frame['rho_ba']), 2)
        self.m_rho_ea     = round(np.mean(data_frame['rho_ea']), 2)
        self.m_rho_ga     = round(np.mean(data_frame['rho_ga']), 2)
        self.m_rho_eb     = round(np.mean(data_frame['rho_eb']), 2)
        self.m_rho_gb     = round(np.mean(data_frame['rho_gb']), 2)
        self.m_rho_ge     = round(np.mean(data_frame['rho_ge']), 2)

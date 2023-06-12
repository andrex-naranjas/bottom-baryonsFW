#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
---------------------------------------------------------------
 Authors: A. Ramirez-Morales (andres.ramirez.morales@cern.ch)
          H. Garcia-Tecocoatzi
---------------------------------------------------------------
"""
# data visalization module
import sys
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import pandas as pd
# framework includes
from bottomfw.common import data_utils as du
from bottomfw.baryons import bottom_states as bs
    

def plot(sample,name,xlabel,quark, states, workpath='.'):
    # sample plots
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.hist(sample, 100, density=True, label = 'Sampling')
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    y = norm.pdf(x, np.mean(sample), np.std(sample))
    p, alpha = du.normal_test(sample, alpha=0.05, verbose=False)        
    plt.plot(x,y,label='Gaussian fit')
    plt.text(0.15, 0.9,'$\\mu$={}, $\\sigma$={}'.format(round(np.mean(sample),1), round(np.std(sample),4)),
             ha='center', va='center', transform=ax.transAxes)
    plt.text(0.15, 0.8,'$p_{{val}}$={}, $\\alpha$={}'.format(round(p,3), alpha),
             ha='center', va='center', transform=ax.transAxes)
    plt.legend(frameon=False)
    plt.xlabel(xlabel)
    plt.ylabel('Arbitrary Units')
    plt.title(quark+' mesons')
    plt.savefig(workpath+'/plots/'+quark+'_bootstrap_'+name+'_'+states+'.pdf')
    plt.close()


def paper_tables_results(baryons, di_three_quark='diquark', decay_width=False,
                         asymmetric=False, prev_params=False, workpath='.', batch_number=None):
    # make statitical summary csv files for later use in latex and plots
    
    masses_df,decays_df=None,None
    n_decay_samples = 0
    di_label=''
    if di_three_quark == 'diquark':
        di_label = 'diquark_'
        if batch_number is not None:
            di_label = '_diquark'  
        # get the original quantum numbers and experimental data
        state,sum_mass,J_tot,S_tot,L_tot,I_tot,SU_tot,HO_n,SL,ModEx = bs.states_mass_diquark(baryons)
    else:
        state,sum_mass,J_tot,S_tot,L_tot,I_tot,SU_tot,HO_n,SL,ModEx = bs.states_mass(baryons)
        
    if batch_number is None:
        masses_df = pd.read_csv(workpath + "/tables/masses_states_" + di_label + baryons + ".csv")
        if decay_width and di_three_quark!='diquark':
            decays_df = pd.read_csv(workpath + "/tables/decays_states_" + di_label+baryons + ".csv")
            n_decay_samples = len(decays_df.index)
    else: # merge results from batch jobs
        all_files = glob.glob(os.path.join(workpath+"/batch_results"+di_label+"/"+baryons+"/mass_states/", "*.csv"))
        df_from_each_file = (pd.read_csv(f) for f in all_files)
        masses_df = pd.concat(df_from_each_file, ignore_index=True)
        if decay_width and di_three_quark!='diquark':
            all_files = glob.glob(os.path.join(workpath+"/batch_results"+di_label+"/"+baryons+"/decay_states/", "*.csv"))
            df_from_each_file = (pd.read_csv(f) for f in all_files)
            decays_df = pd.concat(df_from_each_file, ignore_index=True)
            n_decay_samples = len(decays_df.index)

    n_states  = len(masses_df.columns)
    n_samples = len(masses_df.index)

    print(n_samples, di_label)
    baryons_name = baryons
    if prev_params: baryons_name+='_previous'
    if di_three_quark=='diquark': di_label='diquark_'

    if not os.path.exists(workpath+"/tables/"):
        os.mkdir(workpath+"/tables/")    
    f_summary = open(workpath+'/tables/masses_'+di_label+baryons_name+'_summary.csv', "w")
    print('mass,error_up,error_dn,exp_mass,exp_mass_err,decay,decay_up,decay_dn,J_tot,S_tot,L_tot,SU_tot,ModEx,HO_n', file=f_summary)

    # if masses_df is not None and decays_df is not None:
    #     if n_samples != n_decay_samples and decay_width:
    #         sys.exit('Not the same number of masses of decays!!! something fishy is going on. bye')
    
    quantile_dn = int(n_samples*0.025)#1587)   #int(np.floor(N*0.1587))
    quantile_up = int(n_samples*0.975)#8413)+1 #int(np.floor(N*0.8413))    
    quantile_dn = int(n_samples*0.1587)
    quantile_up = int(n_samples*0.8413)
    
    for i in range(n_states):
        sorted_masses = np.sort(np.array(masses_df[str(i)+'_state']))
        # print(len(sorted_masses))
        mass = np.mean(sorted_masses)
        decay = 0
        if decay_width and di_three_quark!='diquark':
            sorted_decays = np.sort(np.array(decays_df[str(i)+'_state']))
            decay = np.mean(sorted_decays)
        
        # asymmetric error calculation via 68%(95%) quantile method
        error_up, error_dn = 0, 0
        if n_samples>1:
            error_up = sorted_masses[quantile_up-1] - mass
            error_dn = sorted_masses[quantile_dn-1] - mass

        up_decay,dn_decay=0,0
        if decay_width and di_three_quark!='diquark':
            if not np.isnan(decay):
                if n_decay_samples>1:
                    up_decay = sorted_decays[quantile_up-1] - decay
                    dn_decay = sorted_decays[quantile_dn-1] - decay
                
        exp_mass_val,exp_mass_err_val = du.exp_mass_val(baryons, J_tot[i], S_tot[i], L_tot[i], ModEx[i], SU_tot[i])
        
        if not asymmetric or n_samples<=1:                
            print(str(round(mass))+','+str(round(error_up))+','+str(round(error_dn))+','+str(exp_mass_val)+','\
                  + str(exp_mass_err_val) +','+ str(round(decay,1))+','+str(round(up_decay,1))+','+str(round(dn_decay,1))\
                  +','+str(J_tot[i])+','+str(S_tot[i])+','+str(L_tot[i])+','+str(SU_tot[i])+','+str(ModEx[i])+','+str(HO_n[i]), file=f_summary)             
        else:
            if not np.isnan(up_decay):                
                print(str(round(mass))+','+str(round(error_up))+','+str(round(error_dn))+','+str(exp_mass_val)+','\
                      + str(exp_mass_err_val) +','+ str(round(decay,1))+','+str(round(up_decay,1))+','+str(round(dn_decay,1))\
                      +','+str(J_tot[i])+','+str(S_tot[i])+','+str(L_tot[i])+','+str(SU_tot[i])+','+str(ModEx[i])+','+str(HO_n[i]), file=f_summary)
            else:
                print(str(round(mass))+','+ str(round((error_up+error_dn)/2))+','+ str(exp_mass_val)+','+ str(exp_mass_err_val)+','\
                      +str(J_tot[i])+','+str(S_tot[i])+','+str(L_tot[i])+','+str(ModEx[i])+','+str(HO_n[i]), file=f_summary)

                
def decay_indi_tables_results(baryons, decay_type="strong", asymmetric=False, prev_params=False, workpath='.', batch_number=None):
    '''
    Method to make statistical summary csv files for later use in latex tables for individual EM decay channels
    '''
    if decay_type=="strong":
        decay_name = ""
        second_name =""
    elif decay_type=="electro":
        decay_name = "em"
        second_name = "_"
    
    baryons_name = baryons
    if not os.path.exists(workpath+"/tables/"):
        os.mkdir(workpath+"/tables/")
    f_indi = open(workpath+'/tables/decays_indi_'+decay_name+second_name+baryons_name+'_summary.csv', "w")    
    state,sum_mass,J_tot,S_tot,L_tot,I_tot,SU_tot,HO_n,SL,ModEx = bs.states_mass(baryons)

    if baryons == "omegas" or baryons=="sigmas" or baryons=="cascades":            
        n_states = 9 # we only have up to P-wave states for electro CHECK!!
    else:
        n_states = 8 # we only have up to P-wave states for electro CHECK!!

    for i in range(len(state)):
        if (decay_type == "electro" and i >= n_states):
            break

        decay_indi_df = None
        if batch_number is None:
            decay_indi_df = pd.read_csv(workpath+"/tables/decays_indi"+second_name+decay_name+"/decays_state_"+str(i)+"_"+baryons+".csv")
        else: # merge results from batch jobs
            all_files = glob.glob(os.path.join(workpath+"/batch_results/"+baryons+"/decays_indi"+second_name+decay_name+"/state_"+str(i)+"/", "*.csv"))
            df_from_each_file = (pd.read_csv(f) for f in all_files)
            decay_indi_df = pd.concat(df_from_each_file, ignore_index=True)
            
        n_channels  = len(decay_indi_df.columns)
        n_samples   = len(decay_indi_df.index)
        # print(n_channels, n_samples)
        # input()
        quantile_dn = int(n_samples*0.025)#1587)   #int(np.floor(N*0.1587))
        quantile_up = int(n_samples*0.975)#8413)+1 #int(np.floor(N*0.8413))
            
        if(i==0): # print the header only once
            decays_header = ''
            for k in range(n_channels):
                decays_header+='decay_'+str(k)+',dec_up_'+str(k)+',dec_dn_'+str(k)+','
            print(decays_header+'dec_tot,dec_tot_up,dec_tot_dn,J_tot,S_tot,L_tot,SU_tot,ModEx,HO_n', file=f_indi)

        decays_value = ''
        total_decay, total_decay_up, total_decay_dn = 0,0,0
        for k in range(n_channels):
            sorted_decays_indi = np.sort(np.array(decay_indi_df[str(k)+'_channel']))
            decay_indi = np.mean(sorted_decays_indi)
            up_decay,dn_decay=0,0
            if not np.isnan(decay_indi):
                if n_samples>1:
                    up_decay = abs(sorted_decays_indi[quantile_up-1] - decay_indi)
                    dn_decay = abs(sorted_decays_indi[quantile_dn-1] - decay_indi)
            else:
                up_decay = np.nan
                dn_decay = np.nan
                
            decays_value += str(round(decay_indi,1))+','+str(round(up_decay,1))+','+str(round(dn_decay,1))+','
            total_decay += decay_indi
            total_decay_up += up_decay * up_decay
            total_decay_dn += dn_decay * dn_decay

        tot_str = str(round(total_decay,1))+','+str(round(np.sqrt(total_decay_up),1))+','+str(round(np.sqrt(total_decay_dn),1))
        quantum_str = str(J_tot[i])+','+str(S_tot[i])+','+str(L_tot[i])+','+str(SU_tot[i])+','+str(ModEx[i]+','+str(HO_n[i]))
        print(decays_value+tot_str+','+quantum_str, file=f_indi)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-----------------------------------------------------------------------
 code to obtain uncertainties of quarkonium mass spectrum via bootstrap
 authors: A. Ramirez-Morales (andres.ramirez.morales@cern.ch) and
          H. Garcia-Tecocoatzi
 ----------------------------------------------------------------------
"""
import sys
from iminuit import Minuit
import numpy as np
import datetime
# framework modules
from decays.decay_wrapper import decay
import decays.decay_utils as du
from decays.decay_width import DecayWidths


width = decay()

def model(v, gamma):
    return v*gamma

def least_squares(gamma):
    # y_var_0 = sigma_0 # best sigma_0 (omega=9.87) (cascade'=2.54) (cascades_anti3=5.35)
    # yvar_0 = y_var_0*np.ones(2)
    # yvar = y_errors_exp
    # yvar_2 = np.power(yvar_0, 2) + np.power(yvar, 2)
    yvar_2 = 1.0
    pred_m = model(param_v, gamma)
    yval_2 = np.power( (pred_m - exp_m), 2)
    return np.sum( np.divide(yval_2, yvar_2) )


def fit(least_squares):
    m = Minuit(least_squares, gamma=0)
    m.errordef=Minuit.LEAST_SQUARES
    m.migrad()
    return m

def sample_gauss(mu, sigma):
    return np.random.normal(mu, sigma, 10000)

def random(sample):
    return np.random.choice(sample, size=None)

def decay_coef():
    
    decay_test = DecayWidths()
    coef_decay=([])
    gamma = 1
    m1_in,m2_in,m3_in = 4926,466,300
    k_prim = 4848
    
    # baryons = 'omegas'
    # ModEx_val= 'lam'
    # #MassA = 3.0656*1000 # santopinto uses experimental values to calculate the width, which I think is wrong

    # MassA = 3.0074*1000 # gauss_3000
    # SA_qm = 1/2
    # LA_qm = 1
    # JA_qm = 1/2
    # SL_qm = 1
    # # get the sum of angular variables and h.o. integrals, gamma = 1
    # decay_coef = decay_test.total_decay_width(baryons, k_prim, MassA, SA_qm, LA_qm, JA_qm, SL_qm, ModEx_val, gamma_val=1, m1=m1_in, m2=m2_in, m3=m3_in)
    # coef_decay = np.append(coef_decay, decay_coef)

    # MassA = 3.049*1000 # gauss_3050
    # JA_qm = 1/2
    # SA_qm = 3/2
    # LA_qm = 1
    # SL_qm = 1
    # # get the sum of angular variables and h.o. integrals, gamma = 1
    # decay_coef = decay_test.total_decay_width(baryons, k_prim, MassA, SA_qm, LA_qm, JA_qm, SL_qm, ModEx_val, gamma_val=1,  m1=m1_in, m2=m2_in, m3=m3_in)
    # coef_decay = np.append(coef_decay, decay_coef)

    # MassA = 3.034*1000 # gauss_3065
    # JA_qm = 3/2
    # SA_qm = 1/2
    # LA_qm = 1
    # SL_qm = 1
    # # get the sum of angular variables and h.o. integrals, gamma = 1
    # decay_coef = decay_test.total_decay_width(baryons, k_prim, MassA, SA_qm, LA_qm, JA_qm, SL_qm, ModEx_val, gamma_val=1, m1=m1_in, m2=m2_in, m3=m3_in)
    # coef_decay = np.append(coef_decay, decay_coef)
    
    # MassA = 3.0764*1000 # gauss_3090
    # JA_qm = 3/2
    # SA_qm = 3/2
    # LA_qm = 1
    # SL_qm = 1
    # # get the sum of angular variables and h.o. integrals, gamma = 1
    # decay_coef = decay_test.total_decay_width(baryons, k_prim, MassA, SA_qm, LA_qm, JA_qm, SL_qm, ModEx_val, gamma_val=1, m1=m1_in, m2=m2_in, m3=m3_in)
    # coef_decay = np.append(coef_decay, decay_coef)

    # MassA = 3.122*1000 # gauss_3120  was not in the global fit
    # JA_qm = 5/2
    # SA_qm = 3/2
    # LA_qm = 1
    # SL_qm = 1
    # # get the sum of angular variables and h.o. integrals, gamma = 1
    # decay_coef = decay_test.total_decay_width(baryons, k_prim, MassA, SA_qm, LA_qm, JA_qm, SL_qm, ModEx_val, gamma_val=1, m1=m1_in, m2=m2_in, m3=m3_in)
    # coef_decay = np.append(coef_decay, decay_coef)

    # cascades prime
    # baryons = 'cascades'
    # ModEx_val= 'grd'
    # MassA = 2.6405*1000 # gauss_2645
    # JA_qm = 3/2
    # SA_qm = 3/2
    # LA_qm = 0
    # SL_qm = 1
    # # get the sum of angular variables and h.o. integrals, gamma = 1
    # decay_coef = decay_test.total_decay_width(baryons, k_prim, MassA, SA_qm, LA_qm, JA_qm, SL_qm, ModEx_val, gamma_val=1, m1=m1_in, m2=m2_in, m3=m3_in)
    # coef_decay = np.append(coef_decay, decay_coef)

    # ModEx_val= 'lam'
    # MassA = 2.9616*1000 # gauss_2965
    # JA_qm = 3/2
    # SA_qm = 3/2
    # LA_qm = 1
    # SL_qm = 1
    # # get the sum of angular variables and h.o. integrals, gamma = 1
    # decay_coef = decay_test.total_decay_width(baryons, k_prim, MassA, SA_qm, LA_qm, JA_qm, SL_qm, ModEx_val, gamma_val=1, m1=m1_in, m2=m2_in, m3=m3_in)
    # coef_decay = np.append(coef_decay, decay_coef)

    # sigmas
    baryons = 'sigmas'
    ModEx_val= 'grd'
    MassA = 5.8131*1000  # gauss_5813 = 5813.1
    JA_qm = 1/2
    SA_qm = 1/2
    LA_qm = 0
    SL_qm = 1
    # get the sum of angular variables and h.o. integrals, gamma = 1
    decay_coef = decay_test.total_decay_width(baryons, k_prim, MassA, SA_qm, LA_qm, JA_qm, SL_qm, ModEx_val, gamma_val=1, m1=m1_in, m2=m2_in, m3=m3_in)
    coef_decay = np.append(coef_decay, decay_coef)

    ModEx_val= 'grd' 
    MassA = 5.837*1000 # gauss_5837 = 5837.0
    JA_qm = 3/2
    SA_qm = 3/2
    LA_qm = 0
    SL_qm = 1
    # get the sum of angular variables and h.o. integrals, gamma = 1
    decay_coef = decay_test.total_decay_width(baryons, k_prim, MassA, SA_qm, LA_qm, JA_qm, SL_qm, ModEx_val, gamma_val=1, m1=m1_in, m2=m2_in, m3=m3_in)
    coef_decay = np.append(coef_decay, decay_coef) # not included in the fit

    
    # lambdas
    # baryons = 'lambdas'
    # ModEx_val= 'lam'
    # MassA = 2.6162*1000  # gauss_2592
    # JA_qm = 1/2
    # SA_qm = 1/2
    # LA_qm = 1
    # SL_qm = 0
    # # get the sum of angular variables and h.o. integrals, gamma = 1
    # decay_coef = decay_test.total_decay_width(baryons, k_prim, MassA, SA_qm, LA_qm, JA_qm, SL_qm, ModEx_val, gamma_val=1, m1=m1_in, m2=m2_in, m3=m3_in)
    # coef_decay = np.append(coef_decay, decay_coef)
    
    # ModEx_val= 'lam'
    # MassA = 3.025*1000  # gauss_2880
    # JA_qm = 5/2
    # SA_qm = 1/2
    # LA_qm = 2
    # SL_qm = 0
    # # get the sum of angular variables and h.o. integrals, gamma = 1
    # decay_coef = decay_test.total_decay_width(baryons, k_prim, MassA, SA_qm, LA_qm, JA_qm, SL_qm, ModEx_val, gamma_val=1, m1=m1_in, m2=m2_in, m3=m3_in)
    # coef_decay = np.append(coef_decay, decay_coef) # not included in the fit
    
    # # cascades_anti3
    # baryons = 'cascades_anti3'
    # ModEx_val= 'lam'
    # MassA = 2.7878*1000 # gauss_2792
    # JA_qm = 1/2
    # SA_qm = 1/2
    # LA_qm = 1
    # SL_qm = 0
    # # get the sum of angular variables and h.o. integrals, gamma = 1
    # decay_coef = decay_test.total_decay_width(baryons, k_prim, MassA, SA_qm, LA_qm, JA_qm, SL_qm, ModEx_val, gamma_val=1, m1=m1_in, m2=m2_in, m3=m3_in)
    # coef_decay = np.append(coef_decay, decay_coef)

    # ModEx_val= 'lam'
    # MassA = 2.8149*1000 # gauss_2818
    # JA_qm = 3/2
    # SA_qm = 1/2
    # LA_qm = 1
    # SL_qm = 0
    # # get the sum of angular variables and h.o. integrals, gamma = 1
    # decay_coef = decay_test.total_decay_width(baryons, k_prim, MassA, SA_qm, LA_qm, JA_qm, SL_qm, ModEx_val, gamma_val=1, m1=m1_in, m2=m2_in, m3=m3_in)
    # coef_decay = np.append(coef_decay, decay_coef)
    
    return coef_decay


# arrays to store the sampled parameters
sampled_gamma=([])

# input parameters
param_v = decay_coef()

# start bootstrap
start = datetime.datetime.now()

# gaussian pdf with the measured value and with experimental uncertainty decay widths
# Omega states
gauss_3000 = sample_gauss(4.5,0.7)   # OK paper table
gauss_3050 = sample_gauss(0.8,0.223) # OK paper table
gauss_3065 = sample_gauss(3.5, 0.447213595) # OK paper table
gauss_3090 = sample_gauss(8.7, 1.280624847)
gauss_3120 = sample_gauss(60, 10)

# cascades prime states
gauss_2645 = sample_gauss(2.25, 0.41) # OK paper table
gauss_2965 = sample_gauss(14.1, 1.6) # OK paper table

# sigma states
gauss_5813 = sample_gauss(5.0, 0.5)
gauss_5837 = sample_gauss(10.4, 2.0)

# lambda
gauss_2592 = sample_gauss(2.6, 0.6)  # OK paper table
gauss_2880 = sample_gauss(5.6, 0.8)  # OK paper table

# cascades anti3plet
gauss_2792 = sample_gauss(9.5, 2.0) # OK paper table
gauss_2818 = sample_gauss(2.48, 0.5) # OK paper table

# insert optimization lines here decay_untils.py

# construct the simulated sampling distribution (bootstrap technique)
for _ in range(100):
    # measured and quark-sum sampled masses
    exp_m = np.array([
        # # omegas
        # random(gauss_3000),
        # random(gauss_3050),
        # random(gauss_3065),
        # random(gauss_3090),
        # random(gauss_3120),
        # # cascades prime
        # random(gauss_2645),
        # random(gauss_2965),
        # sigmas
        random(gauss_5813),
        random(gauss_5837), # not in the fit
        # # lambdas
        # random(gauss_2592),
        # random(gauss_2880),
        # # cascades anti3
        # random(gauss_2792),
        # random(gauss_2818),
        
    ])   
    # perform the parameter fitting (via minimizing squared distance)
    m = fit(least_squares)
    sampled_gamma = np.append(sampled_gamma, np.sqrt(m.values['gamma']))

    
sampled = sampled_gamma.mean()


print(sampled)
print(np.mean(sampled_gamma), np.std(sampled_gamma), 'valor chido')
print(round(sampled, 2), 'sampled', round(np.std(sampled_gamma), 2))
print(round(np.sqrt(sampled), 2), 'sampled root squared')

# gauss_3000 = sample_gauss(4.6,0.7)
# gauss_3090 = sample_gauss(8.7, 1.3)
# gauss_3120 = sample_gauss(60, 26)

# # cascades states
# gauss_2965 = sample_gauss(14.1, 1.6)

# # sigma states
# gauss_2800 = sample_gauss(75, 60)

# # lambda
# gauss_2592 = sample_gauss(2.6, 0.6)

# # cascades anti3plet
# gauss_2792 = sample_gauss(9.5, 2.0)
# gauss_2818 = sample_gauss(2.48, 0.5)

# #sampled = 424
# sampled=4
# print(round(sampled*sampled*param_v[0],  2), sampled, param_v[0] ,  'random(gauss_3000)') #random(gauss_3000)
# print(round(sampled*sampled*param_v[1],  2), sampled, param_v[1] ,  'random(gauss_3050)') #random(gauss_3050)
# print(round(sampled*sampled*param_v[2],  2), sampled, param_v[2] ,  'random(gauss_3065)') #random(gauss_3065)
# print(round(sampled*sampled*param_v[3],  2), sampled, param_v[3] ,  'random(gauss_3090)') #random(gauss_3090)
# print(round(sampled*sampled*param_v[4],  2), sampled, param_v[4] ,  'random(gauss_2645)') #random(gauss_2645)
# print(round(sampled*sampled*param_v[5],  2), sampled, param_v[5] ,  'random(gauss_2965)') #random(gauss_2965)
# print(round(sampled*sampled*param_v[6],  2), sampled, param_v[6] ,  'random(gauss_2455)') #random(gauss_2455)
# print(round(sampled*sampled*param_v[7],  2), sampled, param_v[7] ,  'random(gauss_2520)') #random(gauss_2520)
# print(round(sampled*sampled*param_v[8],  2), sampled, param_v[8] ,  'random(gauss_2592)') #random(gauss_2592)
# print(round(sampled*param_v[9],  2), sampled, param_v[9] ,  'random(gauss_2880)') #random(gauss_2880)
# print(round(sampled*param_v[10], 2), sampled, param_v[10],  'random(gauss_2792)') #random(gauss_2792)                                                             
# print(round(sampled*param_v[10], 2), sampled, param_v[11],  'random(gauss_2818)') #random(gauss_2818)
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        



end = datetime.datetime.now()
elapsed_time = end - start
print("Elapsed total time = " + str(elapsed_time))

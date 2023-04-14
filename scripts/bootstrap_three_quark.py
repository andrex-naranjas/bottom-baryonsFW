#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
------------------------------------------------------------------------------
 Script to obtain uncertainties of heavy mass spectrum and widhts via bootstrap
 Authors: A. Ramirez-Morales (andres.ramirez.morales@cern.ch) and
          H. Garcia-Tecocoatzi
-----------------------------------------------------------------------------
"""
import sys
import os
from iminuit import Minuit
import numpy as np
import datetime
import pandas as pd
# framework modules
from bottomfw.baryons import data_preparation as dp
from bottomfw.baryons.bottom_three_quark import BottomThreeQuark


if len(sys.argv) <= 1:
    sys.exit('Provide bottom states group name. Try again!')

#states = 'omega' # All, omega, cascades, sigmaLamb
run_baryons = sys.argv[1]

# for running batch jobs with htcondor
batch_number = None
if len(sys.argv) == 4:
    batch_number = sys.argv[2]
    workpath = sys.argv[3]
else:
    workpath = os.getcwd()

# input parameters
param_v,param_w,param_x,param_y,param_z,param_q1,param_q2,param_q3,\
    param_is_rho,param_is_lam,param_is_omega,param_is_cascade,param_is_sigma = dp.fetch_data_extended()


def model(q1, q2, q3, is_rho, is_lam, is_omega, is_cascade, is_sigma, v, w, x, y, z, m1, m2, m3, k, a, b, e, g):
    return q1*m1 + q2*m2 + q3*m3 + \
        v*k*np.sqrt(1./(is_rho*(is_omega*m2 + is_cascade*((m2+m3)/2) + is_sigma*m3 ) + \
                        is_lam*(is_omega*((3*m2*m1)/(2.*m2+m1)) + is_cascade*((1.5*(m2+m3)*m1)/(m1+m2+m3)) + is_sigma*((3.*m3*m1)/(2.*m3+m1)) ) )) + \
                  w*a + x*b + y*e + z*g


def least_squares(m1, m2, m3, k, a, b, e, g):
    # y_var_0 = sigma_0 # best sigma_0=15.42
    # yvar_0 = y_var_0*np.ones(19)
    # yvar = y_errors_exp
    # yvar_2 = np.power(yvar_0, 2) + np.power(yvar, 2)
    yvar_2 = 0.001
    pred_m = model(param_q1, param_q2, param_q3, param_is_rho, param_is_lam,
                   param_is_omega, param_is_cascade, param_is_sigma, param_v,
                   param_w, param_x, param_y, param_z,
                   m1, m2, m3, k, a, b, e, g)
    yval_2 = np.power( (pred_m - exp_m), 2)
    return np.sum( np.divide(yval_2, yvar_2) )


def fit(least_squares):
    m = Minuit(least_squares, m1=1, m2=1, m3=1, k=0, a=0, b=0, e=0, g=0)#1400, m2=300, m3=250, k=0, a=0, b=0, e=0, g=0)
    m.limits['m1'] = (4000, 6000)
    m.limits['m2'] = (400, 470)
    m.limits['m3'] = (250, 300)
    m.errordef=Minuit.LEAST_SQUARES
    m.migrad()
    return m

def sample_gauss(mu, sigma):
    return np.random.normal(mu, sigma, 10000)

def random(sample, random_n=1):
    #return np.mean(resample(sample, replace=False, n_samples=1, random_state=random_n))
    return np.random.choice(sample, size=None)


# arrays to store the sampled parameters
sampled_k,sampled_a,sampled_b,sampled_e,sampled_g = ([]),([]),([]),([]),([])
sampled_m1,sampled_m2,sampled_m3 = ([]),([]),([])

# arrays to store sampled correlation coeficients
rho_m2m1,rho_m3m1,rho_km1,rho_am1,rho_bm1,rho_em1,rho_gm1 = ([]),([]),([]),([]),([]),([]),([])
rho_m3m2,rho_km2,rho_am2,rho_bm2,rho_em2,rho_gm2,rho_km3  = ([]),([]),([]),([]),([]),([]),([])
rho_am3,rho_bm3,rho_em3,rho_gm3,rho_ak,rho_bk,rho_ek      = ([]),([]),([]),([]),([]),([]),([])
rho_gk, rho_ba, rho_ea, rho_ga, rho_eb,rho_gb,rho_ge      = ([]),([]),([]),([]),([]),([]),([])

# start bootstrap
start = datetime.datetime.now()

sigma_model = 10**2 # to be obtained with optimization (Li.Jin)
# gaussian pdf with the measured value and with experimental and model(sigma_model) uncertainties
# Omega states
gauss_6061 = sample_gauss(6045.2, np.power((1.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 2702), PDG
gauss_6316 = sample_gauss(6315.6, np.power((1.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 2767), PDG
gauss_6330 = sample_gauss(6330.3, np.power((1.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 3015), PDG
gauss_6340 = sample_gauss(6339.7, np.power((1.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 3044), PDG
gauss_6350 = sample_gauss(6349.8, np.power((1.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 3051), PDG
# Cascade b sextet
gauss_5935 = sample_gauss(5935.0, np.power((1.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 2461), PDG Average
gauss_5953 = sample_gauss(5953.8, np.power((1.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 2796), PDG Average
gauss_6328 = sample_gauss(6227.9, np.power((1.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 2832), PDG Average
# Sigma b
gauss_5813 = sample_gauss(5813.1, np.power((1.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 2453), PDG Average
gauss_5837 = sample_gauss(5837.0, np.power((1.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 2517), PDG Average
gauss_6097 = sample_gauss(6096.9, np.power((1.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 2819), PDG Average
# Lambda b
gauss_5617 = sample_gauss(5619.6, np.power((1.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 2283), PDG
gauss_5912 = sample_gauss(5912.2, np.power((1.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 2649), PDG
gauss_5920 = sample_gauss(5920.1, np.power((1.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 2685), PDG
gauss_6146 = sample_gauss(6146.2, np.power((1.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 2283), PDG
gauss_6152 = sample_gauss(6152.5, np.power((1.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 2649), PDG
gauss_6070 = sample_gauss(6072.3, np.power((1.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 2685), PDG
# Cascades anti-3-plet
gauss_5794 = sample_gauss(5794.5, np.power((1.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 2570), PDG Average
gauss_6100 = sample_gauss(6100.0, np.power((1.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 2635), PDG Average
gauss_6327 = sample_gauss(6327.0, np.power((1.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 2934), LHCb PhysRevLett
gauss_6333 = sample_gauss(6333.0, np.power((1.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 2941), LHCb PhysRevLett

# plug here the sigma_0 optimization lines from data_utils.py

# construct the simulated sampling distribution (bootstrap technique)
for i in range(10): # max 10000 with decays included, computationally expensive
    #if(states=='All'):
    exp_m = np.array([ # measured baryon masses
        # omegas
        random(gauss_6061),
        random(gauss_6316),
        random(gauss_6330),
        random(gauss_6340),
        random(gauss_6350),
        # Cascade 
        random(gauss_5935),
        random(gauss_5953),
        random(gauss_6328),
        # Sigma b         
        random(gauss_5813),
        random(gauss_5837),
        random(gauss_6097),
        # Lambda b
        random(gauss_5617),
        random(gauss_5912),
        random(gauss_5920),
        # random(gauss_6146),
        # random(gauss_6152),
        # Cascades
        random(gauss_5794),
        random(gauss_6100),
        # random(gauss_6327),
        # random(gauss_6333)
    ])    
    # perform the parameter fitting (via minimizing squared distance)
    m = fit(least_squares)

    sampled_m1 = np.append(sampled_m1, m.values['m1'])
    sampled_m2 = np.append(sampled_m2, m.values['m2'])
    sampled_m3 = np.append(sampled_m3, m.values['m3'])

    sampled_k = np.append(sampled_k, m.values['k'])
    sampled_a = np.append(sampled_a, m.values['a'])
    sampled_b = np.append(sampled_b, m.values['b'])
    sampled_e = np.append(sampled_e, m.values['e'])
    sampled_g = np.append(sampled_g, m.values['g'])

    # if states != 'omega':
    # correlation matrix
    corr = m.covariance.correlation()

    rho_m2m1 = np.append(rho_m2m1, corr['m2','m1'])
    rho_m3m1 = np.append(rho_m3m1, corr['m3','m1'])
    rho_km1  = np.append(rho_km1,  corr['k','m1'])
    rho_am1  = np.append(rho_am1,  corr['a','m1'])
    rho_bm1  = np.append(rho_bm1,  corr['b','m1'])
    rho_em1  = np.append(rho_em1,  corr['e','m1'])
    rho_gm1  = np.append(rho_gm1,  corr['g','m1'])

    rho_m3m2  = np.append(rho_m3m2, corr['m3','m2'])
    rho_km2   = np.append(rho_km2 , corr['k','m2'])
    rho_am2   = np.append(rho_am2 , corr['a','m2'])
    rho_bm2   = np.append(rho_bm2 , corr['b','m2'])
    rho_em2   = np.append(rho_em2 , corr['e','m2'])
    rho_gm2   = np.append(rho_gm2 , corr['g','m2'])

    rho_km3  = np.append(rho_km3, corr['k','m3'])
    rho_am3  = np.append(rho_am3, corr['a','m3'])
    rho_bm3  = np.append(rho_bm3, corr['b','m3'])
    rho_em3  = np.append(rho_em3, corr['e','m3'])
    rho_gm3  = np.append(rho_gm3, corr['g','m3'])

    rho_ak  = np.append(rho_ak, corr['a','k'])
    rho_bk  = np.append(rho_bk, corr['b','k'])
    rho_ek  = np.append(rho_ek, corr['e','k'])
    rho_gk  = np.append(rho_gk, corr['g','k'])

    rho_ba  = np.append(rho_ba, corr['b','a'])
    rho_ea  = np.append(rho_ea, corr['e','a'])
    rho_ga  = np.append(rho_ga, corr['g','a'])

    rho_eb  = np.append(rho_eb, corr['e','b'])
    rho_gb  = np.append(rho_gb, corr['g','b'])
    
    rho_ge  = np.append(rho_ge, corr['g','e'])
    

# save bootstrap results
df = pd.DataFrame({"M1" : sampled_m1,"M2" : sampled_m2,"M3" : sampled_m3,
                   "K" : sampled_k,   "A" : sampled_a,
                   "B": sampled_b,    "E" : sampled_e, "G" : sampled_g})

if batch_number is None:
    if not os.path.exists(workpath+"/tables/"):
        os.makedirs(workpath+"/tables/")        
    df.to_csv(workpath+"/tables/bootstrap_param_"+run_baryons+".csv", index=False)
else:
    if not os.path.exists(workpath+"/batch_results/"+run_baryons+"/parameters/"):
        os.makedirs(workpath+"/batch_results/"+run_baryons+"/parameters/")
    df.to_csv(workpath+"/batch_results/"+run_baryons+"/parameters/"+str(batch_number)+".csv", index=False)

# create dictionaries
param   = {'q1':param_q1, 'q2':param_q2, 'q3':param_q3,'is_rho':param_is_rho, 'is_lam':param_is_lam,'is_omega':param_is_omega,
           'is_cascade':param_is_cascade, 'is_sigma':param_is_sigma,'V':param_v, 'W':param_w, 'X':param_x, 'Y':param_y, 'Z':param_z}
sampled = {'sampled_m1':sampled_m1,'sampled_m2':sampled_m2,'sampled_m3':sampled_m3,'sampled_k':sampled_k,
           'sampled_a':sampled_a, 'sampled_b':sampled_b, 'sampled_e':sampled_e, 'sampled_g':sampled_g}

corr_mat_ext ={'rho_m2m1':rho_m2m1, 'rho_m3m1':rho_m3m1, 'rho_km1':rho_km1, 'rho_am1':rho_am1, 'rho_bm1':rho_bm1, 'rho_em1':rho_em1, 'rho_gm1':rho_gm1,
               'rho_m3m2':rho_m3m2, 'rho_km2':rho_km2, 'rho_am2':rho_am2, 'rho_bm2':rho_bm2, 'rho_em2':rho_em2, 'rho_gm2':rho_gm2, 'rho_km3':rho_km3,
               'rho_am3':rho_am3, 'rho_bm3':rho_bm3, 'rho_em3':rho_em3, 'rho_gm3':rho_gm3, 'rho_ak':rho_ak, 'rho_bk':rho_bk, 'rho_ek':rho_ek,
               'rho_gk':rho_gk, 'rho_ba':rho_ba, 'rho_ea':rho_ea, 'rho_ga':rho_ga, 'rho_eb':rho_eb, 'rho_gb':rho_gb, 'rho_ge':rho_ge}

df = pd.DataFrame(corr_mat_ext)
if batch_number is None:
    if not os.path.exists(workpath+"/tables/"):
        os.makedirs(workpath+"/tables/")
    df.to_csv(workpath+"/tables/bootstrap_correlation_"+run_baryons+".csv", index=False)
else:
    if not os.path.exists(workpath+"/batch_results/"+run_baryons+"/correlation/"):
        os.makedirs(workpath+"/batch_results/"+run_baryons+"/correlation/")
    df.to_csv(workpath+"/batch_results/"+run_baryons+"/correlation/"+str(batch_number)+".csv", index=False)





# calculate the masses and decays using the bootstrap simulation above
# results = BottomThreeQuark(baryons=run_baryons, params=param, sampled=sampled, corr_mat=corr_mat_ext, asymmetric=True,
#                            decay_width=True, bootstrap_width=True, batch_number=batch_number, workpath=workpath)
# results.fetch_values()

# print('Getting paper results for:', run_baryons)
# #input()
# # omegas,cascades,sigmas,lambdas,cascades_anti3
# results.paper_results_predictions(bootstrap=True, bootstrap_width=True, prev_params=False) # all running for paper
# #results.paper_results_predictions(baryons=run_baryons,        bootstrap=True, bootstrap_width=True, prev_params=False, decay_width=True) # all running for paper
# # avoid time consuming decay widths calculations
# # results.paper_results_predictions(baryons=run_baryons,        bootstrap=True, bootstrap_width=False, prev_params=False, decay_width=False)
# # nominal results, expected to be the same as the previous paper, very important check
# # results.paper_results_predictions(baryons=run_baryons,         bootstrap=False, bootstrap_width=False, prev_params=True, decay_width=True)
# end = datetime.datetime.now()
# elapsed_time = end - start
# print("Elapsed total time = " + str(elapsed_time))


results = BottomThreeQuark(baryons=run_baryons, params=param, sampled=sampled, corr_mat=corr_mat_ext, asymmetric=True,
                           decay_width=False, bootstrap_width=False, decay_width_em=True, bootstrap_width_em=False, batch_number=batch_number, workpath=workpath)
results.fetch_values()

print('Getting paper results for:', run_baryons)
#input()
# omegas,cascades,sigmas,lambdas,cascades_anti3
results.paper_results_predictions(bootstrap=True, bootstrap_width=False, prev_params=False) # all running for paper
#results.paper_results_predictions(baryons=run_baryons,        bootstrap=True, bootstrap_width=True, prev_params=False, decay_width=True) # all running for paper
# avoid time consuming decay widths calculations
# results.paper_results_predictions(baryons=run_baryons,        bootstrap=True, bootstrap_width=False, prev_params=False, decay_width=False)
# nominal results, expected to be the same as the previous paper, very important check
# results.paper_results_predictions(baryons=run_baryons,         bootstrap=False, bootstrap_width=False, prev_params=True, decay_width=True)
end = datetime.datetime.now()
elapsed_time = end - start
print("Elapsed total time = " + str(elapsed_time))

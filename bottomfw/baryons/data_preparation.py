#!/usr/bin/env python3
# - * - coding: utf-8 - * -
"""
 Authors: A. Ramirez-Morales (andres.ramirez.morales@cern.ch)
          H. Garcia-Tecocoatzi
"""
import numpy as np

def fetch_data_extended():

    w_om=1.732050808 # omega_lam
    w_ls=1.732050808 # lambda&sigma_lam
    wcas=1.732050808 # casc_lam
    wcar=1.732050808 # wcas_rho

    # hamiltonian coeficients
    # from total: 3(in)-2(out)-2(in) for omega and 2(in)-2(out)-3(in) for cascade
    param_v          = np.array([0.00, w_om, w_om, w_om, w_om,   0.00, 0.00, wcas,   0.00, 0.00, w_ls,   0.00, w_ls, w_ls, w_ls, w_ls,   0.00, wcas, wcas, wcas]) # coef infront kprim_c
    param_w          = np.array([0.75, 0.75, 0.75, 3.75, 3.75,   0.75, 3.75, 3.75,   0.75, 3.75, 0.75,   0.75, 0.75, 0.75, 0.75, 0.75,   0.75, 0.75, 0.75, 0.75]) # coef infront A
    param_x          = np.array([0.00, -1.0, 0.50, -2.5, -1.0,   0.00, 0.00, 1.50,   0.00, 0.00, -1.0,   0.00, -1.0, 0.50, -1.5, 1.00,   0.00, 0.50, -1.5, 1.00]) # coef infront B
    param_y          = np.array([0.00, 0.00, 0.00, 0.00, 0.00,   0.75, 0.75, 0.75,   2.00, 2.00, 2.00,   0.00, 0.00, 0.00, 0.00, 0.00,   0.75, 0.75, 0.75, 0.75]) # coef infront E
    param_z          = np.array([10/3, 10/3, 10/3, 10/3, 10/3,   10/3, 10/3, 10/3,   10/3, 10/3, 10/3,   4/3 , 4/3 , 4/3 , 4/3 , 4/3 ,   4/3 , 4/3 , 4/3 , 4/3 ]) # coef infront G    
    # quark content
    param_q1         = np.array([1.00, 1.00, 1.00, 1.00, 1.00,   1.00, 1.00, 1.00,   1.00, 1.00, 1.00,   1.00, 1.00, 1.00, 1.00, 1.00,   1.00, 1.00, 1.00, 1.00]) # coef m1 heavy quark
    param_q2         = np.array([2.00, 2.00, 2.00, 2.00, 2.00,   1.00, 1.00, 1.00,   0.00, 0.00, 0.00,   0.00, 0.00, 0.00, 0.00, 0.00,   1.00, 1.00, 1.00, 1.00]) # coef m2 strange quark
    param_q3         = np.array([0.00, 0.00, 0.00, 0.00, 0.00,   1.00, 1.00, 1.00,   2.00, 2.00, 2.00,   2.00, 2.00, 2.00, 2.00, 2.00,   1.00, 1.00, 1.00, 1.00]) # coef m3 light quark
    # harmonic oscillator mode
    param_is_rho     = np.array([9.99, 0.00, 0.00, 0.00, 0.00,   9.99, 9.99, 0.00,   9.99, 9.99, 0.00,   9.99, 0.00, 0.00, 0.00, 0.00,   9.99, 0.00, 0.00, 0.00]) # coef kp
    param_is_lam     = np.array([9.99, 1.00, 1.00, 1.00, 1.00,   9.99, 9.99, 1.00,   9.99, 9.99, 1.00,   9.99, 1.00, 1.00, 1.00, 1.00,   9.99, 1.00, 1.00, 1.00]) # coef kp
    # kind of baryon
    param_is_omega   = np.array([1.00, 1.00, 1.00, 1.00, 1.00,   0.00, 0.00, 0.00,   0.00, 0.00, 0.00,   0.00, 0.00, 0.00, 0.00, 0.00,   0.00, 0.00, 0.00, 0.00]) # coef kp
    param_is_cascade = np.array([0.00, 0.00, 0.00, 0.00, 0.00,   1.00, 1.00, 1.00,   0.00, 0.00, 0.00,   0.00, 0.00, 0.00, 0.00, 0.00,   1.00, 1.00, 1.00, 1.00]) # coef kp
    param_is_sigma   = np.array([0.00, 0.00, 0.00, 0.00, 0.00,   0.00, 0.00, 0.00,   1.00, 1.00, 1.00,   1.00, 1.00, 1.00, 1.00, 1.00,   0.00, 0.00, 0.00, 0.00]) # coef kp

    
    param_v          = np.array([0.00, w_om, w_om, w_om, w_om,   0.00, 0.00, wcas,   0.00, 0.00, w_ls,   0.00, w_ls, w_ls,   0.00, wcas]) # coef infront kprim_c
    param_w          = np.array([0.75, 0.75, 0.75, 3.75, 3.75,   0.75, 3.75, 3.75,   0.75, 3.75, 0.75,   0.75, 0.75, 0.75,   0.75, 0.75]) # coef infront A
    param_x          = np.array([0.00, -1.0, 0.50, -2.5, -1.0,   0.00, 0.00, 1.50,   0.00, 0.00, -1.0,   0.00, -1.0, 0.50,   0.00, 0.50]) # coef infront B
    param_y          = np.array([0.00, 0.00, 0.00, 0.00, 0.00,   0.75, 0.75, 0.75,   2.00, 2.00, 2.00,   0.00, 0.00, 0.00,   0.75, 0.75]) # coef infront E
    param_z          = np.array([10/3, 10/3, 10/3, 10/3, 10/3,   10/3, 10/3, 10/3,   10/3, 10/3, 10/3,   4/3 , 4/3 , 4/3 ,   4/3 , 4/3 ]) # coef infront G    
    # quark content
    param_q1         = np.array([1.00, 1.00, 1.00, 1.00, 1.00,   1.00, 1.00, 1.00,   1.00, 1.00, 1.00,   1.00, 1.00, 1.00,   1.00, 1.00]) # coef m1 heavy quark
    param_q2         = np.array([2.00, 2.00, 2.00, 2.00, 2.00,   1.00, 1.00, 1.00,   0.00, 0.00, 0.00,   0.00, 0.00, 0.00,   1.00, 1.00]) # coef m2 strange quark
    param_q3         = np.array([0.00, 0.00, 0.00, 0.00, 0.00,   1.00, 1.00, 1.00,   2.00, 2.00, 2.00,   2.00, 2.00, 2.00,   1.00, 1.00]) # coef m3 light quark
    # harmonic oscillator mode
    param_is_rho     = np.array([9.99, 0.00, 0.00, 0.00, 0.00,   9.99, 9.99, 0.00,   9.99, 9.99, 0.00,   9.99, 0.00, 0.00,   9.99, 0.00]) # coef kp
    param_is_lam     = np.array([9.99, 1.00, 1.00, 1.00, 1.00,   9.99, 9.99, 1.00,   9.99, 9.99, 1.00,   9.99, 1.00, 1.00,   9.99, 1.00]) # coef kp
    # kind of baryon
    param_is_omega   = np.array([1.00, 1.00, 1.00, 1.00, 1.00,   0.00, 0.00, 0.00,   0.00, 0.00, 0.00,   0.00, 0.00, 0.00,   0.00, 0.00]) # coef kp
    param_is_cascade = np.array([0.00, 0.00, 0.00, 0.00, 0.00,   1.00, 1.00, 1.00,   0.00, 0.00, 0.00,   0.00, 0.00, 0.00,   1.00, 1.00]) # coef kp
    param_is_sigma   = np.array([0.00, 0.00, 0.00, 0.00, 0.00,   0.00, 0.00, 0.00,   1.00, 1.00, 1.00,   1.00, 1.00, 1.00,   0.00, 0.00]) # coef kp

    return param_v,param_w,param_x,param_y,param_z,param_q1,param_q2,param_q3,\
        param_is_rho,param_is_lam,param_is_omega,param_is_cascade,param_is_sigma



def fetch_data_diquark(fit_type="trad"):

    w_om=1.732050808 # omega_lam
    w_ls=1.732050808 # lambda&sigma_lam
    wcas=1.732050808 # casc_lam
    wcar=1.732050808 # wcas_rho

    # hamiltonian coeficients
    param_v          = np.array([0.00, w_om, w_om, w_om, w_om,   0.00, 0.00, wcas,   0.00, 0.00, w_ls,   0.00, w_ls, w_ls,    0.00, wcas]) # coef infront kprim_c  w_ls, w_ls,
    param_w          = np.array([0.75, 0.75, 0.75, 3.75, 3.75,   0.75, 3.75, 3.75,   0.75, 3.75, 0.75,   0.75, 0.75, 0.75,    0.75, 0.75]) # coef infront A        0.75, 0.75,
    param_x          = np.array([0.00, -1.0, 0.50, -2.5, -1.0,   0.00, 0.00, 1.50,   0.00, 0.00, -1.0,   0.00, -1.0, 0.50,    0.00, 0.50]) # coef infront B        -1.5, 1.00,
    param_y          = np.array([0.00, 0.00, 0.00, 0.00, 0.00,   0.75, 0.75, 0.75,   2.00, 2.00, 2.00,   0.00, 0.00, 0.00,    0.75, 0.75]) # coef infront E        0.00, 0.00,
    param_z          = np.array([10/3, 10/3, 10/3, 10/3, 10/3,   10/3, 10/3, 10/3,   10/3, 10/3, 10/3,   4/3 , 4/3 , 4/3 ,    4/3 , 4/3 ]) # coef infront G        4/3 , 4/3 ,
        
    # kind of baryon
    if fit_type=="trad":
        param_is_omega   = np.array([1.00, 1.00, 1.00, 1.00, 1.00,   0.00, 0.00, 0.00,   0.00, 0.00, 0.00,   0.00, 0.00, 0.00,  0.00, 0.00]) # coef kp
        param_is_cascade = np.array([0.00, 0.00, 0.00, 0.00, 0.00,   1.00, 1.00, 1.00,   0.00, 0.00, 0.00,   0.00, 0.00, 0.00,  1.00, 1.00]) # coef kp
        param_is_sigma   = np.array([0.00, 0.00, 0.00, 0.00, 0.00,   0.00, 0.00, 0.00,   1.00, 1.00, 1.00,   1.00, 1.00, 1.00,  0.00, 0.00]) # coef kp
        return param_v,param_w,param_x,param_y,param_z,param_is_omega,param_is_cascade,param_is_sigma
            
    param_is_omega     = np.array([1.00, 1.00, 1.00, 1.00, 1.00,   0.00, 0.00, 0.00,   0.00, 0.00, 0.00,   0.00, 0.00, 0.00,   0.00, 0.00]) # coef kp            0.00, 0.00,
    param_is_cascade_p = np.array([0.00, 0.00, 0.00, 0.00, 0.00,   1.00, 1.00, 1.00,   0.00, 0.00, 0.00,   0.00, 0.00, 0.00,   0.00, 0.00]) # coef kp            0.00, 0.00,
    param_is_sigma     = np.array([0.00, 0.00, 0.00, 0.00, 0.00,   0.00, 0.00, 0.00,   1.00, 1.00, 1.00,   0.00, 0.00, 0.00,   0.00, 0.00]) # coef kp            0.00, 0.00,
    param_is_lambda    = np.array([0.00, 0.00, 0.00, 0.00, 0.00,   0.00, 0.00, 0.00,   0.00, 0.00, 0.00,   1.00, 1.00, 1.00,   0.00, 0.00]) # coef kp            1.00, 1.00,
    param_is_cascade   = np.array([0.00, 0.00, 0.00, 0.00, 0.00,   0.00, 0.00, 0.00,   0.00, 0.00, 0.00,   0.00, 0.00, 0.00,   1.00, 1.00]) # coef kp            0.00, 0.00,
                
    #return param_v,param_w,param_x,param_y,param_z,param_is_omega,param_is_cascade,param_is_sigma
    return param_v, param_w, param_x, param_y, param_z, param_is_omega, param_is_cascade_p, param_is_sigma, param_is_lambda, param_is_cascade

    
def hamiltonian_quantum_factors(state, sum_mass, J_tot, S_tot, L_tot, I_tot, SU_tot, HO_n, ModEx):
    # compute the (quantum) coeficients that multiplied the fitted/bootstrap parameters
    # this function needs to be fed directly with quantum numbers storeds in a numpy.array
    
    param_v,param_w,param_x,param_y,param_z,mass_sum = [],[],[],[],[],[]
    
    for i in range(len(sum_mass)):
        n_ho = HO_n[i]
        omega_ho=0
        omega_mx=0
        if(ModEx[i]=='lam' or ModEx[i]=='rpl'):
            if(state[i]=='omg'):   omega_ho=0.058892513
            elif(state[i]=='cas'): omega_ho=0.062695053
            elif(state[i]=='tri'): omega_ho=0.062695053
            elif(state[i]=='sig'): omega_ho=0.068087711
            elif(state[i]=='Lam'): omega_ho=0.068087711
        elif(ModEx[i]=='rho' or ModEx[i]=='rpr'):
            if(state[i]=='omg'):   omega_ho=0.081649658
            elif(state[i]=='cas'): omega_ho=0.089742361
            elif(state[i]=='tri'): omega_ho=0.089742361
            elif(state[i]=='sig'): omega_ho=0.100843897
            elif(state[i]=='Lam'): omega_ho=0.100843897
        if(ModEx[i]=='mix'):
            n_ho = 1
            if(state[i]=='omg'):
                omega_ho=0.058892513
                omega_mx=0.081649658
            elif(state[i]=='cas'):
                omega_ho=0.062695053
                omega_mx=0.089742361
            elif(state[i]=='tri'):
                omega_ho=0.062695053
                omega_mx=0.089742361
            elif(state[i]=='sig'):
                omega_ho=0.068087711
                omega_mx=0.100843897
            elif(state[i]=='Lam'):
                omega_ho=0.068087711
                omega_mx=0.100843897
                
        param_v.append(n_ho * omega_ho + n_ho * omega_mx)    # coef infront kprim 
        param_w.append((S_tot[i] + 1) * S_tot[i]) # coef infront A
        param_x.append(0.5 * ( (J_tot[i] + 1) * J_tot[i] - (L_tot[i] + 1) * L_tot[i] - (S_tot[i] + 1) * S_tot[i] ) ) # coef infront B
        param_y.append((I_tot[i] + 1) * I_tot[i]) # coef infront E
        param_z.append(SU_tot[i])               # coef infront G
        mass_sum.append(sum_mass[i])            # sum of masses

    param_v = np.array(param_v)
    param_w = np.array(param_w)
    param_x = np.array(param_x)
    param_y = np.array(param_y)
    param_z = np.array(param_z)
    mass_sum= np.array(mass_sum)
    
    return param_v,param_w,param_x,param_y,param_z,mass_sum


def reduced_mass(state, ModEx, HO_n):    
    # calculate parameters for all the states
    param_v = 0
    param_q1 = 1
    param_q2,param_q3 = 0,0
    param_is_rho,param_is_lam = 0,0
    param_is_omega,param_is_cascade,param_is_sigma=0,0,0
    const = np.sqrt(3.)
    
    # set HO excitation mode
    if(ModEx=='lam' or ModEx=='rpl'):
        param_is_lam = 1
        param_v      = HO_n * const
    elif(ModEx=='rho' or ModEx=='rpr'):
        param_is_rho = 1
        param_v      = HO_n * const
    elif(ModEx=='mix'):
        param_is_lam = 1
        param_is_rho = 1
        param_v      = const

    # set baryon type
    if(state=='omg'):
        param_q2 = 2
        param_q3 = 0
        param_is_omega = 1
    elif(state=='cas' or state=='tri'):
        param_q2 = 1
        param_q3 = 1
        param_is_cascade = 1
    elif(state=='sig' or state=='Lam'):
        param_q2 = 0
        param_q3 = 2
        param_is_sigma = 1
    
    return param_v,param_q1,param_q2,param_q3, param_is_rho,param_is_lam,\
        param_is_omega, param_is_cascade, param_is_sigma


def hamiltonian_quantum_factors_extended(state, J_tot, S_tot, L_tot, I_tot, SU_tot, HO_n, ModEx):
    """
    This method computes the (quantum) coeficients that multiplied the fitted/bootstrap parameters
    This function needs to be fed directly with quantum numbers stored in a numpy.array
    """
    param_v,param_w,param_x,param_y,param_z = [],[],[],[],[]
    param_q1,param_q2,param_q3 = [],[],[]
    param_is_rho,param_is_lam = [],[]
    param_is_omega,param_is_cascade,param_is_sigma = [],[],[]
    
    for i in range(len(J_tot)):

        param_v_i,\
            param_q1_i,param_q2_i, param_q3_i,\
            param_is_rho_i, param_is_lam_i, param_is_omega_i,\
            param_is_cascade_i,param_is_sigma_i = reduced_mass(state[i], ModEx[i], HO_n[i])
                
        param_v.append(param_v_i)    # coef infront kprim
        param_q1.append(param_q1_i)
        param_q2.append(param_q2_i) 
        param_q3.append(param_q3_i)
        param_is_rho.append(param_is_rho_i)
        param_is_lam.append(param_is_lam_i)     
        param_is_omega.append(param_is_omega_i)
        param_is_cascade.append(param_is_cascade_i)
        param_is_sigma.append(param_is_sigma_i)

        param_w.append((S_tot[i] + 1) * S_tot[i]) # coef infront A
        param_x.append(0.5 * ( (J_tot[i] + 1) * J_tot[i] - (L_tot[i] + 1) * L_tot[i] - (S_tot[i] + 1) * S_tot[i] ) ) # coef infront B
        param_y.append((I_tot[i] + 1) * I_tot[i]) # coef infront E
        param_z.append(SU_tot[i])               # coef infront G


    param_v = np.array(param_v)

    param_q1         = np.array(param_q1)
    param_q2         = np.array(param_q2)         
    param_q3         = np.array(param_q3)         
    param_is_rho     = np.array(param_is_rho)
    param_is_lam     = np.array(param_is_lam)     
    param_is_omega   = np.array(param_is_omega)
    param_is_cascade = np.array(param_is_cascade)
    param_is_sigma   = np.array(param_is_sigma)
    
    param_w = np.array(param_w)
    param_x = np.array(param_x)
    param_y = np.array(param_y)
    param_z = np.array(param_z)
    
    return param_v, param_q1, param_q2, param_q3, param_is_rho, param_is_lam,\
        param_is_omega, param_is_cascade, param_is_sigma,\
        param_w,param_x,param_y,param_z


def reduced_flavor(state, ModEx, HO_n):
    # calculate parameters for all the states
    param_v = 0
    param_q1 = 1
    param_q2,param_q3 = 0,0
    param_is_rho,param_is_lam = 0,0
    param_is_omega,param_is_cascade_p,param_is_sigma,param_is_lambda,param_is_cascade=0,0,0,0,0
    const = np.sqrt(3.)
    
    # set HO excitation mode
    if(ModEx=='lam' or ModEx=='rpl'):
        param_is_lam = 1
        param_v      = HO_n * const
    elif(ModEx=='rho' or ModEx=='rpr'):
        param_is_rho = 1
        param_v      = HO_n * const
    elif(ModEx=='mix'):
        param_is_lam = 1
        param_is_rho = 1
        param_v      = const

    # set baryon type
    if(state=='omg'):
        param_q2 = 2
        param_q3 = 0
        param_is_omega = 1
    elif(state=='cas'):
        param_q2 = 1
        param_q3 = 1
        param_is_cascade_p = 1
    elif(state=='sig'):
        param_q2 = 0
        param_q3 = 2
        param_is_sigma = 1
    elif(state=='Lam'):
        param_q2 = 0
        param_q3 = 2
        param_is_lambda = 1
    elif(state=='tri'):
        param_q2 = 1
        param_q3 = 1
        param_is_cascade = 1
    
    return param_v,param_q1,param_q2,param_q3, param_is_rho,param_is_lam,\
        param_is_omega, param_is_cascade_p, param_is_sigma, param_is_lambda, param_is_cascade

def hamiltonian_quantum_factors_flavor(state, J_tot, S_tot, L_tot, I_tot, SU_tot, HO_n, ModEx):
    """
    This method computes the (quantum) coeficients that multiplied the fitted/bootstrap parameters
    This function needs to be fed directly with quantum numbers stored in a numpy.array
    """
    param_v,param_w,param_x,param_y,param_z = [],[],[],[],[]
    param_q1,param_q2,param_q3 = [],[],[]
    param_is_rho,param_is_lam = [],[]
    param_is_omega,param_is_cascade_p,param_is_sigma,param_is_lambda,param_is_cascade = [],[],[],[],[]
    
    for i in range(len(J_tot)):

        param_v_i,\
            param_q1_i,param_q2_i, param_q3_i,\
            param_is_rho_i, param_is_lam_i, param_is_omega_i,\
            param_is_cascade_p_i, param_is_sigma_i, param_is_lambda_i, param_is_cascade_i = reduced_flavor(state[i], ModEx[i], HO_n[i])
                
        param_v.append(param_v_i)    # coef infront kprim
        param_q1.append(param_q1_i)
        param_q2.append(param_q2_i)         
        param_q3.append(param_q3_i)         
        param_is_rho.append(param_is_rho_i)
        param_is_lam.append(param_is_lam_i)     
        param_is_omega.append(param_is_omega_i)
        param_is_cascade_p.append(param_is_cascade_p_i)
        param_is_sigma.append(param_is_sigma_i)
        param_is_lambda.append(param_is_lambda_i)
        param_is_cascade.append(param_is_cascade_i)

        param_w.append((S_tot[i] + 1) * S_tot[i]) # coef infront A
        param_x.append(0.5 * ( (J_tot[i] + 1) * J_tot[i] - (L_tot[i] + 1) * L_tot[i] - (S_tot[i] + 1) * S_tot[i] ) ) # coef infront B
        param_y.append((I_tot[i] + 1) * I_tot[i]) # coef infront E
        param_z.append(SU_tot[i])               # coef infront G


    param_v = np.array(param_v)

    param_q1         = np.array(param_q1)
    param_q2         = np.array(param_q2)         
    param_q3         = np.array(param_q3)         
    param_is_rho     = np.array(param_is_rho)
    param_is_lam     = np.array(param_is_lam)     
    param_is_omega   = np.array(param_is_omega)
    param_is_cascade = np.array(param_is_cascade)
    param_is_sigma   = np.array(param_is_sigma)
    
    param_w = np.array(param_w)
    param_x = np.array(param_x)
    param_y = np.array(param_y)
    param_z = np.array(param_z)
    
    return param_v, param_q1, param_q2, param_q3, param_is_rho, param_is_lam,\
        param_is_omega, param_is_cascade_p, param_is_sigma, param_is_lambda, param_is_cascade,\
        param_w,param_x,param_y,param_z

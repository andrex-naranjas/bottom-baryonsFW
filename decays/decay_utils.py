#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
---------------------------------------------------------------
 Code to calcualte heavy-baryon decay widths
 authors: A. Ramirez-Morales (andres.ramirez.morales@cern.ch)
          H. Garcia-Tecocoatzi
 ---------------------------------------------------------------
'''

import numpy as np

pion_mass        = 0.13725 * 1000
kaon_mass        = 0.49564 * 1000
eta_mass         = 0.54786 * 1000
eta_p_mass       = 0.95778 * 1000
meson_rho_mass   = 0.77518 * 1000
kaon_s_mass      = 0.89555 * 1000
small_omega_mass = 0.78266 * 1000
phi_mass         = 1.01946 * 1000        
proton_mass      = 0.93891 * 1000
N_star_1520_mass = 1.51500 * 1000
N_star_1535_mass = 1.53000 * 1000
N_star_1680_mass = 1.68500 * 1000
N_star_1720_mass = 1.72000 * 1000
Delta_mass       = 1.23200 * 1000
Lambda_8_mass    = 1.11568 * 1000
Lambda_star_mass = 1.51900 * 1000
Xi_8_mass        = 1.31820 * 1000
Xi_10_mass       = 1.53370 * 1000
Sigma_8_mass     = 1.11932 * 1000
Sigma_10_mass    = 1.38460 * 1000
# bottom hadrons
lambda_mass      = 5.61960 * 1000
xi_p_mass        = 5.93502 * 1000
xi_mass          = 5.79700 * 1000
xi_s_mass        = 6.07800 * 1000
sigma_mass       = 5.81056 * 1000
sigma_s_mass     = 5.83032 * 1000
omega_mass       = 6.04520 * 1000
omega_s_mass     = 6.09300 * 1000
B0_mass          = 5.27966 * 1000
Bs_mass          = 5.36692 * 1000
B_star_mass      = 5.32471 * 1000


def append_dic(baryons, state1, state2,state3,
               state4, state5, state6, state7):
    baryons.append(state1)
    baryons.append(state2)
    baryons.append(state3)
    baryons.append(state4)
    baryons.append(state5)
    baryons.append(state6)
    baryons.append(state7)
    return baryons

def state_labels(baryon, ModEx, decPr, L_tot):
    
    if(baryon==1):
        baryon_name = "omega"
        if(decPr==1):
            decPr_name = "Xi+K"
        elif(decPr==2):
            decPr_name = "Xi'+K"
        elif(decPr==3):
            decPr_name = "Xi*+K"
        elif(decPr==4):
            decPr_name = "Xi+K*"
        elif(decPr==5):
            decPr_name = "Xi'+K*"
        elif(decPr==6):
            decPr_name = "Xi*'+K*"
        elif(decPr==7):
            decPr_name = "Omg+eta"
        elif(decPr==8):
            decPr_name = "Omg*+eta"
        elif(decPr==9):
            decPr_name = "Omg+phi"
        elif(decPr==10):
            decPr_name = "Omg*+phi"
        elif(decPr==11):
            decPr_name = "Omg+eta'"
        elif(decPr==12):
            decPr_name = "Omg*+eta'"
        elif(decPr==13):
            decPr_name = "Xi8+B"
        elif(decPr==14):
            decPr_name = "Xi10+B"

    elif(baryon==2 or baryon==5):
        if(baryon==2): baryon_name = 'cas_6'
        if(baryon==5): baryon_name = 'cas_3'
        if(decPr==1):
            decPr_name = "lam+K"
        elif(decPr==2):
            decPr_name = "Xi+Pi"
        elif(decPr==3):
            decPr_name = "Xi'+Pi"
        elif(decPr==4):
            decPr_name = "Xi*+Pi"
        elif(decPr==5):
            decPr_name = "Sig+K"
        elif(decPr==6):
            decPr_name = "Sig*+K"
        elif(decPr==7):
            decPr_name = "Xi+eta"
        elif(decPr==8):
            decPr_name = "Lam+K*"
        elif(decPr==9):
            decPr_name = "Xi+rho"
        elif(decPr==10):
            decPr_name = "Xi'+rho"
        elif(decPr==11):
            decPr_name = "Xi*+rho"
        elif(decPr==12):
            decPr_name = "Sig+K*"
        elif(decPr==13):
            decPr_name = "Sig*+K*"
        elif(decPr==14):
            decPr_name = "Xi'+eta"
        elif(decPr==15):
            decPr_name = "Xi*+eta"
        elif(decPr==16):
            decPr_name = "Xi+eta'"
        elif(decPr==17):
            decPr_name = "Xi'+eta'"
        elif(decPr==18):
            decPr_name = "Xi*+eta'"
        elif(decPr==19):
            decPr_name = "Xi+omg"
        elif(decPr==20):
            decPr_name = "Xi'+omg"
        elif(decPr==21):
            decPr_name = "Xi*+omg"
        elif(decPr==22):
            decPr_name = "Xi+phi"
        elif(decPr==23):
            decPr_name = "Xi'+phi"
        elif(decPr==24):
            decPr_name = "Xi*+phi"
        elif(decPr==25):
            if(baryon==2): decPr_name = "Sigma_8+B"
            if(baryon==5): decPr_name = "Lambda_8+B"
        elif(decPr==26):
            if(baryon==2): decPr_name = "Xi_8+Bs"
            if(baryon==5): decPr_name = "Lambda_8+B*"
        elif(decPr==27):
            if(baryon==2): decPr_name = "Sigma_8+B*"
            if(baryon==5): decPr_name = "Sigma_8+B"
        elif(decPr==28):
            if(baryon==2): decPr_name = "Sigma_10+B"
            if(baryon==5): decPr_name = "Lambda_8*+B"

    elif(baryon==3):
        baryon_name = "sigma"
        if(decPr==1):
            decPr_name = "Sig+Pi"
        elif(decPr==2):
            decPr_name = "Sig*+Pi"
        elif(decPr==3):
            decPr_name = "lam+Pi"
        elif(decPr==4):
            decPr_name = "Sig+Eta"
        elif(decPr==5):
            decPr_name = "Xi+K"
        elif(decPr==6):
            decPr_name = "Sig+rho"
        elif(decPr==7):
            decPr_name = "Sig*+rho"
        elif(decPr==8):
            decPr_name = "Lam+rho"
        elif(decPr==9):
            decPr_name = "Sig*+eta"
        elif(decPr==10):
            decPr_name = "Sig+eta'"
        elif(decPr==11):
            decPr_name = "Sig*+eta'"
        elif(decPr==12):
            decPr_name = "Xi'+K"
        elif(decPr==13):
            decPr_name = "Xi*+K"
        elif(decPr==14):
            decPr_name = "Xi+K*"
        elif(decPr==15):
            decPr_name = "Xi'+K*"
        elif(decPr==16):
            decPr_name = "Xi*+K*"
        elif(decPr==17):
            decPr_name = "Sig+omg"
        elif(decPr==18):
            decPr_name = "Sig*+omg"
        elif(decPr==19):
            decPr_name = "N+B"
        elif(decPr==20):
            decPr_name = "Sigma_8+Bs"
        elif(decPr==21):
            decPr_name = "N+B*"
        elif(decPr==22):
            decPr_name = "Delta+B"
        elif(decPr==23):
            decPr_name = "N*(1520)+B"
        elif(decPr==24):
            decPr_name = "N*(1535)+B"
        elif(decPr==25):
            decPr_name = "N*(1680)+B"
        elif(decPr==26):
            decPr_name = "N*(1720)+B"
            
    elif(baryon==4):
        baryon_name = 'lamda'
        if(decPr==1):
            decPr_name = "Sig+Pi"
        elif(decPr==2):
            decPr_name = "Sig*+Pi"
        elif(decPr==3):
            decPr_name = "lam+eta"
        elif(decPr==4):
            decPr_name = "Sig+rho"
        elif(decPr==5):
            decPr_name = "Sig*+rho"
        elif(decPr==6):
            decPr_name = "Lamb+eta'"
        elif(decPr==7):
            decPr_name = "Lamb+omg"
        elif(decPr==8):
            decPr_name = "Xi+K"
        elif(decPr==9):
            decPr_name = "Xi'+K"
        elif(decPr==10):
            decPr_name = "Xi*+K"
        elif(decPr==11):
            decPr_name = "Xi+K*"
        elif(decPr==12):
            decPr_name = "Xi'+K*"
        elif(decPr==13):
            decPr_name = "Xi*+K*"
        elif(decPr==14):
            decPr_name = "N+B"
        elif(decPr==15):
            decPr_name = "N+B*"

    if(ModEx==0):   ModEx_name ='Ground'
    elif(ModEx==1):
        if(L_tot==1):
            ModEx_name ='P-Wave-lam'
        elif(L_tot==2):
            ModEx_name ='D-Wave-lam'
    elif(ModEx==2):
        if(L_tot==1):
            ModEx_name ='P-Wave-rho'
        elif(L_tot==2):
            ModEx_name ='D-Wave-rho'
    elif(ModEx==3): ModEx_name ='Radial-lam'
    elif(ModEx==4): ModEx_name ='Radial-rho'
    elif(ModEx==5): ModEx_name ='Mixed'

    return baryon_name,ModEx_name,decPr_name


def asymmetric_decay_indi_error(list_array_decays):

    indi_up, indi_dn = ([]),([])    
    for i in range(len(list_array_decays[0])): # range == no.decays of a specific channel
        indi_channel = list_array_decays[:,i]
        indi_mean = np.mean(indi_channel)
        boot_size = indi_channel.size # array size == bootstrap size
        sort_indi_channel = np.sort(indi_channel)
        quantile_dn = int(boot_size*0.1587)   #int(np.floor(N*0.1587))
        quantile_up = int(boot_size*0.8413)+1 #int(np.floor(N*0.8413))
        indi_up = np.append(indi_up, sort_indi_channel[quantile_up-1] - indi_mean)
        indi_dn = np.append(indi_dn, sort_indi_channel[quantile_dn-1] - indi_mean)

    return indi_up, indi_dn
        

def latex_decay_label(baryon, decPr):

    if(baryon==1 or baryon=='omegas'):
        baryon_name = "omega"
        if(decPr==1):
            decPr_name = ("$\Xi_{b} K$", xi_mass, kaon_mass)
        elif(decPr==2):                                                                         
            decPr_name = ("$\Xi'_{b} K$", xi_p_mass, kaon_mass)
        elif(decPr==3):                              
            decPr_name = ("$\Xi^{*}_{b} K$", xi_s_mass, kaon_mass)
        elif(decPr==4):
            decPr_name = ("$\Xi_{b} K^{*}$", xi_mass, kaon_s_mass)      
        elif(decPr==5):                                                                         
            decPr_name = ("$\Xi'_{b}K^{*}$", xi_p_mass, kaon_s_mass)    
        elif(decPr==6):                                                                         
            decPr_name = ("$\Xi^{*}_{b} K^{*}$", xi_s_mass, kaon_s_mass)    
        elif(decPr==7):                                                                         
            decPr_name = ("$\Omega_{b} \eta$", omega_mass, eta_mass)      
        elif(decPr==8):                                                                         
            decPr_name = ("$\Omega^{*}_{b} \eta$", omega_s_mass, eta_mass)    
        elif(decPr==9):                                                                         
            decPr_name = ("$\Omega_{b} \phi$", omega_mass, phi_mass)      
        elif(decPr==10):
            decPr_name = ("$\Omega^{*}_{b} \phi$", omega_s_mass, phi_mass)    
        elif(decPr==11):
            decPr_name = ("$\Omega_{b} \eta'$", omega_mass, eta_p_mass)    
        elif(decPr==12):
            decPr_name = ("$\Omega^{*}_{b} \eta'$", omega_s_mass, eta_p_mass)  
        elif(decPr==13):
            decPr_name = ("$\Xi_{8} B$", Xi_8_mass, B0_mass)        
        elif(decPr==14):
            decPr_name = ("$\Xi_{10} B$", Xi_10_mass, B0_mass)

    elif(baryon==2 or baryon==5 or baryon=='cascades' or baryon=='cascades_anti3'):
        if(baryon==2): baryon_name = 'cas_6'
        if(baryon==5): baryon_name = 'cas_3'        
        if(decPr==1):
            decPr_name = ("$\Lambda_{b} K$", lambda_mass, kaon_mass)
        elif(decPr==2):             
            decPr_name = ("$\Xi_{b} \pi$", xi_mass, pion_mass)
        elif(decPr==3):                                                                                           
            decPr_name = ("$\Xi'_{b} \pi$", xi_p_mass, pion_mass)
        elif(decPr==4):                                                                                           
            decPr_name = ("$\Xi^{*}_{b} \pi$", xi_s_mass, pion_mass)
        elif(decPr==5):                                                                                           
            decPr_name = ("$\Sigma_{b} K$", sigma_mass, kaon_mass)
        elif(decPr==6):                                                                                           
            decPr_name = ("$\Sigma^{*}_{b} K$", sigma_s_mass, kaon_mass)
        elif(decPr==7):                                                                                           
            decPr_name = ("$\Xi_{b} \eta$", xi_mass, eta_mass)
        elif(decPr==8):                                                                                           
            decPr_name = ("$\Lambda_{b} K^{*}$", lambda_mass, kaon_s_mass)
        elif(decPr==9):                                                                                           
            decPr_name = ("$\Xi_{b} \\rho$", xi_mass, meson_rho_mass)
        elif(decPr==10):                                                                                          
            decPr_name = ("$\Xi'_{b} \\rho$", xi_p_mass, meson_rho_mass)
        elif(decPr==11):                                                                                          
            decPr_name = ("$\Xi^{*}_{b} \\rho$", xi_s_mass, meson_rho_mass)
        elif(decPr==12):                                                                                          
            decPr_name = ("$\Sigma_{b} K^{*}$", sigma_mass, kaon_s_mass)
        elif(decPr==13):                          
            decPr_name = ("$\Sigma^{*}_{b} K^{*}$", sigma_s_mass, kaon_s_mass)
        elif(decPr==14):                        
            decPr_name = ("$\Xi'_{b} \eta$", xi_p_mass, eta_mass)
        elif(decPr==15):
            decPr_name = ("$\Xi^{*}_{b} \eta$", xi_s_mass, eta_mass)
        elif(decPr==16):
            decPr_name = ("$\Xi_{b} \eta'$", xi_mass, eta_p_mass)
        elif(decPr==17):
            decPr_name = ("$\Xi'_{b} \eta'$", xi_p_mass, eta_p_mass)                     
        elif(decPr==18):
            decPr_name = ("$\Xi^{*}_{b} \eta'$", xi_s_mass, eta_p_mass)                     
        elif(decPr==19):
            decPr_name = ("$\Xi_{b} \omega$", xi_mass, small_omega_mass)
        elif(decPr==20):                                          
            decPr_name = ("$\Xi'_{b} \omega$", xi_p_mass, small_omega_mass)               
        elif(decPr==21):
            decPr_name = ("$\Xi^{*}_{b} \omega$", xi_s_mass, small_omega_mass)               
        elif(decPr==22):                                                                                          
            decPr_name = ("$\Xi_{b} \phi$", xi_mass, phi_mass)                         
        elif(decPr==23):
            decPr_name = ("$\Xi'_{b} \phi$", xi_p_mass, phi_mass)                       
        elif(decPr==24):
            decPr_name = ("$\Xi^{*}_{b} \phi$", xi_s_mass, phi_mass)
        elif(decPr==25 and baryon=='cascades'):
            decPr_name = ("$\Sigma_{8} B$", Sigma_8_mass, B0_mass) #//Sigma_8_lam+B   ^{\lambda}
        elif(decPr==25 and baryon=='cascades_anti3'):
            decPr_name = ("$\Lambda_{8} B$", Lambda_8_mass, B0_mass) # //Lambda_8_rho D    # ^{\\rho}
        elif(decPr==26 and baryon=='cascades'):
            decPr_name = ("$\Xi_{8} D_{s}$", Xi_8_mass, Bs_mass) # //Xi_8_lam+Bs      ^{\lambda}     
        elif(decPr==26 and baryon=='cascades_anti3'):
            decPr_name = ("$\Lambda_{8} D^{*}$", Lambda_8_mass, B_star_mass) # //Lambda_8_rho D # ^{\\rho}
        elif(decPr==27 and baryon=='cascades'):          
            decPr_name = ("$\Sigma_{8} D^{*}$", Sigma_8_mass, B_star_mass) #  //Sigma_8_lam+B* # ^{\lambda}
        elif(decPr==27 and baryon=='cascades_anti3'):
            decPr_name = ("$\Sigma_{8} B$", Sigma_8_mass, B0_mass) #  //Sigma_8_rho D    # ^{\\rho}
        elif(decPr==28 and baryon=='cascades'):
            decPr_name = ("$\Sigma_{10} B$", Sigma_10_mass, B0_mass) #  //Sigma_10_lam+B   #^{\lambda}
        elif(decPr==28 and baryon=='cascades_anti3'):
            decPr_name = ("$\Lambda_{8}^{*} B$",  Lambda_star_mass, B0_mass) # //Lambda*_8_rho D #,\\rho
            
    elif(baryon==3 or baryon=='sigmas'):
        baryon_name = 'sigma'
        if(decPr==1):
            decPr_name = ("$\Sigma_{b} \pi$",  sigma_mass, pion_mass)
        elif(decPr==2):
            decPr_name = ("$\Sigma^{*}_{b} \pi$", sigma_s_mass, pion_mass)                 
        elif(decPr==3):
            decPr_name = ("$\Lambda_{b} \pi$", lambda_mass, pion_mass)                 
        elif(decPr==4):
            decPr_name = ("$\Sigma_{b} \eta$", sigma_mass, eta_mass)                  
        elif(decPr==5):
            decPr_name = ("$\Xi_{b} K$", xi_mass, kaon_mass)
        elif(decPr==6):
            decPr_name = ("$\Sigma_{b}\\rho$",  sigma_mass, meson_rho_mass)         
        elif(decPr==7):                                                                             
            decPr_name = ("$\Sigma^{*}_{b}\\rho$", sigma_s_mass, meson_rho_mass)       
        elif(decPr==8):                                                                             
            decPr_name = ("$\Lambda_{b}\\rho$", lambda_mass, meson_rho_mass)        
        elif(decPr==9):                                                                             
            decPr_name = ("$\Sigma^{*}_{b}\eta$", sigma_s_mass, eta_mass)             
        elif(decPr==10):                                                                            
            decPr_name = ("$\Sigma_{b}\eta'$", sigma_mass, eta_p_mass)             
        elif(decPr==11):                                                                            
            decPr_name = ("$\Sigma^{*}_{b}\eta'$",  sigma_s_mass, eta_p_mass)           
        elif(decPr==12):                                                                            
            decPr_name = ("$\Xi'_{b}K$", xi_p_mass, kaon_mass)               
        elif(decPr==13):                                                                            
            decPr_name = ("$\Xi^{*}_{b}K$", xi_s_mass, kaon_mass)               
        elif(decPr==14):                                                                            
            decPr_name = ("$\Xi_{b} K^{*}$", xi_mass, kaon_s_mass)               
        elif(decPr==15):                                                                            
            decPr_name = ("$\Xi'_{b} K^{*}$",  xi_p_mass, kaon_s_mass)             
        elif(decPr==16):                                                                            
            decPr_name = ("$\Xi^{*}_{b} K^{*}$",  xi_s_mass, kaon_s_mass)             
        elif(decPr==17):                                                                            
            decPr_name = ("$\Sigma_{b}\omega$", sigma_mass, small_omega_mass)       
        elif(decPr==18):                                                                            
            decPr_name = ("$\Sigma^{*}_{b}\omega$", sigma_s_mass, small_omega_mass)     
        elif(decPr==19):                                                                            
            decPr_name = ("$N B$", proton_mass, B0_mass) #  //N+B                  #^{\lambda}
        elif(decPr==20):                                                                            
            decPr_name = ("$\Sigma_{8} B_{s}$", Sigma_8_mass, Bs_mass) # //Sigma_8+Bs           # ^{\lambda}
        elif(decPr==21):                                                                            
            decPr_name = ("$N B^{*}$", proton_mass, B_star_mass) # //N+B*                #^{\lambda}
        elif(decPr==22):                                                                            
            decPr_name = ("$\Delta B$", Delta_mass, B0_mass) #  //Delta+B        
        elif(decPr==23):                                                                            
            decPr_name = ("$N^{*}_{1} B$", N_star_1520_mass, B0_mass) # //N*(1520)+B     # (1520,\lambda
        elif(decPr==24):                                                                              
            decPr_name = ("$N^{*}_{2} B$", N_star_1535_mass, B0_mass) #  //N*(1535)+     # (1535,\lambda
        elif(decPr==25):                                                                              
            decPr_name = ("$N^{*}_{3} B$", N_star_1680_mass, B0_mass) # //N*(1680)+B      # (1680,\lambda
        elif(decPr==26):                                                                             
            decPr_name = ("$N^{*}_{4} B$", N_star_1720_mass, B0_mass) # //N*(1720)+B     # (1720,\lambda)
            
    elif(baryon==4 or baryon=='lambdas'):
        baryon_name = 'lamda'
        if(decPr==1):
            decPr_name = ("$\Sigma_{b} \pi$", sigma_mass, pion_mass)        
        elif(decPr==2):                                                                   
            decPr_name = ("$\Sigma^{*}_{b} \pi$", sigma_s_mass, pion_mass)        
        elif(decPr==3):                                                                   
            decPr_name = ("$\Lambda_{b} \eta$", lambda_mass, eta_mass)         
        elif(decPr==4):                                                                   
            decPr_name = ("$\Sigma_{b}\\rho$", sigma_mass, meson_rho_mass)    
        elif(decPr==5):                                                                   
            decPr_name = ("$\Sigma^{*}\\rho$", sigma_s_mass, meson_rho_mass)  
        elif(decPr==6):                                                                   
            decPr_name = ("$\Lambda_{b}\eta'$", lambda_mass, eta_p_mass)       
        elif(decPr==7):                                                                   
            decPr_name = ("$\Lambda_{b}\omega$", lambda_mass, small_omega_mass) 
        elif(decPr==8):                                                                   
            decPr_name = ("$\Xi_{b} K$", xi_mass, kaon_mass)            
        elif(decPr==9):                                                                   
            decPr_name = ("$\Xi'_{b} K$", xi_p_mass, kaon_mass)          
        elif(decPr==10):                                                                  
            decPr_name = ("$\Xi^{*}_{b} K$", xi_s_mass, kaon_mass)          
        elif(decPr==11):                                                                  
            decPr_name = ("$\Xi_{b} K^{*}$", xi_mass, kaon_s_mass)          
        elif(decPr==12):                                                                  
            decPr_name = ("$\Xi'_{b} K^{*}$", xi_p_mass, kaon_s_mass)        
        elif(decPr==13):
            decPr_name = ("$\Xi^{*}_{b} K^{*}$", xi_s_mass, kaon_s_mass)        
        elif(decPr==14):
            decPr_name = ("$N B$", proton_mass, B0_mass) #  //N+B  
        elif(decPr==15):
            decPr_name = ("$N B^{*}$", proton_mass, B0_mass) #  //N+B  not 15 decays
    return decPr_name


def print_row_latex(mass_a, masses_b, masses_c, state_name, state_decays, errors_up, errors_dn, f_out):
    nstate=len(state_decays)

    no_errors = False # for no bootstrap for decay widths
    if(errors_up is None or errors_dn is None):
        no_errors = True
        
    no_errors = True
    cons_energy_count = 0
    print(state_name, end='',file=f_out)
    for i in range(nstate):
        value=0
        if(state_decays[i]==0.0):
            if mass_a < masses_b[i] + masses_c[i]: # check if the width is zero for energy conservation
                cons_energy_count+=1
                value = '$0$'
                print(value,"  &", end='', file=f_out)
            else:
                value = '-'
                print(value,"  &", end='', file=f_out)
        else:
            value = round(state_decays[i],1)
            if not no_errors:
                error_up = abs(round(errors_up[i],1))
                error_dn = abs(round(errors_dn[i],1))
                print("$",value,"$  &  ", end='', file=f_out)
                # print("$",value,"_{-",error_dn, "}^{+",error_up,"}$  &  ", end='', file=f_out)
            else:
                print(value,"  &", end='', file=f_out)

    # total decay width           
    if np.sum(state_decays) != 0.0:
        sum_value = round(np.sum(state_decays),1)
        if not no_errors:
            error_sum_up = round(np.sqrt(np.sum(np.square(errors_up))),1)
            error_sum_dn = (-1.)*round(np.sqrt(np.sum(np.square(errors_dn))),1)
            print("$",sum_value,"_{",error_sum_dn, "}^{+",error_sum_up,"}$ \\\\", file=f_out)
        else:
            print(sum_value," \\\\", file=f_out)            
    else:
        if cons_energy_count != 0:
            sum_value = '$0$'
            print(sum_value," \\\\", file=f_out)
        else:
            sum_value = '-'
            print(sum_value," \\\\", file=f_out)    
    

def print_header_latex(name_states, f_out):
    nNames = len(name_states)    
    print("\\begin{tabular}{c |", end='',file=f_out)
    for i in range(nNames-2):
        print("  p{0.58cm}", end='',file=f_out)
        
    print("p{0.75cm}} \hline \hline", file=f_out)
    for i in range(nNames-1): print(name_states[i]," & ", end='',file=f_out)
    print(name_states[nNames-1]," \\\\ \hline", file=f_out)

def print_bottom_latex(baryons,f_decay):
    print('\hline \hline', file=f_decay)
    print('\end{tabular}', file=f_decay)
    #print("\caption{Decay widths in MeV, for states: $", baryons,"$}",file=f_decay)
    #print("\label{tab:gordo}", file=f_decay)


def decay_masses(baryons, decPr):
    # fetch mass of the decay products        
    pion_mass   = 0.140
    kaon_mass   = 0.493
    lambda_mass = 2.286
    xi_mass     = 2.469
    xi_p_mass   = 2.578
    xi_s_mass   = 2.645
    sigma_mass  = 2.455
    sigma_s_mass= 2.518 # 2.520
    eta_mass    = 0.548
    
    if(baryons=='omegas'):
        if(decPr==1):     return xi_mass,   kaon_mass
        elif(decPr==2):   return xi_p_mass, kaon_mass
        elif(decPr==3):   return xi_s_mass, kaon_mass
    elif(baryons=='cascades'):
        if(decPr==1):     return lambda_mass, kaon_mass
        elif(decPr==2):   return xi_mass,     pion_mass
        elif(decPr==3):   return xi_p_mass,   pion_mass
        elif(decPr==4):   return xi_s_mass,   pion_mass
        elif(decPr==5):   return sigma_mass,  kaon_mass
        elif(decPr==6):   return sigma_s_mass,kaon_mass
        elif(decPr==7):   return xi_mass,     eta_mass
    elif(baryons=='sigmas'):
        if(decPr==1):     return sigma_mass,   pion_mass
        elif(decPr==2):   return sigma_s_mass, pion_mass
        elif(decPr==3):   return lambda_mass,  pion_mass
        elif(decPr==4):   return sigma_mass,   eta_mass
        elif(decPr==5):   return xi_mass,      kaon_mass
    elif(baryons=='lambdas'):
        if(decPr==1):     return sigma_mass,   pion_mass
        elif(decPr==2):   return sigma_s_mass, pion_mass
        elif(decPr==3):   return lambda_mass,  eta_mass
    elif(baryons=='cascades_anti3'):
        if(decPr==1):     return lambda_mass, kaon_mass
        elif(decPr==2):   return xi_mass,     pion_mass
        elif(decPr==3):   return xi_p_mass,   pion_mass
        elif(decPr==4):   return xi_s_mass,   pion_mass
        elif(decPr==5):   return sigma_mass,  kaon_mass
        elif(decPr==6):   return sigma_s_mass,kaon_mass
        elif(decPr==7):   return xi_mass,     eta_mass

def baryon_symbol(baryons="omegas"):
    if baryons=="omegas":
        return "\Omega"
    elif baryons=="sigmas":
        return "\Sigma"
    elif baryons=="lambdas":
        return "\Lambda"
    elif baryons=="cascades":
        return "\Xi'"
    else:
        return "\Xi"

def baryon_quarks(baryons="omegas"):
    if baryons=="omegas":
        return "ssb"
    elif baryons=="sigmas":
        return "nnb"
    elif baryons=="lambdas":
        return "nnb"
    elif baryons=="cascades":
        return "snb'"
    else:
        return "snb"

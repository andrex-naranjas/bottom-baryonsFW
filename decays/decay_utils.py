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
            decPr_name = "Xi8+D"
        elif(decPr==14):
            decPr_name = "Xi10+D"

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
            if(baryon==2): decPr_name = "Sigma_8+D"
            if(baryon==5): decPr_name = "Lambda_8+D"
        elif(decPr==26):
            if(baryon==2): decPr_name = "Xi_8+Ds"
            if(baryon==5): decPr_name = "Lambda_8+D*"
        elif(decPr==27):
            if(baryon==2): decPr_name = "Sigma_8+D*"
            if(baryon==5): decPr_name = "Sigma_8+D"
        elif(decPr==28):
            if(baryon==2): decPr_name = "Sigma_10+D"
            if(baryon==5): decPr_name = "Lambda_8*+D"

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
            decPr_name = "Xi'+K'"
        elif(decPr==13):
            decPr_name = "Xi*+K'"
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
            decPr_name = "N+D"
        elif(decPr==20):
            decPr_name = "Sigma_8+Ds"
        elif(decPr==21):
            decPr_name = "N+D*"
        elif(decPr==22):
            decPr_name = "Delta+D"
        elif(decPr==23):
            decPr_name = "N*(1520)+D"
        elif(decPr==24):
            decPr_name = "N*(1535)+D"
        elif(decPr==25):
            decPr_name = "N*(1680)+D"
        elif(decPr==26):
            decPr_name = "N*(1720)+D"
            
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
            decPr_name = "N+D"
        elif(decPr==15):
            decPr_name = "N+D*"

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
        

def latex_decay_label(baryon,decPr):

    if(baryon==1 or baryon=='omegas'):
        baryon_name = "omega"
        if(decPr==1):
            decPr_name = "$\Xi_{c} K$"
        elif(decPr==2):
            decPr_name = "$\Xi'_{c} K$"
        elif(decPr==3):
            decPr_name = "$\Xi^{*}_{c} K$"
        elif(decPr==4):
            decPr_name = "$\Xi_{c} K^{*}$"
        elif(decPr==5):
            decPr_name = "$\Xi'_{c}K^{*}$"
        elif(decPr==6):
            decPr_name = "$\Xi^{*}_{c} K^{*}$"
        elif(decPr==7):
            decPr_name = "$\Omega_{c} \eta$"
        elif(decPr==8):
            decPr_name = "$\Omega^{*}_{c} \eta$"
        elif(decPr==9):
            decPr_name = "$\Omega_{c} \phi$"
        elif(decPr==10):
            decPr_name = "$\Omega^{*}_{c} \phi$"
        elif(decPr==11):
            decPr_name = "$\Omega_{c} \eta'$"
        elif(decPr==12):
            decPr_name = "$\Omega^{*}_{c} \eta'$"
        elif(decPr==13):
            decPr_name = "$\Xi_{8} D$"
        elif(decPr==14):
            decPr_name = "$\Xi_{10} D$"            

    elif(baryon==2 or baryon==5 or baryon=='cascades' or baryon=='cascades_anti3'):
        if(baryon==2): baryon_name = 'cas_6'
        if(baryon==5): baryon_name = 'cas_3'        
        if(decPr==1):
            decPr_name = "$\Lambda_{c} K$"
        elif(decPr==2):
            decPr_name = "$\Xi_{c} \pi$"
        elif(decPr==3):
            decPr_name = "$\Xi'_{c} \pi$"
        elif(decPr==4):
            decPr_name = "$\Xi^{*}_{c} \pi$"
        elif(decPr==5):
            decPr_name = "$\Sigma_{c} K$"
        elif(decPr==6):
            decPr_name = "$\Sigma^{*}_{c} K$"
        elif(decPr==7):
            decPr_name = "$\Xi_{c} \eta$"
        elif(decPr==8):
            decPr_name = "$\Lambda_{c} K^{*}$"
        elif(decPr==9):
            decPr_name = "$\Xi_{c} \\rho$"
        elif(decPr==10):
            decPr_name = "$\Xi'_{c} \\rho$"
        elif(decPr==11):
            decPr_name = "$\Xi^{*}_{c} \\rho$"
        elif(decPr==12):
            decPr_name = "$\Sigma_{c} K^{*}$"
        elif(decPr==13):
            decPr_name = "$\Sigma^{*}_{c} K^{*}$"
        elif(decPr==14):
            decPr_name = "$\Xi'_{c} \eta$"
        elif(decPr==15):
            decPr_name = "$\Xi^{*}_{c} \eta$"
        elif(decPr==16):
            decPr_name = "$\Xi_{c} \eta'$"
        elif(decPr==17):
            decPr_name = "$\Xi'_{c} \eta'$"
        elif(decPr==18):
            decPr_name = "$\Xi^{*}_{c} \eta'$"
        elif(decPr==19):
            decPr_name = "$\Xi_{c} \omega$"
        elif(decPr==20):
            decPr_name = "$\Xi'_{c} \omega$"
        elif(decPr==21):
            decPr_name = "$\Xi^{*}_{c} \omega$"
        elif(decPr==22):
            decPr_name = "$\Xi_{c} \phi$"
        elif(decPr==23):
            decPr_name = "$\Xi'_{c} \phi$"
        elif(decPr==24):
            decPr_name = "$\Xi^{*}_{c} \phi$"
        elif(decPr==25 and baryon=='cascades'):
            decPr_name = "$\Sigma_{8} D$" # ^{\lambda}
        elif(decPr==25 and baryon=='cascades_anti3'):
            decPr_name = "$\Lambda_{8} D$" # ^{\\rho}
        elif(decPr==26 and baryon=='cascades'):
            decPr_name = "$\Xi_{8} D_{s}$" # ^{\lambda}
        elif(decPr==26 and baryon=='cascades_anti3'):
            decPr_name = "$\Lambda_{8} D^{*}$" # ^{\\rho}
        elif(decPr==27 and baryon=='cascades'):
            decPr_name = "$\Sigma_{8} D^{*}$" # ^{\lambda}
        elif(decPr==27 and baryon=='cascades_anti3'):
            decPr_name = "$\Sigma_{8} D$" # ^{\\rho}
        elif(decPr==28 and baryon=='cascades'):
            decPr_name = "$\Sigma_{10} D$" #^{\lambda}
        elif(decPr==28 and baryon=='cascades_anti3'):
            decPr_name = "$\Lambda_{8}^{*} D$" #,\\rho
            
    elif(baryon==3 or baryon=='sigmas'):
        baryon_name = 'sigma'
        if(decPr==1):
            decPr_name = "$\Sigma_{c} \pi$"
        elif(decPr==2):
            decPr_name = "$\Sigma^{*}_{c} \pi$"
        elif(decPr==3):
            decPr_name = "$\Lambda_{c} \pi$"
        elif(decPr==4):
            decPr_name = "$\Sigma_{c} \eta$"
        elif(decPr==5):
            decPr_name = "$\Xi_{c} K$"
        elif(decPr==6):
            decPr_name = "$\Sigma_{c}\\rho$"
        elif(decPr==7):
            decPr_name = "$\Sigma^{*}_{c}\\rho$"
        elif(decPr==8):
            decPr_name = "$\Lambda_{c}\\rho$"
        elif(decPr==9):
            decPr_name = "$\Sigma^{*}_{c}\eta$"
        elif(decPr==10):
            decPr_name = "$\Sigma_{c}\eta'$"
        elif(decPr==11):
            decPr_name = "$\Sigma^{*}_{c}\eta'$"
        elif(decPr==12):
            decPr_name = "$\Xi'_{c}K'$"
        elif(decPr==13):
            decPr_name = "$\Xi^{*}_{c}K'$"
        elif(decPr==14):
            decPr_name = "$\Xi K^{*}$"
        elif(decPr==15):
            decPr_name = "$\Xi' K^{*}$"
        elif(decPr==16):
            decPr_name = "$\Xi^{*}_{c} K^{*}$"
        elif(decPr==17):
            decPr_name = "$\Sigma\omega$"
        elif(decPr==18):
            decPr_name = "$\Sigma^{*}_{c}\omega$"
        elif(decPr==19):
            decPr_name = "$N D$" #^{\lambda}
        elif(decPr==20):
            decPr_name = "$\Sigma_{8} D_{s}$" # ^{\lambda}
        elif(decPr==21):
            decPr_name = "$N D^{*}$" #^{\lambda}
        elif(decPr==22):
            decPr_name = "$\Delta D$"
        elif(decPr==23):
            decPr_name = "$N^{*}_{1} D$" # (1520),\lambda
        elif(decPr==24):                                 
            decPr_name = "$N^{*}_{2} D$" # (1535),\lambda
        elif(decPr==25):                                 
            decPr_name = "$N^{*}_{3} D$" # (1680),\lambda
        elif(decPr==26):                                 
            decPr_name = "$N^{*}_{4} D$"  # (1720,\lambda)
            
    elif(baryon==4 or baryon=='lambdas'):
        baryon_name = 'lamda'
        if(decPr==1):
            decPr_name = "$\Sigma_{c} \pi$"
        elif(decPr==2):
            decPr_name = "$\Sigma^{*}_{c} \pi$"
        elif(decPr==3):
            decPr_name = "$\Lambda_{c} \eta$"
        elif(decPr==4):
            decPr_name = "$\Sigma_{c}\\rho$"
        elif(decPr==5):
            decPr_name = "$\Sigma^{*}\\rho$"
        elif(decPr==6):
            decPr_name = "$\Lambda_{c}\eta'$"
        elif(decPr==7):
            decPr_name = "$\Lambda_{c}\omega$"
        elif(decPr==8):
            decPr_name = "$\Xi_{c} K$"
        elif(decPr==9):
            decPr_name = "$\Xi'_{c} K$"
        elif(decPr==10):
            decPr_name = "$\Xi^{*}_{c} K$"
        elif(decPr==11):
            decPr_name = "$\Xi_{c} K^{*}$"
        elif(decPr==12):
            decPr_name = "$\Xi'_{c} K^{*}$"
        elif(decPr==13):
            decPr_name = "$\Xi^{*}_{c} K^{*}$"
        elif(decPr==14):
            decPr_name = "$N D$"
        elif(decPr==15):
            decPr_name = "$N D^{*}$"

    return decPr_name

def print_row_latex(state_name, state_decays, errors_up, errors_dn, f_out):
    nstate=len(state_decays)

    no_errors = False # for no bootstrap for decay widths
    if(errors_up is None or errors_dn is None):
        no_errors = True
        
    no_errors = True
    print(state_name, end='',file=f_out)
    for i in range(nstate):
        value=0
        if(state_decays[i]==0.0):
            value = '$\\dagger$'
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

                
    if np.sum(state_decays) != 0.0:
        sum_value = round(np.sum(state_decays),1)
        if not no_errors:
            error_sum_up = round(np.sqrt(np.sum(np.square(errors_up))),1)
            error_sum_dn = (-1.)*round(np.sqrt(np.sum(np.square(errors_dn))),1)
            print("$",sum_value,"_{",error_sum_dn, "}^{+",error_sum_up,"}$ \\\\", file=f_out)
        else:
            print(sum_value," \\\\", file=f_out)            
    else:
        sum_value = '$\\dagger$'
        print(sum_value," \\\\", file=f_out)    
    

def print_header_latex(name_states, f_out):
    nNames = len(name_states)    
    print("\\begin{tabular}{c |", end='',file=f_out)
    for i in range(nNames-2):
        print("  p{0.35cm}", end='',file=f_out)
        
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

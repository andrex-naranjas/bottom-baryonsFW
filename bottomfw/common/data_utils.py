#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
---------------------------------------------------------------
 Authors: A. Ramirez-Morales (andres.ramirez.morales@cern.ch)
          H. Garcia-Tecocoatzi
---------------------------------------------------------------
"""
# data utils module
import numpy as np
from scipy.stats import norm, normaltest, shapiro, chisquare, kstest
from statsmodels.stats.diagnostic import lilliefors


def find_best_sigma_0():
    # function to find the best model error to obtain chi-2/NOF = 1
    # not usable as it is, need to be plugged in bottom_bootstrap and uncomment the sigma_0 lines in least_squares method
    exp_m = np.array([
        # omegas
        np.mean(gauss_2695), np.mean(gauss_2770), np.mean(gauss_3000),
        np.mean(gauss_3090),
        np.mean(gauss_3120),
        # cascades six-plet
        np.mean(gauss_2578), np.mean(gauss_2645),
        np.mean(gauss_2965),
        np.mean(gauss_3055),
        np.mean(gauss_3080),
        # sigma
        np.mean(gauss_2453), np.mean(gauss_2518), np.mean(gauss_2801),
        # lambda         
        np.mean(gauss_2286), np.mean(gauss_2592), np.mean(gauss_2628),
        # cascade anti-3plet
        np.mean(gauss_2469), np.mean(gauss_2792), np.mean(gauss_2816)
    ])

    for i in range(5000):
        # find by eye sigma_0
        y_errors_exp =np.array([1.7, 2.0, 0.22, 0.5, 1.0, 0.9, 0.7, 0.3, 0.4, 2.8, 0.9, 2.8, 20, 0.14, 0.28, 0.19, 1.77, 1.4, 2.07])
        # for diquark
        # y_errors_exp =np.array([1.7, 2.00, 0.22, 0.50, 1.00, 0.90, 0.71, 0.33,0.90, 2.80, 20.1, 0.14, 0.28, 0.19, 1.77, 1.40, 2.07])
        sigma_0 = i*0.01
        m = fit(least_squares)
        if(0.9 < m.fval/(19-m.nfit) < 1.1):
            print(m.fval/((19-m.nfit)), round(m.values['m1']), round(m.values['m2']), round(m.values['m3']), sigma_0)
            print(round(m.values['k']), round(m.values['a']),  round(m.values['b']),   round(m.values['e']), round(m.values['g']) )
            print('********************************************************')
    input()
    print("End scan!")


def linear_algebra_check(name, workpath='.'):
    #name='All'
    from oct2py import octave
    A = octave.load(workpath+'/octave/matrix_'+ name + '.mat')
    b = octave.load(workpath+'/octave/vector_'+ name + '.mat')
    return octave.linsolve(A,b)


def normal_test(sample,alpha,verbose):
    # hypothesis test: null hypothesis, the data is gaussian distributed

    # Shapiro-Wilk
    stat, p = shapiro(sample)
    if verbose:
        if p > alpha: print('Shapiro this is Gaussian', p)
        else:         print('Shapiro this is NOT Gaussian', p)

    # chisquare
    stat, p = chisquare(sample)
    if verbose:
        if p > alpha: print('Chisquare this is Gaussian', p)
        else:         print('Chisquare this is NOT Gaussian', p)

    # lilliefors
    stat, p = lilliefors(sample)
    if verbose:
        if p > alpha: print('Lilliefors this is Gaussian', p)
        else:         print('Lilliefors this is NOT Gaussian', p)

    # kolmogorov
    stat, p = kstest(sample, 'norm')
    if verbose:
        if p > alpha: print('Kolmogorov this is Gaussian', p)
        else:         print('Kolmogorov this is NOT Gaussian', p)

    # Angostino
    k2, p = normaltest(sample)
    if verbose:
        if p > alpha: print('Angostino this is Gaussian', p)
        else:         print('Angostino this is NOT Gaussian', p)    

    return p,alpha


def names_values(name):
    #if name=='All':
    quantum = [
    # omega
    '$\\vert ssc,1/2,1/2,0,0,10/3 \\rangle$ &',   '$\\vert ssc,3/2,3/2,0,0,10/3 \\rangle$ &',   '$\\vert ssc,1/2,1/2,1,0,10/3 \\rangle$ &',
    '$\\vert ssc,1/2,3/2,1,0,10/3 \\rangle$ &',   '$\\vert ssc,3/2,1/2,1,0,10/3 \\rangle$ &',   '$\\vert ssc,1/2,1/2,0,0,10/3 \\rangle$ &',
    '$\\vert ssc,5/2,3/2,0,1,10/3 \\rangle$ &',
    # cascades six-plet                   
    '$\\vert suc,1/2,1/2,0,1/2,10/3 \\rangle$ &', '$\\vert suc,3/2,3/2,0,1/2,10/3 \\rangle$ &', '$\\vert suc,1/2,3/2,1,1/2,10/3 \\rangle$ &',
    '$\\vert suc,3/2,1/2,1,1/2,10/3 \\rangle$ &', '$\\vert suc,3/2,3/2,1,1/2,10/3 \\rangle$ &', '$\\vert suc,1/2,1/2,1,1/2,10/3 \\rangle$ &',
    '$\\vert suc,3/2,1/2,1,1/2,10/3 \\rangle$ &',
    # sigma
    '$\\vert uuc,1/2,1/2,0,1,10/3 \\rangle$ &',   '$\\vert uuc,3/2,3/2,0,1,10/3 \\rangle$ &',   '$\\vert uuc,1/2,1/2,1,1,10/3 \\rangle$ &',
    # lambda
    '$\\vert udc,1/2,1/2,0,0,4/3  \\rangle$ &',   '$\\vert udc,1/2,1/2,1,0,10/3 \\rangle$ &',   '$\\vert udc,3/2,1/2,1,0,4/3  \\rangle$ &',
    # cascades 3plet
    '$\\vert ssc,1/2,1/2,0,1/2,4/3 \\rangle$ &',  '$\\vert ssc,1/2,1/2,1,1/2,4/3 \\rangle$ &',  '$\\vert ssc,3/2,1/2,1,1/2,10/3 \\rangle$ &',
    '$\\vert ssc,3/2,1/2,1,1/2,10/3 \\rangle$ &']
    
    old = np.array([
        # omegas
        2702.4, 2767.0, 3015.8, #3044.5,# 3051.6,
        3080.4,
        3140.2,
        # cascades six-plet
        2570.1, 2634.8, #2934.1, 2941.2,
        2969.9,
        3060.2,
        3096.07,
        # sigma
        2453.1, 2517.7, 2819.0,
        # lambda
        2283.7, 2649.7, 2685.6,
        # cascades triplet
        2461.2, 2796.5, 2832.4#, 2980.2
        ])
    exp = np.array([
        # omegas
        2695.0, 2765.9, 3000.4, #3050.2,# 3065.6,
        3090.2,
        3119.1,
        # cascades six-plet
        2578.0, 2645.97, #2923.1,#2938.6,#
        2964.9,
        3055.9,
        3078.6,
        # sigma            
        2453.5, 2518.1, 2800.0,
        # lambda C
        2286.5, 2592.3, 2628.1,
        # cascade anti-3-plet
        2469.42, 2793.3, 2818.5#, 2968.6
    ])
    delta_exp = np.array([
        # omegas
        1.7, 2.0, 0.22, 0.13,# 0.28, 0.5,
        1.0,
        # cascades six-plet
        0.9, 0.71, 0.4, 0.3,# 3.3, 0.4,
        2.8,
        # sigma            
        0.9,2.8,20,
        # lambda
        0.14, 0.28,0.19,
        # cascade 3 plet
        1.8,1.4,2.1#, 3.3
    ])
    
    return quantum, old, exp, delta_exp


def wave_label(S_tot, J_tot, L_tot):
    spin_label='^{'+str(int(2*S_tot+1))+'}'    
    jtot_string = ''
    if J_tot   == 0.5: jtot_string = '1/2'
    elif J_tot == 1.5: jtot_string = '3/2'
    elif J_tot == 2.5: jtot_string = '5/2'
    elif J_tot == 3.5: jtot_string = '7/2'            
    jtot_label='_{'+jtot_string+'}'
    if L_tot==0:
        return "$"+spin_label + 'S' + jtot_label + '$'
    elif L_tot==1:
        return "$"+spin_label + 'P' + jtot_label + '$'
    elif L_tot==2:
        return "$"+spin_label + 'D' + jtot_label + '$'

def flavor_label(baryons='omegas'):
    if baryons=='omegas' or baryons=='cascades' or baryons=='sigmas':
        return '${\\bf {6}}_{\\rm f}$' # 6_{sf}
    elif baryons=='lambdas' or baryons=='cascades_anti3':
        return '${\\bf {3}}_{\\rm c}$' # 3_{c}
    else:
        return '$x$'
    
def baryon_name(baryons):
    if baryons=='omegas':
        return "$\Omega_{b}(ssb)$"
    elif baryons=='sigmas':
        return "$\Sigma_{b}(nnb)$"
    elif baryons=='cascades':
        return "$\Xi'_{b}(snb)$"
    elif baryons=='cascades_anti3':
        return "$\Xi_{b}(snb)$"
    elif baryons=='lambdas':
        return "$\Lambda_{b}(snb)$"
    else:
        return ''

def baryon_name_plots(baryons):
    if baryons=='omegas':
        return "$\Omega_{b}$"
    elif baryons=='sigmas':
        return "$\Sigma_{b}$"
    elif baryons=='cascades':
        return "$\Xi'_{b}$"
    elif baryons=='cascades_anti3':
        return "$\Xi_{b}$"
    elif baryons=='lambdas':
        return "$\Lambda_{b}$"
    else:
        return ''

def name_quantum_state(baryons, J_tot, S_tot, L_tot, ModEx, SU_tot):

    if baryons=='omegas':
        quarks = 'ssb'
    elif baryons=='cascades' or baryons=='cascades_anti3':
        quarks = 'snb'
    elif baryons=='lambdas' or baryons=='sigmas':
        quarks = 'nnb'
    else:
        quarks = 'x'
        
    if J_tot==0.5:
        j_name = '1/2'
    elif J_tot==1.5:
        j_name = '3/2'
    elif J_tot==2.5:
        j_name = '5/2'
    elif J_tot==3.5:
        j_name = '7/2'
    else:
        j_name = 'x'

    if S_tot==0.5:
        s_name = '1/2'
    elif S_tot==1.5:
        s_name = '3/2'
    else:
        s_name = 'x'

    if L_tot==0 and (ModEx=='lam' or ModEx=='rho' or ModEx=='grd'):
        l_lam = '0'
        l_rho = '0'
    elif L_tot==1 and ModEx=='lam':
        l_lam = '1'
        l_rho = '0'
    elif L_tot==1 and ModEx=='rho':
        l_lam = '0'
        l_rho = '1'
    elif L_tot==2 and ModEx=='lam':
        l_lam = '2'
        l_rho = '0'
    elif L_tot==2 and ModEx=='rho':
        l_lam = '0'
        l_rho = '2'
    elif L_tot==0 and ModEx=='rpl':
        l_lam = '0'
        l_rho = '0'
    elif L_tot==0 and ModEx=='rpr':
        l_lam = '0'
        l_rho = '0'
    elif L_tot==2 and ModEx=='mix':
        l_lam = '1'
        l_rho = '1'
    elif L_tot==1 and ModEx=='mix':
        l_lam = '1'
        l_rho = '1'
    elif L_tot==0 and ModEx=='mix':
        l_lam = '1'
        l_rho = '1'
    else:
        l_lam = 'x'
        l_rho = 'x'
        
    kappa_lam, kappa_rho='0','0'
    if L_tot==0 and ModEx=='rpl':
        kappa_lam = '1' # '1_{\\lambda}'
    elif L_tot==0 and ModEx=='rpr':
        kappa_rho = '1' #'1_{\\rho}'

    if SU_tot==float(10/3):
        flavour = '{\\bf {6}}_{\\rm f}' # 6_{sf}
    elif SU_tot==float(4/3):
        flavour = '{\\bf {3}}_{\\rm c}' # 3_{c}
    else:
        flavour = 'x'

    return '$\\vert l_{\\lambda}\\!\\!=\\!'+ l_lam + ', l_{\\rho}\\!\\!=\\!'+l_rho + ', k_{\\lambda}\\!\\!=\\!' + kappa_lam + ', k_{\\rho}\\!\\!=\\!'+kappa_rho +' \\rangle$ &'
        
    # return '$\\vert '+ quarks +','+ j_name + ',' + s_name +',' + l_lam + ',' + l_rho + ',' + kappa + ',' + flavour +' \\rangle$ &'


def exp_mass_val(baryons, J_tot, S_tot, L_tot, ModEx, SU_tot):
    
    if baryons=='omegas':
        if J_tot==1/2 and S_tot==1/2 and L_tot==0 and ModEx=='grd' and SU_tot==float(10/3):
            return 6045.2, 1.0
        elif J_tot==1/2 and S_tot==1/2 and L_tot==1 and ModEx=='lam' and SU_tot==float(10/3):
            return 6315.6, 1.0
        elif J_tot==3/2 and S_tot==1/2 and L_tot==1 and ModEx=='lam' and SU_tot==float(10/3):
            return 6330.3, 1.0
        elif J_tot==1/2 and S_tot==3/2 and L_tot==1 and ModEx=='lam' and SU_tot==float(10/3):
            return 6339.7, 1.0
        elif J_tot==3/2 and S_tot==3/2 and L_tot==1 and ModEx=='lam' and SU_tot==float(10/3):
            return 6349.8, 1.0
        else:
            return 0,0
    elif baryons=='cascades':
        if J_tot==1/2 and S_tot==1/2 and L_tot==0 and ModEx=='grd' and SU_tot==float(10/3):
            return 5935.02, 1.0
        elif J_tot==3/2 and S_tot==3/2 and L_tot==0 and ModEx=='grd' and SU_tot==float(10/3):
            return 5953.82, 1.0
        elif J_tot==5/2 and S_tot==3/2 and L_tot==1 and ModEx=='lam' and SU_tot==float(10/3):
            return 6227.85, 1.0
        else:
            return 0, 0    
    if baryons=='sigmas':
        if J_tot==1/2 and S_tot==1/2 and L_tot==0 and ModEx=='grd' and SU_tot==float(10/3):
            return 5813.1, 1.0
        elif J_tot==3/2 and S_tot==3/2 and L_tot==0 and ModEx=='grd' and SU_tot==float(10/3):
            return 5837.03, 1.0
        elif J_tot==1/2 and S_tot==1/2 and L_tot==1 and ModEx=='lam' and SU_tot==float(10/3):
            return 6096.9, 1.0
        else:
            return 0, 0
    elif baryons=='lambdas':
        if J_tot==1/2 and S_tot==1/2 and L_tot==0 and ModEx=='grd' and SU_tot==float(4/3):
            return 5619.6, 1.0
        elif J_tot==1/2 and S_tot==1/2 and L_tot==1 and ModEx=='lam' and SU_tot==float(4/3):
            return 5912.2, 1.0
        elif J_tot==3/2 and S_tot==1/2 and L_tot==1 and ModEx=='lam' and SU_tot==float(4/3):
            return 5920.1, 1.0
        elif J_tot==3/2 and S_tot==1/2 and L_tot==2 and ModEx=='lam' and SU_tot==float(4/3):
            return 6146.2, 1.0
        elif J_tot==5/2 and S_tot==1/2 and L_tot==2 and ModEx=='lam' and SU_tot==float(4/3):
            return 6152.5, 1.0
        elif J_tot==5/2 and S_tot==1/2 and L_tot==2 and ModEx=='lam' and SU_tot==float(4/3):
            return 6072.3, 1.0
        else:
            return 0,0
    elif baryons=='cascades_anti3':
        if J_tot==1/2 and S_tot==1/2 and L_tot==0 and ModEx=='grd' and SU_tot==float(4/3):
            return 5794.45, 1.0
        elif J_tot==3/2 and S_tot==1/2 and L_tot==1 and ModEx=='lam' and SU_tot==float(4/3):
            return 6100.00, 1.0
        elif J_tot==3/2 and S_tot==1/2 and L_tot==2 and ModEx=='lam' and SU_tot==float(4/3):
            return 6327.00, 1.0
        elif J_tot==5/2 and S_tot==1/2 and L_tot==2 and ModEx=='lam' and SU_tot==float(4/3):
            return 6333.00, 1.0        
        else:
            return 0, 0


def exp_mass_width(baryons, J_tot, S_tot, L_tot, ModEx, SU_tot): # latex names

    if baryons=='omegas':
        if J_tot==1/2 and S_tot==1/2 and L_tot==0 and ModEx=='grd' and SU_tot==float(10/3):
            return '$6045.2\\pm 1.2$' , '$\\approx 0$'
        elif J_tot==1/2 and S_tot==1/2 and L_tot==1 and ModEx=='lam' and SU_tot==float(10/3):
            return '$6315.6\\pm 0.6$', '$<4.2$'
        elif J_tot==3/2 and S_tot==1/2 and L_tot==1 and ModEx=='lam' and SU_tot==float(10/3):
            return '$6339.7\\pm 0.6$' , '$<1.8$'
        elif J_tot==1/2 and S_tot==3/2 and L_tot==1 and ModEx=='lam' and SU_tot==float(10/3):
            return '$6330.3\\pm 0.6$' , '$<4.7$'
        elif J_tot==3/2 and S_tot==3/2 and L_tot==1 and ModEx=='lam' and SU_tot==float(10/3):
            return '$6349.8\\pm 0.6$' , '$<3.2$'
        else:
            return '$\\dagger$', '$\\dagger$'
    elif baryons=='cascades':
        if J_tot==1/2 and S_tot==1/2 and L_tot==0 and ModEx=='grd' and SU_tot==float(10/3):
            return '$5935.02\\pm 0.05$' , '$<0.08$'
        elif J_tot==3/2 and S_tot==3/2 and L_tot==0 and ModEx=='grd' and SU_tot==float(10/3):
            return '$5952.3\\pm 0.6$' , '$0.90\pm 0.18$'
        elif J_tot==5/2 and S_tot==3/2 and L_tot==1 and ModEx=='lam' and SU_tot==float(10/3):
            return '$6227.9\\pm 0.9$' , '$19.9\pm 2.6$'
        else:
            return '$\\dagger$', '$\\dagger$'
    if baryons=='sigmas':
        if J_tot==1/2 and S_tot==1/2 and L_tot==0 and ModEx=='grd' and SU_tot==float(10/3):
            return '$5813.1\\pm 0.3$' , '$5.0\\pm 0.5$'
        elif J_tot==3/2 and S_tot==3/2 and L_tot==0 and ModEx=='grd' and SU_tot==float(10/3):
            return '$5832.5\\pm 0.5$' , '$9.9\\pm 0.9$'
        elif J_tot==1/2 and S_tot==1/2 and L_tot==1 and ModEx=='lam' and SU_tot==float(10/3):
            return '$6096.9\\pm 1.8$' , '$30\\pm 7$'
        else:
            return '$\\dagger$', '$\\dagger$'
    elif baryons=='lambdas':
        if J_tot==1/2 and S_tot==1/2 and L_tot==0 and ModEx=='grd' and SU_tot==float(4/3):
            return '$5619.60\\pm 0.17$' , '$\\approx 0$'
        elif J_tot==1/2 and S_tot==1/2 and L_tot==1 and ModEx=='lam' and SU_tot==float(4/3):
            return '$5912.19\\pm 0.17$' , '$<0.25$'
        elif J_tot==3/2 and S_tot==1/2 and L_tot==1 and ModEx=='lam' and SU_tot==float(4/3):
            return '$5920.09\\pm 0.17$' , '$<0.19$'
        elif J_tot==3/2 and S_tot==1/2 and L_tot==2 and ModEx=='lam' and SU_tot==float(4/3):
            return '$6146.2\\pm 0.4$' , '$2.9\\pm 1.3$'
        elif J_tot==5/2 and S_tot==1/2 and L_tot==2 and ModEx=='lam' and SU_tot==float(4/3):
            return '$6152.5\\pm 0.4$' , '$2.1\\pm 0.9$'
        #elif J_tot==5/2 and S_tot==1/2 and L_tot==2 and ModEx=='lam' and SU_tot==float(4/3):
        #    return '$6072.3\\pm 1.0$' , '$xx$'
        else:
            return '$\\dagger$', '$\\dagger$'
    elif baryons=='cascades_anti3':
        if J_tot==1/2 and S_tot==1/2 and L_tot==0 and ModEx=='grd' and SU_tot==float(4/3):
            return '$5797.5\\pm 0.6$' , '$\\approx 0$'
        elif J_tot==3/2 and S_tot==1/2 and L_tot==1 and ModEx=='lam' and SU_tot==float(4/3):
            return '$6100.3\\pm 0.6$' , '$<1.9$'
        elif J_tot==3/2 and S_tot==1/2 and L_tot==2 and ModEx=='lam' and SU_tot==float(4/3):
            return '$6327.4\\pm 1.6$' , '$1\\pm5$'
        elif J_tot==5/2 and S_tot==1/2 and L_tot==2 and ModEx=='lam' and SU_tot==float(4/3):
            return '$6332.7\\pm 0.3$' , '$<1.60$'
        else:
            return '$\\dagger$', '$\\dagger$'
        

def names_values_paper(baryons, nStates):
    # these are for display reasons only
    if baryons=='omegas':        
        quantum = ['$\\vert ssc,1/2,1/2,0,0,10/3 \\rangle$ &',   '$\\vert ssc,3/2,3/2,0,0,10/3 \\rangle$ &',   '$\\vert ssc,1/2,1/2,1,0,10/3 \\rangle$ &',
                   '$\\vert ssc,1/2,3/2,1,0,10/3 \\rangle$ &',   '$\\vert ssc,3/2,1/2,1,0,10/3 \\rangle$ &',   '$\\vert ssc,1/2,1/2,0,0,10/3 \\rangle$ &',
                   '$\\vert ssc,5/2,1/2,0,1,10/3 \\rangle$ &']
        exp       = np.array([2695.0, 2765.9, 3000.4, 3050.2, 3065.6, 3090.2, 3119.1])
        delta_exp = np.array([   1.7,    2.0,   0.22,   0.13,   0.28,    0.5,    1.0])

    if baryons=='cascades':
        quantum = ['$\\vert suc,1/2,1/2,0,1/2,10/3 \\rangle$ &', '$\\vert suc,3/2,3/2,0,1/2,10/3 \\rangle$ &', '$\\vert suc,1/2,3/2,1,1/2,10/3 \\rangle$ &',
                   '$\\vert suc,3/2,1/2,1,1/2,10/3 \\rangle$ &', '$\\vert suc,3/2,3/2,1,1/2,10/3 \\rangle$ &', '$\\vert suc,3/2,3/2,1,1/2,10/3 \\rangle$ &',
                   '$\\vert suc,3/2,3/2,1,1/2,10/3 \\rangle$ &']                   
        exp =       np.array([2578.0, 2645.9, 0.0, 2923.0, 2938.5, 2964.9, 0.0, 3055.9, 3078.6])
        delta_exp = np.array([   0.9,   0.71, 0.0,   0.35,    0.3,   0.33, 0.0,    0.4,    2.8])

    if baryons=='sigmas':
        quantum = ['$\\vert uuc,1/2,1/2,0,1,10/3 \\rangle$ &',   '$\\vert uuc,3/2,3/2,0,1,10/3 \\rangle$ &',   '$\\vert uuc,1/2,1/2,1,1,10/3 \\rangle$ &']
        exp       = np.array([2453.5, 2518.1, 2800.0])
        delta_exp = np.array([   0.9,    2.8,   20.0])        

    if baryons=='lambdas':
        quantum = ['$\\vert udc,1/2,1/2,0,0,4/3  \\rangle$ &',   '$\\vert udc,1/2,1/2,1,0,10/3 \\rangle$ &',   '$\\vert udc,3/2,1/2,1,0,4/3  \\rangle$ &'
        ]
        exp       = np.array([2286.5, 2592.3, 2625.0
        ])
        delta_exp = np.array([  0.14,   0.28,   0.18
        ])

    if baryons=='cascades_anti3':
        quantum = ['$\\vert ssc,1/2,1/2,0,1/2,4/3 \\rangle$ &',  '$\\vert ssc,1/2,1/2,1,1/2,4/3 \\rangle$ &',  '$\\vert ssc,3/2,1/2,1,1/2,10/3 \\rangle$ &' ]
        exp       = np.array([2469.42, 2793.3, 2818.49, 0.0, 2968.6])
        delta_exp = np.array([   1.77,   0.28,    2.07, 0.0,   3.3])

    # fill arrays with zeros for dimension consistency
    if(nStates-len(exp)>0):
        exp       = np.append(exp, np.zeros(nStates-len(exp)) )
        delta_exp = np.append(delta_exp, np.zeros(nStates-len(delta_exp)))
        for i in range(nStates-len(quantum)): quantum.append('$\\vert xxc,x/2,x/2,x,xx,x/3 \\rangle$ &')
        
    return quantum, exp, delta_exp

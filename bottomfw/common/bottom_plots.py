#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
---------------------------------------------------------------
 code to obtain uncertainties of quarkonium mass spectrum
 authors: A. Ramirez-Morales (andres.ramirez.morales@cern.ch)
          H. Garcia-Tecocoatzi
 ---------------------------------------------------------------
'''
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib.offsetbox import AnchoredText
import pandas as pd
# framework includes
from bottomfw.common import data_utils as du


class BottomPlots:
    """
    Class to produce mass spectra plots
    """    
    def __init__(self, baryons, workpath=".", batch_results=False):
        self.m_baryons = baryons
        self.m_workpath = workpath
        self.m_batch = batch_results
        self.m_plot_title = None
        self.m_data = None
        self.m_diquark = None
        self.m_model_flag = None
                
        
    def load_data(self, model_flag):
        '''loads data either diquark and three-quark'''        
        self.m_model_flag = model_flag
        
        if model_flag=='diquark':
            self.m_plot_title =  r'$\bf{Quark\mathrm{-}diquark\,\,model}$'
            self.m_data = pd.read_csv(self.m_workpath+"/tables/masses_diquark_" + self.m_baryons + "_summary.csv")
            self.m_diquark = "_diquark"
        else:
            self.m_plot_title = r'$\mathrm{\bf{Three\mathrm{-}quark\,\,model}}$'
            self.m_data = pd.read_csv(self.m_workpath+"/tables/masses_" + self.m_baryons + "_summary.csv")
            self.m_diquark = ""
        
    def mass_spectrum_plot(self):
        # timesfont = {'fontname':'Times New Roman'} # check this later
        # plt.xlabel('xlabel', **timesfont)
        # plt.ylabel('ylabel', **timesfont)
        baryon_name = du.baryon_name_plots(self.m_baryons)

        fig, ax = plt.subplots()
        #fig.tight_layout()
        tick_positions = ([])
        tick_labels=['$1/2^+$', '$1/2^-$','$3/2^+$', '$3/2^-$', '$5/2^+$', '$5/2^-$', '$7/2^+$']
        ax.set_xticklabels(tick_labels)
        ax.xaxis.label.set_size(15)
        ax.yaxis.label.set_size(15)
        ax.set_xlabel('$J^P$')
        ax.set_ylabel('Mass (MeV)')

        if self.m_baryons=='omegas':
            ax.set_ylim([6000,7050])
            if self.m_model_flag=="diquark":
                ax.set_ylim([6000,6700])
        elif self.m_baryons=='sigmas':
            ax.set_ylim([5700,7000])
            if self.m_model_flag=="diquark":
                ax.set_ylim([5700,6600])
        elif self.m_baryons=='cascades':
            ax.set_ylim([5800,7050])
            if self.m_model_flag=="diquark":
                ax.set_ylim([5500,7000])
        elif self.m_baryons=='lambdas':
            ax.set_ylim([5500,6800])
            if self.m_model_flag=="diquark":
                ax.set_ylim([5500,6400])
        elif self.m_baryons=='cascades_anti3':
            ax.set_ylim([5700,6900])
            if self.m_model_flag=="diquark":
                ax.set_ylim([5700,6500])

        x1 = np.array([0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00])
        x2 = np.array([0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25])
        j_tot  = np.array([0.50, 1.5, 2.5, 3.5])
        parity = np.array(['p', 'n'])
        off_set = 0.05
        pos = 0

        for i in range(len(j_tot)):
            for j in range(len(parity)):
                pos += 1
                pos2 = pos - 2
                if parity[j] =='p':
                    J_PRED = self.m_data.loc[(self.m_data['J_tot']==j_tot[i]) & (self.m_data['L_tot']!=1)]
                    J_EXP  = self.m_data.loc[(self.m_data['J_tot']==j_tot[i]) & (self.m_data['exp_mass']!=0) & (self.m_data['L_tot']!=1)]
                elif parity[j] =='n':
                    J_PRED = self.m_data.loc[(self.m_data['J_tot']==j_tot[i]) & (self.m_data['L_tot']==1)]
                    J_EXP  = self.m_data.loc[(self.m_data['J_tot']==j_tot[i]) & (self.m_data['exp_mass']!=0) & (self.m_data['L_tot']==1)]

                J_PRED_up = (J_PRED['mass']+J_PRED['error_up']).to_numpy()
                J_PRED_dn = (J_PRED['mass']-J_PRED['error_up']).to_numpy()
        
                J_EXP_up = (J_EXP['exp_mass']+J_EXP['exp_mass_err']).to_numpy()
                J_EXP_dn = (J_EXP['exp_mass']-J_EXP['exp_mass_err']).to_numpy()

                x = [x1[pos-1] + off_set*pos, x2[pos-1] + off_set*pos2]
                vx = np.full(len(J_PRED_dn), (x[1]+x[0])*0.5)
                vx_exp = np.full(len(J_EXP_dn), (x[1]+x[0])*0.5)
                if pos!=8: tick_positions = np.append(tick_positions, (x[1]+x[0])*0.5)
            
                # plot predictions
                plt.plot(vx, J_PRED['mass'], linestyle='None', marker='_',  markersize=28.5, mew=0.95, color='b')                
                for k in range(len(J_PRED_dn)):
                    plt.fill_between(x, J_PRED_dn[k], J_PRED_up[k], color='b', linestyle = 'None', alpha=0.5, edgecolor='None')

                # plot experimental data
                plt.plot(vx_exp, J_EXP['exp_mass'], linestyle='None', marker='_',  markersize=28.5, mew=0.95, color='r')
                for k in range(len(J_EXP_dn)):
                    plt.fill_between(x, J_EXP_dn[k], J_EXP_up[k], color='r', linestyle = 'None', alpha=0.5, edgecolor='None')

        blue_patch = mpatches.Patch(color='b', alpha=0.5, linewidth=0)
        blue_line  = mlines.Line2D([], [], color='b', linewidth=0.95, markersize=0)
        red_patch  = mpatches.Patch(color='r', alpha=0.5, linewidth=0)
        red_line   = mlines.Line2D([], [], color='r', linewidth=0.95, markersize=0)
        plt.legend([(blue_patch, blue_line), (red_patch, red_line)], 
                   ['This work pred.', 'PDG data'], loc=4, fancybox=False, title=self.m_plot_title, title_fontsize='medium')

        ax.set_xticks(tick_positions)
        fig.tight_layout()
        at = AnchoredText(baryon_name, prop=dict(size=17.5), frameon=False, loc='upper left')
        at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
        ax.add_artist(at)
        if not os.path.exists(self.m_workpath+"/plots/"):
            os.mkdir(self.m_workpath+"/plots/")
        plt.savefig(self.m_workpath+'/plots/spectra_'+self.m_baryons+self.m_diquark+'.png')

#!/usr/bin/env python3
# - * - coding: utf-8 - * -
"""
---------------------------------------------------------------
 Authors: A. Ramirez-Morales (andres.ramirez.morales@cern.ch)
          H. Garcia-Tecocoatzi
---------------------------------------------------------------
"""
import numpy as np
# framework includes
from decays import decay_utils as du

def decay_table_channels(baryons, channel_widths, sym_errors, up_errors, dn_errors, symm, f_out):
    """
    Individual channel tables
    """
    channel_widths = np.array(channel_widths)

    no_errors = False # for no bootstrap for decay widths
    if(sym_errors is None or up_errors is None or dn_errors is None):
        no_errors = True
        
    name_decays=[]
    name_decays.append('State')
    for i in range(len(channel_widths[0])):
        name_decays.append(du.latex_decay_label(baryons,i+1))
        
    name_decays.append('Tot $\\Gamma$')                   
    du.print_header_latex(name_decays, f_out)    
    
    for i in range(len(channel_widths)):
        if symm and not no_errors:
            error_up = sym_errors[i]
            error_dn = sym_errors[i]
        elif not no_errors:
            error_up = up_errors[i]
            error_dn = dn_errors[i]
        else:
            error_up = None
            error_dn = None
        
        du.print_row_latex(state_name=str(i+1), state_decays=channel_widths[i],
                           errors_up=error_up, errors_dn=error_dn, f_out=f_out)
        
    du.print_bottom_latex(baryons,f_out)

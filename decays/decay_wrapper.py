#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
---------------------------------------------------------------
 Python wrapper for the DecayWidths C++ shared library
 authors: A. Ramirez-Morales (andres.ramirez.morales@cern.ch)
          H. Garcia-Tecocoatzi
 ---------------------------------------------------------------
'''
import ctypes
from ctypes import cdll


class decay(object):
    def __init__(self, workpath="."):
        self.m_lib = cdll.LoadLibrary(workpath+'/decays/DecayWidths/libbottomdecay.so')
        self.obj = self.m_lib.charm_new()

    def decay_width(self, MA_val, MB_val, MC_val, GA_val, SA_val,
                          LA_val, JA_val, SL_val, AL_val, AR_val,
                          baryon, excMode, prodDecay):
        
        self.m_lib.charm_execute.restype = ctypes.c_double
        MA_val = ctypes.c_double(MA_val)
        MB_val = ctypes.c_double(MB_val)
        MC_val = ctypes.c_double(MC_val)
        GA_val = ctypes.c_double(GA_val)
        SA_val = ctypes.c_double(SA_val)
        LA_val = ctypes.c_double(LA_val)
        JA_val = ctypes.c_double(JA_val)
        SL_val = ctypes.c_double(SL_val)
        AL_val = ctypes.c_double(AL_val)
        AR_val = ctypes.c_double(AR_val)
        baryon = ctypes.c_int(baryon)
        excMode = ctypes.c_int(excMode)
        prodDecay = ctypes.c_int(prodDecay)
        
        decay_value = self.m_lib.charm_execute(self.obj, MA_val, MB_val, MC_val, GA_val,
                                        SA_val,   LA_val, JA_val, SL_val, AL_val, AR_val,
                                        baryon,  excMode, prodDecay)
        return decay_value

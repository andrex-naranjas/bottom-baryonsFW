"""
---------------------------------------------------------------
 Authors: A. Ramirez-Morales (andres.ramirez.morales@cern.ch)
          H. Garcia-Tecocoatzi
---------------------------------------------------------------
"""
import ctypes
from ctypes import CDLL
import os


class decay(object):
    """
    Python wrapper for the DecayWidths C++ shared library
    """
    def __init__(self, workpath="."):
        self.workpath = os.path.dirname(os.path.realpath(__file__))
        
    def decay_width(self, MA_avg_val, MB_avg_val, MC_avg_val, MA_val, MB_val, MC_val,
                    GA_val, SA_val,
                    LA_val, JA_val, SL_val, AL_val, AR_val,
                    baryon, excMode, prodDecay):
        """
        Method to convert the python variables to c++ objects
        """
        MA_avg_val = ctypes.c_double(MA_val)
        MB_avg_val = ctypes.c_double(MB_val)
        MC_avg_val = ctypes.c_double(MC_val)

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
        m_lib = ctypes.CDLL(os.path.join(self.workpath+"/DecayWidths", 'libbottomdecay.so'))
        m_lib.bottom_execute.restype = ctypes.c_double
        m_lib.bottom_execute.argtypes = [ctypes.c_double]
        decay_value = m_lib.bottom_execute(MA_avg_val, MB_avg_val, MC_avg_val, MA_val, MB_val, MC_val,
                                           GA_val, SA_val,
                                           LA_val, JA_val, SL_val, AL_val, AR_val,
                                           baryon, excMode, prodDecay)
        return decay_value


    def electro_width(self, MA_val, SA_val, JA_val, LA_val, SlA_val, LlA_val, LrA_val,
                      MB_val,
                      AL_val, AR_val,
                      mbottom_val, mupdown_val, mstrange_val,
                      baryon, excMode, prodDecay):
        """
        Method to convert the python variables to c++ objects
        """
        MA_val  = ctypes.c_double(MA_val)
        SA_val  = ctypes.c_double(SA_val)
        JA_val  = ctypes.c_double(JA_val)
        LA_val  = ctypes.c_double(LA_val)
        SlA_val = ctypes.c_double(SlA_val)
        LlA_val = ctypes.c_double(LlA_val)
        LrA_val = ctypes.c_double(LrA_val)
        MB_val  = ctypes.c_double(MB_val)
        AL_val = ctypes.c_double(AL_val)
        AR_val = ctypes.c_double(AR_val)
        mbottom_val = ctypes.c_double(mbottom_val)
        mupdown_val = ctypes.c_double(mupdown_val)
        mstrange_val = ctypes.c_double(mstrange_val)
        baryon = ctypes.c_int(baryon)
        excMode = ctypes.c_int(excMode)
        prodDecay = ctypes.c_int(prodDecay)
        m_lib = ctypes.CDLL(os.path.join(self.workpath+"/DecayWidths", 'libbottomdecay.so'))
        m_lib.electro_execute.restype = ctypes.c_double
        m_lib.electro_execute.argtypes = [ctypes.c_double]
        decay_value = m_lib.electro_execute(MA_val, SA_val, JA_val, LA_val, SlA_val, LlA_val, LrA_val,
                                            MB_val,
                                            AL_val, AR_val,
                                            mbottom_val, mupdown_val, mstrange_val,
                                            baryon, excMode, prodDecay)
        return decay_value

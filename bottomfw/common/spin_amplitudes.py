#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
---------------------------------------------------------------
 authors: C. A. Vaquera Araujo (vaquera@fisica.ugto.mx)
          A. Rivero (ailierrivero@gmail.com)
          A. Ramirez Morales (andres@knu.ac.kr)
---------------------------------------------------------------
"""
import sympy.physics.mechanics as mech
#mech.init_vprinting()
from sympy import *
from sympy.physics.quantum import TensorProduct
from sympy.physics.matrices import msigma
import numpy as np

class SpinAmplitudes():
    """
    Class to get electromagnetic spin amplitudes
    """        
    def __init__(self, baryons):
        self.m_baryons = baryons # baryons podria ser omegas, sigmas, lambdas, etc
        self.spin_u = Matrix([[1],[0]]) #spin state uparrow
        self.spin_d = Matrix([[0],[1]]) #spin state downarrow
            
    def su2_generator(self, index):
        """
        Method to calculate SU(2) generator
        index == generator index (1,2,3)
        equivalent j(i)
        """
        return Rational(1,2)*msigma(index)

    def identity_matrix(self, dim):
        """
        Method to obtain the identity matrix
        """
        return eye(dim)

    def tensor_product(self, index, dim=2):
        """
        Method to calculate the direct product 2⊗2⊗2 generators
        J(i)
        """
        TP1 = TensorProduct(TensorProduct(self.su2_generator(index), self.identity_matrix(dim)), self.identity_matrix(dim)) 
        TP2 = TensorProduct(TensorProduct(self.identity_matrix(dim), self.su2_generator(index)), self.identity_matrix(dim))
        TP3 = TensorProduct(TensorProduct(self.identity_matrix(dim), self.identity_matrix(dim)), self.su2_generator(index)) 
        return TP1 + TP2 + TP3

    def ladder_operator(self, sign=1, dim1=1, dim2=2):
        """
        Method to calculate ladder operators
        """
        if sign > 0:
            return self.su2_generator(dim1) + I*self.su2_generator(dim2)
        else:
            return self.su2_generator(dim1) - I*self.su2_generator(dim2)

    def ladder_operator_tensor(self, sign=1, dim1=1, dim2=2):
        """
        Method to calculate ladder operators tensor
        """
        if sign > 0:
            return self.tensor_product(dim1) + I*self.tensor_product(dim2)
        else:
            return self.tensor_product(dim1) - I*self.tensor_product(dim2)
   
    def spint_states(self, state_a, state_b, state_c):
        """
        Method to calculate the spin states
        state_a,b,c should be sympy matrices of the form: Matrix([[0],...,[1]])
        k 
        """
        return TensorProduct(TensorProduct(state_a ,state_b), state_c)

    def symmetric_states(self, m_proj):
        """
        Method to calculate symmetric states
        """
        if m_proj in [3/2, 1/2, -1/2, -3/2]:
            i = 3/2 - m_proj
            st = self.spint_states(self.spin_u, self.spin_u, self.spin_u)
            while i > 0:
                v1 = self.ladder_operator_tensor(sign=-1) * st
                st = v1 / sqrt((transpose(v1)*v1)[0])
                i = i - 1
        return st

    def rho_states(self, m_proj):
        """
        Method to calculate the rho states 
        """
        if m_proj in [1/2, -1/2]:
            i = 1/2 - m_proj
            st = (self.spint_states(self.spin_u, self.spin_d, self.spin_u) - self.spint_states(self.spin_d, self.spin_u, self.spin_u))/sqrt(2)
            while i > 0:
                v1 = self.ladder_operator_tensor(sign=-1) * st
                st = v1/sqrt((transpose(v1)*v1)[0])
                i = i-1
        return st

    def lambda_states(self, m_proj):
        """
        Method to calculate the lambda states
        """
        if m_proj in [1/2, -1/2]:
            i = 1/2 - m_proj
            st = (2*self.spint_states(self.spin_u, self.spin_u, self.spin_d) - self.spint_states(self.spin_u, self.spin_d, self.spin_u) -\
                self.spint_states(self.spin_d, self.spin_u, self.spin_u))/sqrt(6)
            while i>0:
                v1 = self.ladder_operator_tensor(sign=-1) * st
                st = v1/sqrt((transpose(v1)*v1)[0])
                i = i - 1
        return st

    def matrix_elements(self, n, x_i, y_j):
        """
        Method to calculate the matrix element for the spin part
        """
        if n==1:
            sm = TensorProduct(TensorProduct(self.ladder_operator(sign=-1, dim1=1, dim2=2), self.identity_matrix(dim=2)), self.identity_matrix(dim=2))  
        else: 
            if n==2:
                sm = TensorProduct(TensorProduct(self.identity_matrix(dim=2), self.ladder_operator(sign=-1, dim1=1, dim2=2)), self.identity_matrix(dim=2))
            else:
                sm = TensorProduct(TensorProduct(self.identity_matrix(dim=2), self.identity_matrix(dim=2)), self.ladder_operator(sign=-1, dim1=1, dim2=2))
        return np.array(conjugate(transpose(x_i)) * sm * y_j).astype(np.float64).flatten()[0]

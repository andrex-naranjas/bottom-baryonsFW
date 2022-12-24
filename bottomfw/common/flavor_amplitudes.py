#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------
 authors: C. A. Vaquera Araujo (vaquera@fisica.ugto.mx)
          A. Rivero (ailierrivero@gmail.com)
          A. Ramirez Morales (andres@knu.ac.kr)
------------------------------------------------------------
"""
import sympy.physics.mechanics as mech
from sympy import *
from sympy.physics.quantum import TensorProduct
import numpy as np

#Values of the masses computed in our work
m_q = 299
m_s = 465
m_b = 4928

#Values of the magnetic moments of each quark according to Dothan 1982
mu_u = (2/3) * 1/(2 * m_q)
mu_d = (-1/3) * 1/(2 * m_q)
mu_s = (-1/3) * 1/(2 * m_s)
mu_b = (-1/3) * 1/(2 * m_b)


class FlavorAmplitudes():
    """
    Class to get electromagnetic flavor amplitudes
    """        
    def __init__(self, baryons):        
        # baryons available: omegas, sigmas, lambdas, cascades_prime, cascades
        self.m_baryons = baryons  
        # up quark state
        self.u = Matrix([[1],[0],[0],[0]]) 
        # down quark state
        self.d = Matrix([[0],[1],[0],[0]]) 
        # strange quark state
        self.s = Matrix([[0],[0],[1],[0]]) 
        # bottom quark state
        self.b = Matrix([[0],[0],[0],[1]]) 
              
    def magnetic_operators(self, index): 
        """
        Method to calculate the Magnetic Operators (mu(i))
        index == operator index (1,2,3)
        """
        return Matrix([[mu_u, 0, 0, 0], [0, mu_d, 0, 0], [0, 0, mu_s, 0], [0, 0, 0, mu_b]])

    def identity_matrix(self, dim):
        """
        Method to obtain the identity matrix
        """
        return eye(dim)

    def tensor_product_operator_1(self, index=1, dim=4): 
        """
        Method to define the direct product 2⊗2⊗2 magnetic operator for mu_1 (Mu(1))
        """
        return TensorProduct(TensorProduct(self.magnetic_operators(index), self.identity_matrix(dim)), self.identity_matrix(dim))

    def tensor_product_operator_2(self, index=2, dim=4):
        """
        Method to define the direct product 2⊗2⊗2 magnetic operator for mu_2 (Mu(2))
        """
        return TensorProduct(TensorProduct(self.identity_matrix(dim), self.magnetic_operators(index)), self.identity_matrix(dim))

    def tensor_product_operator_3(self, index=3, dim=4):
        """
        Method to define the direct product 2⊗2⊗2 magnetic operator for mu_3 (Mu(3))
        """
        return TensorProduct(TensorProduct(self.identity_matrix(dim), self.identity_matrix(dim)), self.magnetic_operators(index))        
    
    def flavor_state(self, state_a, state_b, state_c): 
        """
        Method to calculate the flavor states (Fs)
        state_a,b,c should be sympy matrices of the form: Matrix([[0],...,[1]])
        """
        return TensorProduct(TensorProduct(state_a ,state_b), state_c)
    
    
    #Flavor States of the anti-triplet
    def Lambda_b_0(self):
        """
        Method to define the Lambda_b_0 flavor state
        """
        return (self.flavor_state(self.u, self.d, self.b) - self.flavor_state(self.d, self.u, self.b)) / sqrt(2)

    def Xi_b_0(self):
        """
        Method to define the Xi_b_0 flavor state
        """
        return (self.flavor_state(self.u, self.s, self.b) - self.flavor_state(self.s, self.u, self.b)) / sqrt(2)

    def Xi_b_m(self):
        """
        Method to define the Xi_b_m flavor state
        """
        return (self.flavor_state(self.d, self.s, self.b) - self.flavor_state(self.s, self.d, self.b)) / sqrt(2)


    #Flavor States of the sextet
    def Sigma_b_p(self):
        """
        Method to define the Sigma_b_p flavor state
        """
        return (self.flavor_state(self.u, self.u, self.b))

    def Sigma_b_0(self):
        """
        Method to define the Sigma_b_0 flavor state
        """
        return (self.flavor_state(self.u, self.d, self.b) + self.flavor_state(self.d, self.u, self.b)) / sqrt(2)
  
    def Sigma_b_m(self):
        """
        Method to define the Sigma_b_p flavor state
        """
        return (self.flavor_state(self.d, self.d, self.b))

    def Xi_b_pr_0(self):
        """
        Method to define the Xi_b_pr_0 flavor state
        """
        return (self.flavor_state(self.u, self.s, self.b) + self.flavor_state(self.s, self.u, self.b)) / sqrt(2)

    
    def Xi_b_pr_m(self):
        """
        Method to define the Xi_b_pr_m flavor state
        """
        return (self.flavor_state(self.d,self.s,self.b) + self.flavor_state(self.s, self.d, self.b)) / sqrt(2)

    def Omega_b_m(self):
        """
        Method to define the Omega_b_m flavor state
        """
        return (self.flavor_state(self.s, self.s, self.b))
    
    def flavor_matrix_elements(self, n, x_i, y_j):
        """
        Method to calculate the matrix element for the flavor part
        """
        if n==1:
            Mu_i=TensorProduct(TensorProduct(self.magnetic_operators(index=1), self.identity_matrix(dim=4)), self.identity_matrix(dim=4))  
        else: 
            if n==2:
                Mu_i=TensorProduct(TensorProduct(self.identity_matrix(dim=4), self.magnetic_operators(index=2)), self.identity_matrix(dim=4))
            else:
                Mu_i=TensorProduct(TensorProduct(self.identity_matrix(dim=4), self.identity_matrix(dim=4)), self.magnetic_operators(index=3))
        return np.array(conjugate(transpose(x_i)) * Mu_i * y_j).astype(np.float64).flatten()[0]
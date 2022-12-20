#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
---------------------------------------------------------------
 authors: C. A. Vaquera Araujo (vaquera@fisica.ugto.mx)
          A. Rivero (ailierrivero@gmail.com)
          A. Ramirez Morales (andres@knu.ac.kr)
 ---------------------------------------------------------------
'''

import sympy.physics.mechanics as mech
#mech.init_vprinting()
from sympy import *
from sympy.physics.quantum import TensorProduct

#mu_u, mu_d, mu_s, mu_b=symbols('mu_u, mu_d, mu_s, mu_b',positive=True) # Definitions of the magnetic moments of each quark
mu_u = 1
mu_d = 2
mu_s = 3
mu_b = 4

class FlavorAmplitudes():
    """
    Class to get electromagnetic flavor amplitudes
    """    
    
    def __init__(self, baryons):

        self.m_baryons = baryons  
        # baryons podria ser omegas, sigmas, lambdas, etc
        self.u = Matrix([[1],[0],[0],[0]]) 
        # up quark state
        self.d = Matrix([[0],[1],[0],[0]]) 
        # down quark state
        self.s = Matrix([[0],[0],[1],[0]]) 
        # strange quark state
        self.b = Matrix([[0],[0],[0],[1]]) 
        # bottom quark state
       
 
    def magnetic_operators(self, index): 
        #mu(i)
        """
        Method to calculate the Magnetic Operators
        index == operator index (1,2,3)
        """
        return Matrix([[mu_u, 0, 0, 0],[0, mu_d, 0, 0],[0, 0, mu_s, 0],[0, 0, 0, mu_b]])


    def identity_matrix(self, dim):
        """
        Method to obtain the identity matrix
        """
        return eye(dim)


    def tensor_product_operator_1(self, index=1, dim=4): 
        # Mu(1):
        """
        Method to define the direct product 2⊗2⊗2 magnetic operator for mu_1 
        """
        return TensorProduct(TensorProduct(self.magnetic_operators(index), self.identity_matrix(dim)), self.identity_matrix(dim))

    def tensor_product_operator_2(self, index=2, dim=4): 
        # Mu(2):
        """
        Method to define the direct product 2⊗2⊗2 magnetic operator for mu_2 
        """
        return TensorProduct(TensorProduct(self.identity_matrix(dim), self.magnetic_operators(index)), self.identity_matrix(dim))

    def tensor_product_operator_3(self, index=3, dim=4): 
        # Mu(3):
        """
        Method to define the direct product 2⊗2⊗2 magnetic operator for mu_3 
        """
        return TensorProduct(TensorProduct(self.identity_matrix(dim), self.identity_matrix(dim)), self.magnetic_operators(index))        
    

    def flavor_state(self, state_a, state_b, state_c): 
        # Fs
        """
        Method to calculate the flavor states
        state_a,b,c should be sympy matrices of the form: Matrix([[0],...,[1]])
        """
        return TensorProduct(TensorProduct(state_a ,state_b), state_c)

    
    def Sigma_0(self):
        """
        Method to define the Sigma_0 flavor state
        """
        return (self.flavor_state(self.u,self.d,self.b)+self.flavor_state(self.d,self.u,self.b))/sqrt(2)

    def Lambda_0(self):
        """
        Method to define the Lambda_0 flavor state
        """
        return (self.flavor_state(self.u,self.d,self.b)-self.flavor_state(self.d,self.u,self.b))/sqrt(2)

    def Xi_prime_0(self):
        """
        Method to define the Xi_prime_0 flavor state
        """
        return (self.flavor_state(self.u,self.s,self.b)+self.flavor_state(self.s,self.u,self.b))/sqrt(2)

    def Xi_0(self):
        """
        Method to define the Xi_0 flavor state
        """
        return (self.flavor_state(self.u,self.s,self.b)-self.flavor_state(self.s,self.u,self.b))/sqrt(2)

    def Xi_prime_m(self):
        """
        Method to define the Xi_prime_m flavor state
        """
        return (self.flavor_state(self.d,self.s,self.b)+self.flavor_state(self.s,self.d,self.b))/sqrt(2)

    def Xi_m(self):
        """
        Method to define the Xi_m flavor state
        """
        return (self.flavor_state(self.d,self.s,self.b)-self.flavor_state(self.s,self.d,self.b))/sqrt(2)

    
    def transition_lambda0_sigma0_1(self):
        """
        Method to calculate the matrix element for the flavor part
        """
        return conjugate(transpose(self.Lambda_0()))*self.tensor_product_operator_1()*self.Sigma_0()

    def transition_lambda0_sigma0_2(self):
        """
        Method to calculate the matrix element for the flavor part
        """
        return conjugate(transpose(self.Lambda_0()))*self.tensor_product_operator_2()*self.Sigma_0()    

    def transition_lambda0_sigma0_3(self):
        """
        Method to calculate the matrix element for the flavor part
        """
        return conjugate(transpose(self.Lambda_0()))*self.tensor_product_operator_3()*self.Sigma_0()    


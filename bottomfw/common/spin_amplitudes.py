#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
---------------------------------------------------------------
 authors: A. Rivero ()
          A. Ramirez Morales (andres@knu.ac.kr)
 ---------------------------------------------------------------
'''

import sympy.physics.mechanics as mech
#mech.init_vprinting()
from sympy import *
from sympy.physics.quantum import TensorProduct
from sympy.physics.matrices import msigma

u = Matrix([[1],[0]])
d = Matrix([[0],[1]])

class SpinAmplitudes():
    """
    Class to get electromagnetic spin amplitudes
    """    
    
    def __init__(self, baryons):

        self.m_baryons = baryons
        # baryons podria ser omegas, sigmas, lambdas, etc

    
    def test(self):
        dummy = self.m_baryons


    def su2_generator(self, index):
        """
        Method to calculate SU(2) generator
        index == generator index (1,2,3)
        """
        return Rational(1,2)*msigma(index)


    def identity_matrix(self, dim):
        """
        Method to obtain the identity matrix
        """
        return eye(dim)


    def tensor_product(self, index, dim=2): # J(i):
        """
        Method to calculate the tensor product
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
        

    def spint_states(self, state_a=u, state_b=d, state_c=d): # k
        """
        Method to calculate the spin states
        state_a,b,c should be sympy matrices of the form: Matrix([[0],...,[1]])
        """
        return TensorProduct(TensorProduct(state_a ,state_b), state_c)

    
    def symmetric_states(self, m_proj):
        """
        Method to calculate symmetric states
        """
        if m_proj in [3/2,1/2,-1/2,-3/2]:
            i = 3/2-m_proj
            st = spint_states(u, u, u)
            while i > 0:
                v1 = self.ladder_operator_tensor(sign=-1) * st
                st = v1 / sqrt((transpose(v1)*v1)[0])
                i = i - 1
        return st


    def lambda_states(self, a):
        """
        Method to calculate the lambda states
        """
        if a in [1/2,-1/2]:
            i=1/2-a
            st=(2*k(u,u,d)-k(u,d,u)-k(d,u,u))/sqrt(6)
            while i>0:
                v1=Jm*st
                st=v1/sqrt((transpose(v1)*v1)[0])
                i=i-1
        return st

    def rho_states(self, a):
        """
        Method to calculate the rho states 
        """
        if a in [1/2,-1/2]:
            i=1/2-a
            st=(k(u,d,u)-k(d,u,u))/sqrt(2)
            while i>0:
                v1=Jm*st
                st=v1/sqrt((transpose(v1)*v1)[0])
                i=i-1
        return st


    def matrix_elements(self, n, x, i, y, j):
        """
        Method to calculate the matrix element for the spin part
        """
        if n==1:
            sm = TensorProduct(TensorProduct(self.su2_generator(index), self.identity_matrix(dim)), self.identity_matrix(dim))  
        else: 
            if n==2:
                sm = TensorProduct(TensorProduct(self.identity_matrix(dim), self.su2_generator(index)), self.identity_matrix(dim))
            else:
                sm = TensorProduct(TensorProduct(self.identity_matrix(dim), self.identity_matrix(dim)), self.su2_generator(index))
        return conjugate(transpose(x(i)))*sm*y(j)


    def all_matrix_elements(self):
        """
        Method to calcualte all the matrix elements
        """
        for n in [1,2,3]:
            for x in [s,r,l]:
                for y in  [s,r,l]:
                    if x==s:
                        for i in [3/2,1/2,-1/2,-3/2]:
                            if y==s:
                                for j in [3/2,1/2,-1/2,-3/2]:
                                    print([m(n, x, i, y, j)[0],n,x.__name__,Rational(i),y.__name__,Rational(j)])
                            else:
                                for j in [1/2,-1/2]:
                                    print([m(n, x, i, y, j)[0],n,x.__name__,Rational(i),y.__name__,Rational(j)])
                    else:
                        for i in [1/2,-1/2]:
                            if y==s:
                                for j in [3/2,1/2,-1/2,-3/2]:
                                    print([m(n, x, i, y, j)[0],n,x.__name__,Rational(i),y.__name__,Rational(j)])
                            else:
                                for j in [1/2,-1/2]:
                                    print([m(n, x, i, y, j)[0],n,x.__name__,Rational(i),y.__name__,Rational(j)])

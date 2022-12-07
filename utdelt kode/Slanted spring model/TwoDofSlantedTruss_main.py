# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 21:55:06 2017

@author: bjohau
"""
import math
import numpy as np
import matplotlib.pyplot as plt
import TwoDofSlantedTruss as struct

Ks = 10.0 # Stiffness of slanted element
Kv = 10.0 # Stiffness of vertical element
#Kv = 0.10
Kv =  0.10
LH  = 1.0 # Horizontal length of slanted truss
LV  = -0.2 # Vertical length of slanted truss
F0 = 0.1

problem = struct.SlantedTrussStructure(LH, LV, Ks, Kv, F0)

uVec   = np.zeros(2)
res_Vec   = np.zeros(2)
Lambda = 0.0
arcLength = 0.02

nSteps = 50
maxIter = 30
d_q_prev = np.zeros(2)

lambdaList =[]
u1List = []
u2List = []

for iStep in range(nSteps):

    for iIter in range(maxIter):

        res_Vec = problem.get_residual(uVec,Lambda)
        q_Vec   = problem.get_vecq(Lambda)

        K_mat   = problem.get_Kmat(uVec)

        d_r = np.linalg.solve(K_mat,res_Vec)
        d_q = np.linalg.solve(K_mat,q_Vec)

        dLambda = 0
        if ( iIter == 0 ): # Predictor step
            f = math.sqrt(1.0 + d_q.dot(d_q))
            dLambda = arcLength / f

            dotV = d_q.dot(d_q_prev)
            if ( dotV < 0.0):
                dLambda = -dLambda

            d_q_prev = dLambda * d_q

        else: # Corrector step
            f = (1.0 + d_q.dot(d_q))
            dLambda = -(d_r.dot(d_q)) / f

        delta_uVec = d_r + dLambda * d_q
        uVec = uVec + delta_uVec
        Lambda += dLambda

        print("i = {:2d}{:3d} dotV={:8.1e} dLambda={:10.1e} Lambda={:12.3e} u =[{:12.3e} {:12.3e}],  r = [{:12.3e} {:12.3e}]".format(iStep, iIter, dotV, dLambda, Lambda, uVec[0], uVec[1], res_Vec[0], res_Vec[1]))
        if ( False ):
            lambdaList.append(Lambda)
            u1List.append(uVec[0])
            u2List.append(uVec[1])

        res_Vec = problem.get_residual(uVec,Lambda)
        if ( res_Vec.dot(res_Vec) < 1.0e-15 ):
            break

    lambdaList.append(Lambda)
    u1List.append(uVec[0])
    u2List.append(uVec[1])
    print(" ")

plt.grid()
plt.plot(u1List,lambdaList, label='')
plt.plot(u2List,lambdaList)
plt.show()







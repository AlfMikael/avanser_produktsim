# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 16:43:51 2018

@author: bjohau
"""

import numpy as np
import math

def rot_matrix(theta):
    """
    Return the 2x2 rotation matrix representing a rotation theta
    :param theta:  rotation angle in radians
    :return: Rotation matrix (or tensor)
    """
    s = math.sin(theta)
    c = math.cos(theta)
    R = np.array([[c, -s],
                  [s,  c]])
    return R


def theta_of_deformed(ex, ey, disp_global):
    """

    :param ex: element x coordinate [x1, x2] in undeformed position
    :param ey: element y coordinate [y1, y2] in undeformed position
    :param disp_global:  displacement vector [u1, v1, r1, u2, v2, r2] in global directions
    :return: disp_local_def: displacement vector [u1, v1, r1, u2, v2, r2] in local directions
    """
    eVec12 = np.array([ex[1] - ex[0], ey[1] - ey[0]])  # ex0

    L0 = np.linalg.norm(eVec12)

    ex0 = eVec12 / L0

    # Finner de deformerte rotasjonene
    x1, x2, y1, y2 = ex[0], ex[1], ey[0], ey[1]
    u1, v1, u2, v2 = disp_global[0], disp_global[1], disp_global[3], \
                     disp_global[4]
    r1, r2 = disp_global[2], disp_global[5]

    E_xn = [(x2 + u2) - (x1 + u1),
            (y2 + v2) - (y1 + v1)]
    E_xn = np.array(E_xn)

    Ld = np.linalg.norm(E_xn)  # Nye lengde
    e_xn = E_xn / Ld

    e_yn = [-e_xn[1],
            e_xn[0]]
    e_yn = np.array(e_yn)
    # e_yn = e_yn / np.linalg.norm(e_yn)

    theta = math.atan(e_xn[1] / e_xn[0])
    return theta


def beam2local_def_disp(ex, ey, disp_global):
    """

    :param ex: element x coordinate [x1, x2] in undeformed position
    :param ey: element y coordinate [y1, y2] in undeformed position
    :param disp_global:  displacement vector [u1, v1, r1, u2, v2, r2] in global directions
    :return: disp_local_def: displacement vector [u1, v1, r1, u2, v2, r2] in local directions
    """
    eVec12 = np.array([ex[1] - ex[0], ey[1] - ey[0]])  # ex0

    # L0 = math.sqrt(eVec12 @ eVec12)
    L0 = np.linalg.norm(eVec12)

    # Retningsvektor
    ex0 = eVec12 / L0

    # Henter ut posisjonene
    x1, x2, y1, y2 = ex[0], ex[1], ey[0], ey[1]

    # Henter ut nodeforskyvningene
    u1, v1, u2, v2 = disp_global[0], disp_global[1], disp_global[3], \
                     disp_global[4]

    # Henter ut rotasjonene
    r1, r2 = disp_global[2], disp_global[5]

    # Vektorer for deformert element
    E_xn = [(x2 + u2) - (x1 + u1),
            (y2 + v2) - (y1 + v1)]
    E_xn = np.array(E_xn)

    # Lengde på nytt element
    Ld = np.linalg.norm(E_xn)

    # Retningsvektorer for nytt element
    e_xn = E_xn / Ld

    e_yn = [-e_xn[1],
            e_xn[0]]
    e_yn = np.array(e_yn)  # Konverterer til numpy array

    # Rotasjonsmatriser

    R1 = rot_matrix(r1)
    R2 = rot_matrix(r2)

    t1 = R1 @ ex0
    t2 = R2 @ ex0

    theta1_def = math.asin(e_yn @ t1)
    theta2_def = math.asin(e_yn @ t2)


    def_disp_local = np.array([-0.5 * (Ld - L0),
                               0.0,
                               theta1_def,
                               0.5 * (Ld - L0),
                               0.0,
                               theta2_def])
    #print("Deformert elementvektor", def_disp_local)
    return def_disp_local



def L_deformed(ex, ey, disp_global):
    """

    :param ex: element x coordinate [x1, x2] in undeformed position
    :param ey: element y coordinate [y1, y2] in undeformed position
    :param disp_global:  displacement vector [u1, v1, r1, u2, v2, r2] in global directions
    :return: disp_local_def: displacement vector [u1, v1, r1, u2, v2, r2] in local directions
    """
    eVec12 = np.array([ex[1] - ex[0], ey[1] - ey[0]])  # ex0

    # L0 = math.sqrt(eVec12 @ eVec12)
    L0 = np.linalg.norm(eVec12)

    ex0 = eVec12 / L0

    # Finner de deformerte rotasjonene
    x1, x2, y1, y2 = ex[0], ex[1], ey[0], ey[1]
    u1, v1, u2, v2 = disp_global[0], disp_global[1], disp_global[3], \
                     disp_global[4]
    r1, r2 = disp_global[2], disp_global[5]

    E_xn = [(x2 + u2) - (x1 + u1),
            (y2 + v2) - (y1 + v1)]
    E_xn = np.array(E_xn)

    Ld = np.linalg.norm(E_xn)  # Nye lengde

    return Ld


def beam2corot_Ke_and_Fe(ex,ey,ep, disp_global):
    """
    Compute the stiffness matrix and internal forces for a two dimensional beam element
    relative to deformed configuration.
    
    :param list ex: element x coordinates [x1, x2]
    :param list ey: element y coordinates [y1, y2]
    :param list ep: element properties [E, A, I], E - Young's modulus, A - Cross section area, I - Moment of inertia
    :param list disp_global displacement vector for the element [tx1,ty1,rz1,tx2,ty2,rz2]


    :return mat Ke_global: element stiffness matrix [6 x 6]
    :return mat fe_int_global: element internal forces [6 x 1]
    """
    # Undeformed length and unit vector along element
    eVec12 = np.array([ex[1] - ex[0], ey[1] - ey[0]])
    L0 = math.sqrt(eVec12 @ eVec12)
    eVec12 /= L0

    # Deformed position and unit vector along element
    ex_def = ex + [disp_global[0], disp_global[3]]
    ey_def = ey + [disp_global[1], disp_global[4]]

    # TODO: Quite a bit here
    v_local = beam2local_def_disp(ex, ey, disp_global)

    # Length of element before deformation
    L = L0  # I denne konteksten

    # Material and local forces
    Ke_mat = beam2local_stiff(L0, ep)
    fe_int_local = Ke_mat @ v_local

    # Finding and creating the local geometric stiffness matrix
    Ld = L_deformed(ex, ey, disp_global)

    FMat = np.array([[-fe_int_local[1], fe_int_local[0], 0.0, -fe_int_local[4], fe_int_local[3], 0.0]])
    GMat = np.array([[0.0, -1 / Ld, 0.0, 0.0, 1 / Ld, 0.0]])
    Ke_geo = (FMat.T @ GMat + GMat.T @ FMat) * 0.5

    # Total local stiffness matrix
    Ke_local = Ke_geo + Ke_mat

    # Transform to global stiffness matrix
    Te = beam2corot_Te(ex_def, ey_def)
    fe_int_global = Te @ fe_int_local  # A 6x6 matrix
    Ke_global = Te.T @ Ke_local @ Te  # A 1x6 matrix

    return Ke_global, fe_int_global

    
def beam2corot_Te(ex,ey):
    """
    Compute the transformation matrix for an element
    
    :param list ex: element x coordinates [x1, x2]
    :param list ey: element y coordinates [y1, y2]
    :param list ep: element properties [E, A, I], E - Young's modulus, A - Cross section area, I - Moment of inertia   
    :param list eq: distributed loads, local directions [qx, qy]
    :return mat Te: element transformation from global to local
    """

    n = np.array([ex[1]-ex[0],ey[1]-ey[0]])
    L = np.linalg.norm(n)
    n = n / L  
    
    Te=np.array([
        [ n[0], n[1],  0.,    0.,    0.,   0.],
        [-n[1], n[0],  0.,    0.,    0.,   0.],
        [0.,    0.,    1.,    0.,    0.,   0.],
        [0.,    0.,    0.,   n[0],  n[1],  0.],
        [0.,    0.,    0.,  -n[1],  n[0],  0.],
        [0.,    0.,    0.,    0.,    0.,   1.]
    ])
    

    return Te
    
    
def beam2local_stiff(L:float, ep):  # endret til float
    """
    Compute the stiffness matrix for a two dimensional beam element.
    
    :param list L : element length
    :param list ep: element properties [E, A, I], E - Young's modulus, A - Cross section area, I - Moment of inertia   
    :return mat Kle: element stiffness matrix [6 x 6]
    """

    # Kommentar: L er en float, ikke en list som referert til over

    E=ep[0]
    A=ep[1]
    I=ep[2]
        
    Kle = np.array([
        [E*A/L,              0.,           0.,    -E*A/L,            0.,           0. ],
        [   0.,    12*E*I/L**3.,  6*E*I/L**2.,        0., -12*E*I/L**3.,  6*E*I/L**2. ],
        [   0.,     6*E*I/L**2.,      4*E*I/L,        0.,  -6*E*I/L**2.,     2*E*I/L  ],
        [-E*A/L,             0.,           0.,     E*A/L,            0.,           0. ],
        [   0.,   -12*E*I/L**3., -6*E*I/L**2.,        0.,  12*E*I/L**3., -6*E*I/L**2. ],
        [   0.,     6*E*I/L**2.,      2*E*I/L,        0.,  -6*E*I/L**2.,      4*E*I/L ]
    ])
     
    return Kle


def beam2e(ex, ey, ep, eq=None):
    """
    Compute the linear stiffness matrix for a two dimensional beam element.
    Largely from CALFEM core module

    :param list ex: element x coordinates [x1, x2]
    :param list ey: element y coordinates [y1, y2]
    :param list ep: element properties [E, A, I], E - Young's modulus, A - Cross section area, I - Moment of inertia
    :param list eq: distributed loads, local directions [qx, qy]
    :return mat Ke: element stiffness matrix [6 x 6]
    :return mat fe: element consistent force for distributed load [6 x 1] (if eq!=None)
    """

    n = np.array([ex[1] - ex[0], ey[1] - ey[0]])
    L = np.linalg.norm(n)
    n = n / L

    qx = 0.
    qy = 0.
    if not eq is None:
        qx = eq[0]
        qy = eq[1]

    Kle = beam2local_stiff(L,ep)

    fle = L * np.mat([qx / 2,
                      qy / 2,
                      qy * L / 12,
                      qx / 2,
                      qy / 2,
                      -qy * L / 12]).T

    Te = beam2corot_Te(ex,ey)

    Ke = Te.T @ Kle @ Te
    fe = Te.T @ fle

    if eq is None:
        return Ke
    else:
        return Ke, fe

import math
import numpy as np

class SlantedTrussStructure:

    def __init__(self, LH, LV, Ks, Kv, F0):
        self.LH = LH # Horizontal length
        self.LV = LV  # Vertical height
        self.Ks = Ks  # Stiffness of slanted truss
        self.Kv = Kv  # Stiffness of vertical truss
        self.F0 = F0  # Nominal load, to be scaled with lambda

        self.L0 = math.sqrt(self.LH**2 + self.LV**2) #Inital slanted length


    def get_vecQ(self,Lambda):
        # Return the forces that are applied to the nodes from the external forces

        vecQ = np.array([0, Lambda * self.F0])
        return vecQ

    def get_incremental_load(self):
        # Return the nominal forces that are applied to the nodes from the external forces
        vecq = np.array([0, self.F0])
        return vecq

    def get_internal_forces(self, u):
        # Return the forces that are applied to the nodes from the elements

        # Force from slanted part
        LV_d = self.LV + u[0]
        L_d = math.sqrt(self.LH**2 + LV_d**2)
        N1_slanted = (L_d - self.L0) * self.Ks
        sinAlpha = LV_d / L_d
        N1_vertical = N1_slanted * sinAlpha

        # Force from vertical part
        N2_vertical = ( u[1] - u[0] ) * self.Kv

        vecInternalForce = np.array([-N1_vertical + N2_vertical,
                                     -N2_vertical])
        return vecInternalForce

    def get_residual(self, u, Lambda):
        vecQ = self.get_vecQ(Lambda)
        vecInternal = self.get_vecInternalForces(u,Lambda)
        return (vecQ + vecInternal)

    def get_K_sys(self, u):
        Kmat = np.zeros((2,2))

        # Vertical element

        Kmat[0,0] += self.Kv
        Kmat[0,1] += (-self.Kv)
        Kmat[1,0] += (-self.Kv)
        Kmat[1,1] += self.Kv

        # Slanted element, material stiffness

        LV_d = self.LV + u[0]
        L_d = math.sqrt(self.LH**2 + LV_d**2)
        N1_slanted = (L_d - self.L0) * self.Ks
        sinAlpha = LV_d / L_d
        cosAlpha = self.LH / L_d

        Kmat_vert = sinAlpha * self.Ks * sinAlpha
        Kmat[0,0] += Kmat_vert

        # Geometric stiffness
        Kgeo_vert = cosAlpha * N1_slanted * cosAlpha / L_d
        Kmat[0,0] += Kgeo_vert

        return Kmat





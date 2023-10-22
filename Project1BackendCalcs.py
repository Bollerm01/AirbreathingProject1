





'''
AEEM4063 Project 1 Backend Calculations

Computes the F, mdot, fuel flow f, and efficiencies for a real or ideal cycle HBPR Turbofan

Authors: Matt Boller, Pierce Elliott

University of Cincinnati FS23

'''


# Begin imports
import numpy as np

# Begin computation class

class TFComputation:

    def __init__(self):
        self.M = 0.85 # Mach number
        self.h = 11e3 # Altitude, m
        self.W_TO = 370000 # Max takeoff weight, lbf
        self.W = 0.8*self.W_TO # Cruise weight, lbf
        self.BPR = 10 # bypass ratio
        self.y_c = 1.4 # Gamma cold section
        self.y_h = 1.333 # Gamma hot section (post combustion)
        self.Q = 43100 # Enthalpy of formation for fuel, kJ/kg
        self.pi_f = 1.5 # Fan pressure ratio
        self.pi_c = 36 # Compressor pressure ratio
        self.pi_b = 0.96 # Pressure loss across the combustor
        self.n_i = 0.98 # Inlet efficiency
        self.n_inf_f = 0.89 # Fan polytropic efficiency
        self.n_inf_c = 0.90 # Compressor polytropic efficiency
        self.n_b = 0.99 # Combustor isentropic efficiency
        self.n_inf_t = 0.90 # Turbine polytropic efficiency
        self.n_m = 0.99 # Mechanical efficiency
        self.n_j = 0.99 # Nozzle efficiency
        self.p0 = 0.227e5 # Ambient pressure - ISA table, Pa
        self.S_W = 285 # Wing area, m**2

        # Lift and drag calculations
        L = self.W
        q = self.y_c/2 * self.p0*self.M**2
        CL = L/(q*self.S_W)

        CD = 0.056*CL**2 - 0.004*CL + 0.014
        D = CD*q*self.S_W

        self.T_r = D # Thrust required
                
        return True
        
    def fullCycleCalc(self):


        return [(F/mdot), TSFC, f, thermoEff, propEff, overEff]
     
    def fanCalc(self):


        return True

    def compressCalc(self):
        
        return True
    
    def combustorCalc(self):

        return True
    
    def turbCalc(self):

        return True
    
    def effCalc(self):

        return True

    


        





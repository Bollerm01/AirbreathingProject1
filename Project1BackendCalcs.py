





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

    


        





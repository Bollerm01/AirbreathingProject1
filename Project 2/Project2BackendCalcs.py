import numpy as np


'''
AEEM4063 Project 2 Backend Calculations

Computes compressor and turbine states, angles, and velocities for the multi-stage compressor and turbine design

Authors: Matt Boller, Pierce Elliott

University of Cincinnati FS23

'''

# Begin the code
class TurboMachineryComputation:

    def __init__(self):
        self.BPR = 6.0 # Bypass ratio
        self.mdot = 280.0 # Total mass flow, kg/s
        # Sea level static conditions
        self.pa = 1.01325 # Ambient pressure, bar
        self.Ta = 288.15 # Ambient temperature, K
        self.a = 340.3 # Speed of sound, m/s
        self.rho_a = 1.2250 # Ambient density, kg/m^3
        self.R = 287 # Gas constant, J/(kg*K)
        # Compressor parameters
        self.cpr = 15 # Compressor pressure ratio
        self.n_inf_c = 0.9 # Lower bound of polytropic efficiency
        self.pi_f = 1.4 # Fan pressure ratio
        self.n_inf_f = 0.92 # Fan polytropic efficiency
        self.haller = 0.65 # Lower bound on de Haller criteria
        self.TMach_r = 1.2 # Upper bound on tip Mach number
        self.y_c = 1.4 # Gamma in the cold section
        self.cp_c = 1.005 # Specific heat cold section, kJ/(kg*K)
        self.rr_rt_o = 0.5 # Root to tip ratio at entry to the compressor
        # Turbine parameters
        self.T_04 = 2275 # Turbine inlet temperature
        self.n_inf_t = 0.92 # Turbine polytropic efficiency
        self.phi = 0.78 # Lower bound on flow coefficient
        self.psi = 3.3 # Upper bound on loading coefficient
        self.y_h = 1.333 # Gamma in the hot section
        self.cp_c = 1.148 # Specific heat hot section, kJ/(kg*K)
        # Basic design parameters - constant tip radius assumption
        self.U_t = 350 # Max Tip speed, m/s
        self.C_a = 150 # Inlet axial velocity, m/s

        # Calcs for inlet conditions to the compressor
        # First, calculate fan outlet (assume intake has no loss)
        self.P_02 = self.pi_f*self.pa
        self.T_02 = self.Ta + self.Ta*((self.pi_f)**((self.y_c-1)/(self.y_c*self.n_inf_f)) - 1)

        T2 = self.T_02 - self.C_a**2/(2*self.cp_c*1000)
        self.p2 = self.P_02*(T2/self.T_02)**(self.y_c/(self.y_c-1))  
        self.rho2 = self.p2*1e5/(self.R*T2)

        self.rt = np.sqrt((self.mdot/(self.BPR+1))/(np.pi*self.rho2*self.C_a*(1-self.rr_rt_o**2)))
        self.N = self.U_t/(2*np.pi*self.rt)


        return 
    
    def update(self):
        return self.N, self.rt
    
    def fullturbo():

        return
    
    def fullcompressor():

        return
    
    def fullturbine():

        return
    
    def compressorstage():

        return
    
    def turbinestage():

        return
    
    def velocitytriangle():

        return

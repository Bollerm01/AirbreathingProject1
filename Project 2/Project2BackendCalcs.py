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
        self.lam = 0.98 # Loading coefficient for the first stage
        # Turbine parameters
        self.T_04 = 2275 # Turbine inlet temperature
        self.n_inf_t = 0.92 # Turbine polytropic efficiency
        self.phi = 0.78 # Lower bound on flow coefficient
        self.psi = 3.3 # Upper bound on loading coefficient
        self.y_h = 1.333 # Gamma in the hot section
        self.cp_c = 1.148 # Specific heat hot section, kJ/(kg*K)
        # Basic design parameters - constant tip radius assumption
        self.U_t1 = 350 # Max Tip speed, m/s
        self.C_a = 150 # Inlet axial velocity, m/s

        # Calcs for inlet conditions to the compressor
        # First, calculate fan outlet (assume intake has no loss)
        self.P_02 = self.pi_f*self.pa
        self.T_02 = self.Ta + self.Ta*((self.pi_f)**((self.y_c-1)/(self.y_c*self.n_inf_f)) - 1)

        # Static temperature, pressure, and density at the inlet
        T2 = self.T_02 - self.C_a**2/(2*self.cp_c*1000)
        self.p2 = self.P_02*(T2/self.T_02)**(self.y_c/(self.y_c-1))  
        self.rho2 = self.p2*1e5/(self.R*T2)

        # Tip radius and rotation rate
        self.rt = np.sqrt((self.mdot/(self.BPR+1))/(np.pi*self.rho2*self.C_a*(1-self.rr_rt_o**2)))
        self.N = self.U_t1/(2*np.pi*self.rt)

        # Root and mean radius at the inlet
        self.rr = self.rt*self.rr_rt_o
        self.rm = (self.rr+self.rt)/2

        # Velocity and Mach number relative to the inlet
        V1_t = np.sqrt(self.U_t1**2+self.C_a**2)
        self.M1_t = V1_t/(np.sqrt(self.y_c*self.R*T2))

        # Outlet calculations for the compressor
        self.P_03 = self.P_02*self.cpr
        self.T_03 = self.T_02*(self.cpr)**((self.y_c-1)/(self.y_c*self.n_inf_c))
        # Assume the exit velocity is axial and the same as the initial inlet velocity
        T3 = self.T_03 - self.C_a**2/(2*self.cp_c*1000)
        self.p3 = self.P_03*(T3/self.T_03)**(self.y_c/(self.y_c-1))
        self.rho3 = self.p3*1e5/(self.R*T3)

        A3 = self.mdot/(self.BPR+1)/(self.rho3*self.C_a)
        h3 = A3/(2*np.pi*self.rm)
        self.rt_3 = self.rm + (h3/2)
        self.rr_3 = self.rm - (h3/2)
        
        return 
    
    def update(self):
        return self.N, self.rt, self.rm, self.M1_t, self.rt_3, self.rr_3
    
    def fullturbo():

        return
    
    def fullcompressor(self):
        # Temperature rise through the compressor
        self.delt_Ts = self.T_03-self.T_02
        # Assuming axial velocity remains constant throughout each stage
        # Calculating mean blade speed
        self.U_m = 2*np.pi*self.rm*self.N
        # Calculating the estimated temperature rise per stage
        beta1 = np.arctan(self.U_m/self.C_a)
        V1 = self.C_a/np.cos(beta1)
        # Using de Haller criteria to find V2
        V2 = V1*self.haller
        beta2 = np.arccos(self.C_a/V2)
        T0s_est = self.lam*self.U_m*self.C_a*(np.tan(beta1)-np.tan(beta2))/(self.cp_c*1e3)
        stage_est = np.ceil(self.delt_Ts/T0s_est)        

        return T0s_est, stage_est, self.delt_Ts, self.T_03
    
    def fullturbine():

        return
    
    def compressorstage(self, tr):
        

        return
    
    def turbinestage():

        return
    
    def velocitytriangle():

        return

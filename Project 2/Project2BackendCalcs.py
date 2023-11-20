import numpy as np
import pandas as pd
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
        self.cp_h = 1.148 # Specific heat hot section, kJ/(kg*K)
        # Basic design parameters - constant tip radius assumption
        self.U_t1 = 350 # Max Tip speed, m/s
        self.C_a = 150 # Inlet axial velocity, m/s

        # Other parameters
        self.n_m = 0.99 # assumed mechanical efficiency of 99%
        self.delP_b = 0.06 # assumed burner pressure loss of 6%

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

        # Calculates inlet pressure to the turbine using burner efficiency
        self.P_04 = self.P_03*self.delP_b

        A3 = self.mdot/(self.BPR+1)/(self.rho3*self.C_a)
        h3 = A3/(2*np.pi*self.rm)
        self.rt_3 = self.rm + (h3/2)
        self.rr_3 = self.rm - (h3/2)
        
        return 
    
    def update(self):
        return self.N, self.rt, self.rm, self.rr, self.M1_t, self.rt_3, self.rr_3, self.P_02, self.T_02, self.P_03, self.T_03
        
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
        # Estimated delta T0 per stage
        T0s_est = self.lam*self.U_m*self.C_a*(np.tan(beta1)-np.tan(beta2))/(self.cp_c*1e3)
        # Estimated number of stages
        stage_est = np.ceil(self.delt_Ts/T0s_est)

        alpha1 = 0.0 # radians, inlet blade angle
        
        ################ Stage 1 ################
        # First stage calculation
        delta_T0 = 30 ## Desired temperature rise per stage for the first stage 
        delta_C_w = self.cp_c*1e3*delta_T0/(self.lam*self.U_m)
        # Whirl velocities
        C_w1 = self.C_a*np.tan(alpha1)
        C_w2 = delta_C_w + C_w1
        # Relative angles
        beta1 = np.arctan((self.U_m - C_w1)/self.C_a)
        beta2 = np.arctan((self.U_m - C_w2)/self.C_a)
        # Outlet angle of rotor
        alpha2 = np.arctan(C_w2/self.C_a)
        diffusion = np.cos(beta1)/np.cos(beta2)          
        # Organizing data for tables
        data = np.concatenate((np.rad2deg(np.array([alpha1,alpha2,beta1,beta2])), np.array([0.0,diffusion,0.0])))
        # Outlet stagnation pressure calculation of the first stage
        p031 = self.P_02*(1 + (self.n_inf_c*delta_T0)/self.T_02)**(self.y_c/(self.y_c-1))
        # Inlet tip Mach number
        T2 = self.T_02 - self.C_a**2/(2*self.cp_c*1000)
        M11t = self.C_a/np.cos(beta1)/np.sqrt(self.y_c*self.R*T2)
        # Outlet tip Mach number
        C21 = self.C_a/np.cos(alpha2)
        T021 = self.T_02 + delta_T0
        T21 = T021 - C21**2/(2*self.cp_c*1e3)
        M21t = self.C_a/np.cos(beta2)/np.sqrt(self.y_c*self.R*T21)
        # Approximate Degree of Reaction at the mean radius
        React = 1 - (C_w1+C_w2)/(2*self.U_m)

        meantable = pd.DataFrame(np.round([np.concatenate((data,np.array([T021]),np.array([p031/self.P_02]),np.array([p031]),np.array([M11t]),np.array([M21t]),np.array([React]),np.array([self.lam])))],3), columns=['alpha1','alpha2','beta1','beta2','alpha3','V2/V1','C3/C2','T02','P03/P01','P03','M1t','M2t','Reaction','Loading'])
        # meantable['C3/C2'][0] = 2.0

        ################ Stage 2 ################
        delta_T0 = 29 ## Desired temperature rise for the second stage
        React = 0.7 ## Degree of reaction for the second stage
        self.lam -= 0.03 ## Update loading coefficient
        # Calculate relative blade angles by solving system of eqs
        B1 = delta_T0*self.cp_c*1e3/(self.lam*self.U_m*self.C_a)
        B2 = React*2*self.U_m/self.C_a
        beta2 = np.arctan((B2-B1)/2)
        beta1 = np.arctan(B1+np.tan(beta2))
        # Calculate absolute blade angles
        alpha1 = np.arctan(self.U_m/self.C_a - np.tan(beta1))
        alpha2 = np.arctan(self.U_m/self.C_a - np.tan(beta2))
        # Calculation diffusion for de Haller criteria
        diffusion = np.cos(beta1)/np.cos(beta2)
        # Calculate outlet pressure
        p032 = p031*(1 + (self.n_inf_c*delta_T0)/T021)**(self.y_c/(self.y_c-1))
        # Calculate outlet temperature
        T022 = T021 + delta_T0
        # Calculate inlet Mach number
        C21 = self.C_a/np.cos(alpha1)
        T21 = T021 - C21**2/(2*self.cp_c*1e3)
        M12t = self.C_a/np.cos(beta1)/np.sqrt(self.y_c*self.R*T21)
        # Calculate outlet Mach number
        C22 = self.C_a/np.cos(alpha2)
        T22 = T022 - C22**2/(2*self.cp_c*1e3)
        M22t = self.C_a/np.cos(beta2)/np.sqrt(self.y_c*self.R*T22)
        # Data organization
        data = np.round(np.concatenate((np.rad2deg(np.array([alpha1,alpha2,beta1,beta2])), np.array([0.0,diffusion, 0.0,T022, p032/p031, p032, M12t, M22t, React, self.lam]))),3)
        # Update table
        meantable.loc[len(meantable)] = data

        ################ Stage 3 ################
        delta_T0 = 29 ## Desired temperature rise for the second stage
        React = 0.5 ## Degree of reaction for the second stage
        self.lam -= 0.03 ## Update loading coefficient
        # Calculate relative blade angles by solving system of eqs
        B1 = delta_T0*self.cp_c*1e3/(self.lam*self.U_m*self.C_a)
        B2 = React*2*self.U_m/self.C_a
        beta2 = np.arctan((B2-B1)/2)
        beta1 = np.arctan(B1+np.tan(beta2))
        # Calculate absolute blade angles
        alpha1 = np.arctan(self.U_m/self.C_a - np.tan(beta1))
        alpha2 = np.arctan(self.U_m/self.C_a - np.tan(beta2))
        # Calculation diffusion for de Haller criteria
        diffusion = np.cos(beta1)/np.cos(beta2)
        # Calculate outlet pressure
        p033 = p032*(1 + (self.n_inf_c*delta_T0)/T022)**(self.y_c/(self.y_c-1))
        # Calculate outlet temperature
        T023 = T022 + delta_T0
        # Calculate inlet Mach number
        C22 = self.C_a/np.cos(alpha1)
        T22 = T022 - C22**2/(2*self.cp_c*1e3)
        M13t = self.C_a/np.cos(beta1)/np.sqrt(self.y_c*self.R*T22)
        # Calculate outlet Mach number
        C23 = self.C_a/np.cos(alpha2)
        T23 = T023 - C23**2/(2*self.cp_c*1e3)
        M23t = self.C_a/np.cos(beta2)/np.sqrt(self.y_c*self.R*T23)
        # Data organization
        data = np.round(np.concatenate((np.rad2deg(np.array([alpha1,alpha2,beta1,beta2])), np.array([0.0,diffusion, 0.0,T023, p033/p032, p033, M13t, M23t, React, self.lam]))),3)
        # Update table
        meantable.loc[len(meantable)] = data

        test_p03 = p033
        test_T02 = T023

        ################ Stages 4-16 ################
        for i in range(3, int(stage_est)+1):
            if (self.lam-0.01) > 0.84:
                self.lam -= 0.01

            delta_T0 = 40.0
            React = 0.5

            data, test2_p03, test2_T02, diffusion = self.compressorstage(delta_T0, React, test_p03, test_T02)
            while diffusion < self.haller+0.01:
                delta_T0 -= 0.01
                data, test2_p03, test2_T02, diffusion = self.compressorstage(delta_T0, React, test_p03, test_T02)
            
            test_p03 = test2_p03
            test_T02 = test2_T02
            # Update table
            meantable.loc[len(meantable)] = data

        p01 = test_p03
        T01 = test_T02

        ################ Stage 17 ################
        delta_T0 = ((15.0*self.P_02/p01)**((self.y_c-1)/self.y_c)-1)*T01/self.n_inf_c # Desired temperature rise for the second stage
        React = 0.5 ## Degree of reaction for the second stage        
        # Calculate relative blade angles by solving system of eqs
        B1 = delta_T0*self.cp_c*1e3/(self.lam*self.U_m*self.C_a)
        B2 = React*2*self.U_m/self.C_a
        beta2 = np.arctan((B2-B1)/2)
        beta1 = np.arctan(B1+np.tan(beta2))
        # Calculate absolute blade angles
        alpha1 = np.arctan(self.U_m/self.C_a - np.tan(beta1))
        alpha2 = np.arctan(self.U_m/self.C_a - np.tan(beta2))
        # Calculation diffusion for de Haller criteria
        diffusion = np.cos(beta1)/np.cos(beta2)
        # Calculate outlet pressure
        p03 = p01*(1 + (self.n_inf_c*delta_T0)/T01)**(self.y_c/(self.y_c-1))
        # Calculate outlet temperature
        T02 = T01 + delta_T0
        # Calculate inlet Mach number
        C2 = self.C_a/np.cos(alpha1)
        T1 = T01 - C2**2/(2*self.cp_c*1e3)
        M1t = self.C_a/np.cos(beta1)/np.sqrt(self.y_c*self.R*T1)
        # Calculate outlet Mach number
        C23 = self.C_a/np.cos(alpha2)
        T2 = T02 - C23**2/(2*self.cp_c*1e3)
        M2t = self.C_a/np.cos(beta2)/np.sqrt(self.y_c*self.R*T2)
        # Data organization
        data = np.round(np.concatenate((np.rad2deg(np.array([alpha1,alpha2,beta1,beta2])), np.array([0.0,diffusion, 0.0,T02, p03/p01, p03, M1t, M2t, React, self.lam]))),3)
        # Update table
        meantable.loc[len(meantable)] = data

        # Iterate through table to organize stator outlet angles
        for i in range(0,len(meantable)-1):
            meantable['alpha3'][i] = meantable['alpha1'][i+1]

        # Iterate through table to calculate de Haller values for the stators
        for i in range(0,len(meantable)):
            alpha2 = np.deg2rad(meantable['alpha2'][i])
            alpha3 = np.deg2rad(meantable['alpha3'][i])
            s_diffusion = np.cos(alpha2)/np.cos(alpha3)
            meantable['C3/C2'][i] = np.round(s_diffusion,3)
        


        meantable.index = np.arange(1, len(meantable)+1)
        # meantable.index.name = 'Stage'
        # meantable.reset_index().to_string(index=False)
        return meantable
       
    def fullturbine(self):
        # Sets the rotational speed and mean blade speed
        N = self.N
        Um = 340 #assumed mean blade speed based on experience, m/s
        
        # Calculates the total temperature drop based on a work balance from the compressor (using assummed mech. eff.)
        dT0_turb = (self.cp_c*(self.T_03-self.T_02))/(self.cp_h*self.n_m)
        
        # Uses stage estimation based on constant drop (initial guess) over the stages
        T0s_est = 180 #estimate based on book example, K
        
        # Calculates the temp drop coeff.
        psi_turb = (2*self.cp_h*1e3*T0s_est)/(Um**2)
        # iteration to find the stage drop below the loading coefficient
        T0s_rev = T0s_est
        if np.round(psi_turb, 2) > 3.3:
            while np.round(psi_turb, 2) > 3.3:
                T0s_rev -= 5
                psi_turb = (2*self.cp_h*1e3*T0s_rev)/(Um**2)
        
        stage_est = np.ceil(dT0_turb/T0s_rev)

        ################ Stage 1 ################
        
        return dT0_turb, T0s_rev
    
    def compressorstage(self, delta_T0, React, p01, T01):
        # Calculate relative blade angles by solving system of eqs
        B1 = delta_T0*self.cp_c*1e3/(self.lam*self.U_m*self.C_a)
        B2 = React*2*self.U_m/self.C_a
        beta2 = np.arctan((B2-B1)/2)
        beta1 = np.arctan(B1+np.tan(beta2))
        # Calculate absolute blade angles
        alpha1 = np.arctan(self.U_m/self.C_a - np.tan(beta1))
        alpha2 = np.arctan(self.U_m/self.C_a - np.tan(beta2))
        # Calculation diffusion for de Haller criteria
        diffusion = np.cos(beta1)/np.cos(beta2)
        # Calculate outlet pressure
        p03 = p01*(1 + (self.n_inf_c*delta_T0)/T01)**(self.y_c/(self.y_c-1))
        # Calculate outlet temperature
        T02 = T01 + delta_T0
        # Calculate inlet Mach number
        C2 = self.C_a/np.cos(alpha1)
        T1 = T01 - C2**2/(2*self.cp_c*1e3)
        M1t = self.C_a/np.cos(beta1)/np.sqrt(self.y_c*self.R*T1)
        # Calculate outlet Mach number
        C23 = self.C_a/np.cos(alpha2)
        T2 = T02 - C23**2/(2*self.cp_c*1e3)
        M2t = self.C_a/np.cos(beta2)/np.sqrt(self.y_c*self.R*T2)
        # Data organization
        data = np.round(np.concatenate((np.rad2deg(np.array([alpha1,alpha2,beta1,beta2])), np.array([0.0,diffusion, 0.0,T02, p03/p01, p03, M1t, M2t, React, self.lam]))),3)
        
        return data, p03, T02, diffusion
    
    def velocitytriangle():

        return

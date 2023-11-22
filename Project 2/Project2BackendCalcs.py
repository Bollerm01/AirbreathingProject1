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
        self.U_t1 = 350 # Inlet Tip speed, m/s
        
        self.C_a = 150 # Inlet axial velocity, m/s

        # Other parameters
        self.n_m = 0.99 # assumed mechanical efficiency of 99%
        self.delP_b = 0.00 # assumed burner pressure loss of 0%
        self.n_b = 1.0 # assumed burner eff. of 100%

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

        # Calculates the fuel flow
        FAR = (self.cp_h*(self.T_04) - self.cp_c*(self.T_03))/(43100 - self.cp_h*self.T_03)
        m_fuel = FAR*self.mdoth
        self.mdoth = (self.mdot / (self.BPR + 1)) + m_fuel # Hot Section Flow, kg/s
        
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
        self.lam = self.loadingfactor(1)
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
        # Radius at the inlet to the rotor
        T11 = self.T_02 - self.C_a**2/(2*self.cp_c*1e3)
        p11 = self.P_02*(T11/self.T_02)**(self.y_c/(self.y_c-1))
        rho11 = p11*1e5/(T11*self.R)
        A11 = self.mdot/(self.BPR+1)/(rho11*self.C_a)
        h11 = A11/(2*np.pi*self.rm)
        rt11 = self.rm + (h11/2)
        rr11 = self.rm - (h11/2)
        rr_rt1 = rr11/rt11
        # Radius in between rotor and stator
        p21 = p031*(T21/T021)**(self.y_c/(self.y_c-1))
        rho21 = p21*1e5/(T21*self.R)
        A21 = self.mdot/(self.BPR+1)/(rho21*self.C_a)
        h21 = A21/(2*np.pi*self.rm)
        rt21 = self.rm + (h21/2)
        rr21 = self.rm - (h21/2)
        rr_rt = rr21/rt21

        sizingtable = pd.DataFrame(np.round([np.array([rr_rt1,rr_rt,0.0, h11,h21,0.0])],5),columns=['r_t 1','r_t 2','r_t 3','h1','h2','h3'])

        meantable = pd.DataFrame(np.round([np.concatenate((data,np.array([T021]),np.array([p031/self.P_02]),np.array([p031]),np.array([M11t]),np.array([M21t,C_w1,C_w2]),np.array([React]),np.array([self.lam])))],3), columns=['alpha1','alpha2','beta1','beta2','alpha3','V2/V1','C3/C2','T02','P03/P01','P03','M1t','M2t','Cw1','Cw2','Reaction','Loading'])
        # meantable['C3/C2'][0] = 2.0

        ################ Stage 2 ################
        delta_T0 = 29 ## Desired temperature rise for the second stage
        React = 0.7 ## Degree of reaction for the second stage
        self.lam = self.loadingfactor(2) ## Update loading coefficient
        # Calculate relative blade angles by solving system of eqs
        B1 = delta_T0*self.cp_c*1e3/(self.lam*self.U_m*self.C_a)
        B2 = React*2*self.U_m/self.C_a
        beta2 = np.arctan((B2-B1)/2)
        beta1 = np.arctan(B1+np.tan(beta2))
        # Calculate absolute blade angles
        alpha1 = np.arctan(self.U_m/self.C_a - np.tan(beta1))
        alpha2 = np.arctan(self.U_m/self.C_a - np.tan(beta2))
        # Calculate whirl velocities
        C_w1 = self.C_a*np.tan(alpha1)
        C_w2 = self.C_a*np.tan(alpha2)
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
        data = np.round(np.concatenate((np.rad2deg(np.array([alpha1,alpha2,beta1,beta2])), np.array([0.0,diffusion, 0.0,T022, p032/p031, p032, M12t, M22t,C_w1,C_w2, React, self.lam]))),3)
        # Radius at the inlet to the rotor
        T12 = T021 - C21**2/(2*self.cp_c*1e3)
        p12 = p031*(T12/T021)**(self.y_c/(self.y_c-1))
        rho12 = p12*1e5/(T12*self.R)
        A12 = self.mdot/(self.BPR+1)/(rho12*self.C_a)
        h12 = A12/(2*np.pi*self.rm)
        rt12 = self.rm + (h12/2)
        rr12 = self.rm - (h12/2)
        rr_rt1 = rr12/rt12
        # Radius in between rotor and stator
        p22 = p032*(T22/T022)**(self.y_c/(self.y_c-1))
        rho22 = p22*1e5/(T21*self.R)
        A22 = self.mdot/(self.BPR+1)/(rho22*self.C_a)
        h22 = A22/(2*np.pi*self.rm)
        rt22 = self.rm + (h22/2)
        rr22 = self.rm - (h22/2)
        rr_rt = rr22/rt22

        s_data = np.round(np.array([rr_rt1,rr_rt,0.0, h12,h22,0.0]),5)

        # Update table
        meantable.loc[len(meantable)] = data
        sizingtable.loc[len(sizingtable)] = s_data

        ################ Stage 3 ################
        delta_T0 = 29 ## Desired temperature rise for the second stage
        React = 0.5 ## Degree of reaction for the second stage
        self.lam = self.loadingfactor(3) ## Update loading coefficient
        # Calculate relative blade angles by solving system of eqs
        B1 = delta_T0*self.cp_c*1e3/(self.lam*self.U_m*self.C_a)
        B2 = React*2*self.U_m/self.C_a
        beta2 = np.arctan((B2-B1)/2)
        beta1 = np.arctan(B1+np.tan(beta2))
        # Calculate absolute blade angles
        alpha1 = np.arctan(self.U_m/self.C_a - np.tan(beta1))
        alpha2 = np.arctan(self.U_m/self.C_a - np.tan(beta2))
        # Calculate whirl velocities
        C_w1 = self.C_a*np.tan(alpha1)
        C_w2 = self.C_a*np.tan(alpha2)
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
        data = np.round(np.concatenate((np.rad2deg(np.array([alpha1,alpha2,beta1,beta2])), np.array([0.0,diffusion, 0.0,T023, p033/p032, p033, M13t, M23t,C_w1,C_w2, React, self.lam]))),3)
        # Radius at the inlet to the rotor
        p13 = p032*(T22/T022)**(self.y_c/(self.y_c-1))
        rho13 = p13*1e5/(T22*self.R)
        A13 = self.mdot/(self.BPR+1)/(rho13*self.C_a)
        h13 = A13/(2*np.pi*self.rm)
        rt13 = self.rm + (h13/2)
        rr13 = self.rm - (h13/2)
        rr_rt1 = rr13/rt13
        # Radius in between rotor and stator
        p23 = p033*(T23/T023)**(self.y_c/(self.y_c-1))
        rho23 = p23*1e5/(T23*self.R)
        A23 = self.mdot/(self.BPR+1)/(rho23*self.C_a)
        h23 = A23/(2*np.pi*self.rm)
        rt23 = self.rm + (h23/2)
        rr23 = self.rm - (h23/2)
        rr_rt = rr23/rt23

        s_data = np.round(np.array([rr_rt1,rr_rt,0.0,h13,h23,0.0]),5)

        # Update table
        meantable.loc[len(meantable)] = data
        sizingtable.loc[len(sizingtable)] = s_data

        test_p03 = p033
        test_T02 = T023

        test3_p03 = test_p03

        final_P03 = self.cpr*self.P_02
        ns = 4

        ################ Intermediate Stages ################
        while test3_p03 < final_P03:
            self.lam = self.loadingfactor(ns)
            ns += 1

            delta_T0 = 100.0
            React = 0.5
            data, s_data, test2_p03, test2_T02, diffusion = self.compressorstage(delta_T0, React, test_p03, test_T02)
            while diffusion < self.haller+0.01:
                delta_T0 -= 0.01
                data, s_data, test2_p03, test2_T02, diffusion = self.compressorstage(delta_T0, React, test_p03, test_T02)

            test_p03 = test2_p03
            test_T02 = test2_T02
            # Update table
            meantable.loc[len(meantable)] = data
            sizingtable.loc[len(sizingtable)] = s_data

            # Calculate next estimated max
            delta_T0 = 100.0
            data, s_data, test2_p03, test2_T02, diffusion = self.compressorstage(delta_T0, React, test_p03, test_T02)
            while diffusion < self.haller+0.01:
                delta_T0 -= 0.01
                data, s_data, test2_p03, test2_T02, diffusion = self.compressorstage(delta_T0, React, test_p03, test_T02)
            
            test3_p03 = test2_p03
        
        p01 = test_p03
        T01 = test_T02

        ################ Last Stage ################
        delta_T0 = ((15.0*self.P_02/p01)**((self.y_c-1)/self.y_c)-1)*T01/self.n_inf_c # Desired temperature rise for the second stage
        React = 0.5 ## Degree of reaction for the second stage   
        self.lam = self.loadingfactor(15)     
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
        data = np.round(np.concatenate((np.rad2deg(np.array([alpha1,alpha2,beta1,beta2])), np.array([0.0,diffusion, 0.0,T02, p03/p01, p03, M1t, M2t,C_w1,C_w2, React, self.lam]))),3)
         # Radius at the inlet to the rotor
        C1 = self.C_a/np.cos(alpha1)
        T1 = T01 - C1**2/(2*self.cp_c*1e3)
        p1 = p01*(T1/T01)**(self.y_c/(self.y_c-1))
        rho1 = p1*1e5/(T1*self.R)
        A1 = self.mdot/(self.BPR+1)/(rho1*self.C_a)
        h1 = A1/(2*np.pi*self.rm)
        rt1 = self.rm + (h1/2)
        rr1 = self.rm - (h1/2)
        rr_rt1 = rr1/rt1
        # Radius in between rotor and stator
        p2 = p03*(T2/T02)**(self.y_c/(self.y_c-1))
        rho2 = p2*1e5/(T2*self.R)
        A2 = self.mdot/(self.BPR+1)/(rho2*self.C_a)
        h2 = A2/(2*np.pi*self.rm)
        rt2 = self.rm + (h2/2)
        rr2 = self.rm - (h2/2)
        rr_rt = rr2/rt2
        # Radius outlet of the stator
        T3 = T02 - self.C_a**2/(2*self.cp_c*1e3)
        p3 = p03*(T3/T02)**(self.y_c/(self.y_c-1))
        rho3 = p3*1e5/(T3*self.R)
        A3 = self.mdot/(self.BPR+1)/(rho3*self.C_a)
        h3 = A3/(2*np.pi*self.rm)
        rt3 = self.rm + (h3/2)
        rr3 = self.rm - (h3/2)
        rr_rt3 = rr3/rt3

        s_data = np.round(np.array([rr_rt1,rr_rt,rr_rt3,h1,h2,h3]),5)

        # Update table
        meantable.loc[len(meantable)] = data
        sizingtable.loc[len(sizingtable)] = s_data

        # Iterate through table to organize stator outlet angles
        for i in range(0,len(meantable)-1):
            meantable['alpha3'][i] = meantable['alpha1'][i+1]
            sizingtable['r_t 3'][i] = sizingtable['r_t 1'][i+1]
            sizingtable['h3'][i] = sizingtable['h1'][i+1]

        # Iterate through table to calculate de Haller values for the stators
        for i in range(0,len(meantable)):
            alpha2 = np.deg2rad(meantable['alpha2'][i])
            alpha3 = np.deg2rad(meantable['alpha3'][i])
            s_diffusion = np.cos(alpha2)/np.cos(alpha3)
            meantable['C3/C2'][i] = np.round(s_diffusion,3)
        

        sizingtable.index = np.arange(1, len(sizingtable)+1)
        meantable.index = np.arange(1, len(meantable)+1)
        # meantable.index.name = 'Stage'
        # meantable.reset_index().to_string(index=False)
        return meantable, sizingtable
       



    def fullturbine(self):
        ################ Preliminary Sizing ################
        # Sets the rotational speed and mean blade speed
        N = self.N
        Um = 340 #assumed mean blade speed based on experience, m/s
        lamdaN = 0.05 # assumed standard value
        
        # Calculates the total temperature drop based on a work balance from the compressor (using assummed mech. eff.)
        dT0_turb = (self.cp_c*(self.T_03-self.T_02))/(self.cp_h*self.n_m)
        
        # Uses stage estimation based on constant drop (initial guess) over the stages
        T0s_est = 1000 # K
        
        # Calculates the temp drop coeff.
        psi_turb = (2*self.cp_h*1e3*T0s_est)/(Um**2)
        # iteration to find the stage drop below the loading coefficient
        T0s_rev = T0s_est
        if np.round(psi_turb, 2) > 3.3:
            while np.round(psi_turb, 2) > 3.3:
                T0s_rev -= 5
                psi_turb = (2*self.cp_h*1e3*T0s_rev)/(Um**2)
        
        stage_est = (dT0_turb/T0s_rev)

        ################ Stage 1 ################
        # Assumptions for first stage 
        alpha1 = 0.0
        alpha3 = 0.0

        stg1Angles, stg1meanLambda, stg1meanV, stg1Areas, stg1RTrat, stg1Heights = self.turbinestage(self, alpha1, alpha3, Um, T0s_rev, psi_turb)

        '''# Calculates the B3 and deg. of reaction
        beta3 = np.arctan(np.tan(alpha3) + (1/self.phi))
        Lambda = (2*self.phi*np.tan(beta3)- (psi_turb/2))/2 
        # iterates to find a suitable degree of reaction
        if np.round(Lambda, 3) > 0.420:
            while np.round(Lambda, 3) > 0.420:
                alpha3 += np.deg2rad(5)
                beta3 = np.arctan(np.tan(alpha3) + (1/self.phi))
                Lambda = (2*self.phi*np.tan(beta3)- (psi_turb/2))/2
        
        # Calculates B2 and a2
        beta2 = np.arctan((1/(2*self.phi))*(psi_turb/2 - 2*Lambda))
        alpha2 = np.arctan(np.tan(beta2) + (1/self.phi))

        # Calculates Ca2 and C2
        Ca2 = Um*self.phi
        C2 = Ca2 / np.cos(alpha2)

        # Calculates the isentropic T2', T2 and p2
        T2 = self.T_04 - (C2**2/(2*self.cp_h*1e3))
        T2prime = T2- lamdaN*(C2**2/(2*self.cp_h*1e3))
        stagStaticRat = (self.T_04/T2prime)**(self.y_h/(self.y_h-1))
        CritPrat = 1.853
        if stagStaticRat > CritPrat:
            print('ask the child if he or she is choking')

        P2 = self.P_04/stagStaticRat

        # Calculate the rho2 and A2
        rho2 = (P2*100)/(0.287*T2)      
        A2 = self.mdot/(rho2*Ca2)
        
        # Calculates the Ca1, using Ca2 = Ca3 and C1 = C3
        Ca3 = Ca2
        C3 = Ca3/np.cos(alpha3)
        C1 = C3
        Ca1 = C1 #since alpha1 = 0.0

        # Calculates rho1 and A1
        T1 = self.T_04 - (C1**2/(2*self.cp_h*1e3))
        P1 = self.P_04*(T1/self.T_04)**(self.y_h/(self.y_h-1))
        rho1 = (P1*100)/(0.287*T1)
        A1 = self.mdot/(rho1*Ca1)

        # Calculates the outlet (station 3) conditions
        T_03 = self.T_04 - T0s_rev
        T3 = T_03 - (C3**2)/(2*self.cp_h*1e3)

        # Calculates the P03 from the isentropic eff. (using P04 = P01 and T04 = T01) and P3 from isentropic
        P_03 = self.P_04*(1 - (T0s_rev/(self.n_inf_t)))**(self.y_h/(self.y_h-1))
        P3 = P_03*(T3/T_03)**(self.y_h/(self.y_h-1))

        # Calculates the rho3, A3
        rho3 = (P3*100)/(0.287*T3)
        A3 = self.mdot/(rho3*Ca3)

        # Rework to put at beginning for sizing 
        # Size at inlet and outlet of turbine 
        # Calculate the mean radius, use for calcs to get Um
        # Calculate rm
        rm = Um/(2*np.pi*N)

        # Calculates the h1-h3
        h1 = (N/Um)*A1
        h2 = (N/Um)*A2
        h3 = (N/Um)*A3

        # Calculates the rt/rr
        rtRat1 = (rm + (h1/2))/(rm - (h1/2))
        rtRat2 = (rm + (h2/2))/(rm - (h2/2))
        rtRat3 = (rm + (h3/2))/(rm - (h3/2))

        # # Calculates the Yn (T02 = T04)
        # P_02 = P2 / ((T2/self.T_04)**(self.y_h/(self.y_h - 1))) 
        # print(self.P_04)
        # print(P_02)
        # print(P2)
        # Yn = (self.P_04 - P_02)/(P_02 - P2) '''
        
        # Uses stage estimation based on constant drop (initial guess) over the stages
        T0_remaining = dT0_turb - T0s_rev # K
        T0s_est = T0_remaining

        # Calculates the temp drop coeff.
        psi_turb = (2*self.cp_h*1e3*T0s_est)/(Um**2)
        # iteration to find the stage drop below the loading coefficient
        T0s_rev = T0s_est
        if np.round(psi_turb, 2) > 3.3:
            while np.round(psi_turb, 2) > 3.3:
                T0s_rev -= 5
                psi_turb = (2*self.cp_h*1e3*T0s_rev)/(Um**2)
        
        ################ Stage 2 ################
        # Assumptions for second stage 
        alpha1 = alpha3 # new alpha1 is the previous stage alpha3
        alpha3 = 0.0

        stg2Angles, stg2meanLambda, stg2meanV, stg2Areas, stg2RTrat, stg2Heights = self.turbinestage(self, alpha1, alpha3, Um, T0s_rev, psi_turb)

        ## START HERE WITH CHECKING VALUES AFTER
        return dT0_turb, T0s_rev, stage_est
    


    def compressorstage(self, delta_T0, React, p01, T01):
        # Calculate relative blade angles by solving system of eqs
        B1 = delta_T0*self.cp_c*1e3/(self.lam*self.U_m*self.C_a)
        B2 = React*2*self.U_m/self.C_a
        beta2 = np.arctan((B2-B1)/2)
        beta1 = np.arctan(B1+np.tan(beta2))
        # Calculate absolute blade angles
        alpha1 = np.arctan(self.U_m/self.C_a - np.tan(beta1))
        alpha2 = np.arctan(self.U_m/self.C_a - np.tan(beta2))
        # Calculate whirl velocities
        C_w1 = self.C_a*np.tan(alpha1)
        C_w2 = self.C_a*np.tan(alpha2)
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
        data = np.round(np.concatenate((np.rad2deg(np.array([alpha1,alpha2,beta1,beta2])), np.array([0.0,diffusion, 0.0,T02, p03/p01, p03, M1t, M2t,C_w1,C_w2, React, self.lam]))),3)
        # Radius at the inlet to the rotor
        C1 = self.C_a/np.cos(alpha1)
        T1 = T01 - C1**2/(2*self.cp_c*1e3)
        p1 = p01*(T1/T01)**(self.y_c/(self.y_c-1))
        rho1 = p1*1e5/(T1*self.R)
        A1 = self.mdot/(self.BPR+1)/(rho1*self.C_a)
        h1 = A1/(2*np.pi*self.rm)
        rt1 = self.rm + (h1/2)
        rr1 = self.rm - (h1/2)
        rr_rt1 = rr1/rt1
        # Radius in between rotor and stator
        p2 = p03*(T2/T02)**(self.y_c/(self.y_c-1))
        rho2 = p2*1e5/(T2*self.R)
        A2 = self.mdot/(self.BPR+1)/(rho2*self.C_a)
        h2 = A2/(2*np.pi*self.rm)
        rt2 = self.rm + (h2/2)
        rr2 = self.rm - (h2/2)
        rr_rt = rr2/rt2

        s_data = np.round(np.array([rr_rt1,rr_rt,0.0,h1,h2,0.0]),5)
        
        return data, s_data, p03, T02, diffusion
    
    def turbinestage(self, alpha1, alpha3, Um, T0s_rev, psi_turb):
        '''
            Inputs: alpha1 = stage inlet angle
                    alpha3 = stage exit angle
                    Um = mean blade speed
                    T0s_rev = revised stage temp drop
                    psi_turb = turbine temp drop coeff.
            Outputs:
                    Mean Angles = [alpha1, alpha2, alpha3, B2, B3]
                    Mean Deg. of React. = lambda
                    Mean Velocities = [C1, Ca1, C2, Ca2, V2, C3, Ca3, V3]
                    Pressure Ratio = P_02/P_01
                    Areas = [A1, A2, A3]
                    Root/Tip = [rt1, rt2, rt3]
                    h = [h1, h2, h3]                
        '''
        N = 340 #m/s, based on stress experience
        lambdaN = 0.05 # nozzle loss coefficient based on experience

        # Calculates the B3 and deg. of reaction
        beta3 = np.arctan(np.tan(alpha3) + (1/self.phi))
        Lambda = (2*self.phi*np.tan(beta3)- (psi_turb/2))/2 
        # iterates to find a suitable degree of reaction
        if np.round(Lambda, 3) > 0.420:
            while np.round(Lambda, 3) > 0.420:
                alpha3 += np.deg2rad(5)
                beta3 = np.arctan(np.tan(alpha3) + (1/self.phi))
                Lambda = (2*self.phi*np.tan(beta3)- (psi_turb/2))/2
        
        # Calculates B2 and a2
        beta2 = np.arctan((1/(2*self.phi))*(psi_turb/2 - 2*Lambda))
        alpha2 = np.arctan(np.tan(beta2) + (1/self.phi))

        # Calculates Ca2 and C2
        Ca2 = Um*self.phi
        C2 = Ca2 / np.cos(alpha2)

        # Calculates the isentropic T2', T2 and p2
        T2 = self.T_04 - (C2**2/(2*self.cp_h*1e3))
        T2prime = T2- lambdaN*(C2**2/(2*self.cp_h*1e3))
        stagStaticRat = (self.T_04/T2prime)**(self.y_h/(self.y_h-1))
        CritPrat = 1.853
        if stagStaticRat > CritPrat:
            print('ask the child if he or she is choking')

        P2 = self.T_04/stagStaticRat

        # Calculate the rho2 and A2
        rho2 = (P2*100)/(0.287*T2)      
        A2 = self.mdoth/(rho2*Ca2)
        
        # Calculates the Ca1, using Ca2 = Ca3 and C1 = C3
        Ca3 = Ca2
        C3 = Ca3/np.cos(alpha3)
        C1 = C3
        Ca1 = C1 #since alpha1 = 0.0

        # Calculates rho1 and A1
        T1 = self.T_04 - (C1**2/(2*self.cp_h*1e3))
        P1 = self.P_04*(T1/self.T_04)**(self.y_h/(self.y_h-1))
        rho1 = (P1*100)/(0.287*T1)
        A1 = self.mdoth/(rho1*Ca1)

        # Calculates the outlet (station 3) conditions
        T_03 = self.T_04 - T0s_rev
        T3 = T_03 - (C3**2)/(2*self.cp_h*1e3)

        # Calculates the P03 from the isentropic eff. (using P04 = P01 and T04 = T01) and P3 from isentropic
        P_03 = self.P_04*(1 - (T0s_rev/(self.n_inf_t)))**(self.y_h/(self.y_h-1))
        P3 = P_03*(T3/T_03)**(self.y_h/(self.y_h-1))

        # Calculates the rho3, A3
        rho3 = (P3*100)/(0.287*T3)
        A3 = self.mdoth/(rho3*Ca3)

        # Mean radius calculation
        rm = Um/(2*np.pi*N)

        # Calculates the h1-h3
        h1 = (N/Um)*A1
        h2 = (N/Um)*A2
        h3 = (N/Um)*A3

        # Calculates the rt/rr
        rtRat1 = (rm + (h1/2))/(rm - (h1/2))
        rtRat2 = (rm + (h2/2))/(rm - (h2/2))
        rtRat3 = (rm + (h3/2))/(rm - (h3/2))

        # Calculates the V2 and V3
        V2 = Ca2/np.cos(beta2)
        V3 = ((Um+(C3*np.cos(alpha3)))**2 + Ca3**2)**0.5
        
        # Organizes return statements
        Angles = np.array([alpha1, alpha2, alpha3, beta2, beta3])
        meanLambda = Lambda
        meanV = np.array([C1, Ca1, C2, Ca2, V2, C3, Ca3, V3])
        Areas = np.array([A1, A2, A3])
        RTrat = np.array([rtRat1, rtRat2, rtRat3])
        Heights = np.array([h1, h2, h3])
        
        return Angles, meanLambda, meanV, Areas, RTrat, Heights
    
    def velocitytriangle():

        return
    
    def comp_root_tip(self, meantable,sizingtable,rm,N):
        for i in range(0,len(meantable)):
            h1 = sizingtable['h1'][i]
            h2 = sizingtable['h2'][i]
            h3 = sizingtable['h3'][i]
            rt1 = rm + (h1/2); rr1 = rm - (h1/2)
            rt2 = rm + (h2/2); rr2 = rm - (h2/2)
            rt3 = rm + (h3/2); rr3 = rm - (h3/2)
            
            

        return
    
    def loadingfactor(self,stage):
        # Loading factor calculation curve fit
        a=0.847936849639078; b=0.15697659830655492; c=-0.2412700053204237
        return a + b*np.exp(c*stage)

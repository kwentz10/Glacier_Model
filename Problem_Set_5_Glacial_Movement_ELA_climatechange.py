#Glacier Model 
#Created February 16 2016

from __future__ import division
import numpy as np
from matplotlib import pyplot as plt
from scipy import interpolate


class Glacier:

    def start_conds(self, t_min, t_max, dt, x_min, x_max, dx, zmax, s, h0):
        self.timestep = np.round(np.arange(t_min, t_max+dt, dt),decimals=2)
        self.spacestep = np.arange(x_min, x_max+dx, dx) 
        self.zb = zmax-s*(self.spacestep)
        self.h=np.zeros(len(self.zb))+h0
        self.z=self.zb+self.h
        self.fb = plt.figure(figsize=(10,8)) #figure blueprint
        self.figure1 = plt.subplot(211)
        self.figure2 = plt.subplot(212)
        #inital plot of interpolated ELA line (using temperature anomalies from Marcott et al. 2013 "A reconstruction of regional and global temperature for the past 11,300 years & IPCC WI 2013)
        self.ELA_timestep=[-11000,-10500,-10000,-9500,-9000,-8500,-8000,-7500,
                           -7000,-6500,-6000,-5500,-5000,-4500,-4000,-3500,
                           -3000,-2500,-2000,-1500,-1000,-500,0,100]
        self.ELA_z=[-30, 30, 45, 60, 60, 45, 52.5, 60, 75, 52.5, 55.5, 58.5, 
                    60, 57, 45, 30, 15, 1, -6, 0, -12, -45, 90, 300]
        self.ELA_z=np.array(self.ELA_z)+ELA_mean
        self.ELA_f=interpolate.interp1d(self.ELA_timestep,self.ELA_z, kind='linear') #ELA interpolation function
        self.ELA=self.ELA_f(self.timestep) #interpolate
        self.figure2.plot(self.ELA_timestep, self.ELA_z,'o',self.timestep, self.ELA, color='k')
        plt.ion()
        plt.show()

    def ela_func(self): #determine ELA at each timestep for figure 2
       self.ELA = self.ELA_f(self.timestep[t])
       return self.ELA
       
    def eff_stress_func(self,density_ice,g): #determine effective stress
        self.Ne=0.3*(self.h*density_ice*g)
        return self.Ne
        
    def bas_shear_stress_func(self, density_ice,g,s): #determine basal shear stress
        self.bss = density_ice*g*self.h*np.sin(np.deg2rad(s))   
        return self.bss
    
    def eros_func(self, C1, C2, density_ice,g,s): #combine effective stress, basal shear stress, C1, C2, g, ice density, and s to determine erosion
        self.Eros = C1*(glacier_model.bas_shear_stress_func(density_ice,g,s)**2)*C2/(glacier_model.eff_stress_func(density_ice,g))
        return self.Eros

    def b_func (self, gamma):
        self.b = gamma*(self.z-glacier_model.ela_func())
        return self.b

    def slope_func (self, dx):
        self.slope = np.diff(self.z)/dx
        return self.slope

    def havg_func(self): #average thickness in "edge" cells
        self.havg=(self.h[0:-1]+self.h[1:])/2
        return self.havg

    def q_func(self, A, density_ice, g, dx):
        self.q=A*((density_ice*g*glacier_model.slope_func(dx))**3)*(glacier_model.havg_func()**5)/5
        self.q=np.insert(self.q,0,0)
        self.q=np.append(self.q,0)
        return self.q

    def dqdx_func(self, dx):
        self.dqdx=-(np.diff(glacier_model.q_func(A, density_ice, g, dx))/dx)
        return self.dqdx

    def run (self, gamma, dt, dx, C1, C2, density_ice,g,s):
        global t  #make global so that it can be used in ELA_func as index
        for t in range(len(self.timestep)):
            self.dqdx=glacier_model.dqdx_func(dx)
            self.dhdt=glacier_model.b_func(gamma)-self.dqdx
            self.h+=self.dhdt*dt #add glacier thickness over time
            self.zb=self.zb-(glacier_model.eros_func(C1,C2,density_ice,g,s)*dt) #erosion function
            for i in range(0, len(self.h)): #make sure bottom limit of z does not go below zb
                if self.h[i]<0:
                    self.h[i]=0
            self.z=self.h+self.zb
            if self.timestep[t] % 10 ==0:
                #figure 1 Glacier accumulation       
                self.figure1.clear()
                self.figure1.set_title('Glacier Accumulation & Ablation Over Time')
                self.figure1.set_xlabel('Distance (m)')
                self.figure1.set_ylabel('Elevation (m)')
                self.figure1.set_ylim(2200,3000) #make sure y axis doesn't change
                self.figure1.set_xlim(0, 5000) #make sure x axis doesn't change
                self.figure1.text(2500, 2750, 'Time [yrs]: %d\nELA=%d m' % (self.timestep[t], self.ELA)) #plot text of time
                self.figure1.plot(self.spacestep, self.z, label='Glacier Height', color='b')  
                self.figure1.plot(self.spacestep, self.zb, label='Bedrock Height', color='k')
                self.figure1.plot(self.spacestep, (np.zeros(len(self.spacestep))+self.ELA),  '--r', label='ELA')
                self.figure1.fill_between(self.spacestep, self.zb, self.z, color ='b', interpolate=True)
                self.figure1.fill_between(self.spacestep, 0, self.zb, color='k', interpolate=True)
                self.figure1.legend()
                #figure 2 ELA line
                self.figure2.set_title('ELA')
                self.figure2.set_xlabel('Time Before Present (yrs)')
                self.figure2.set_ylabel('Elevation(m)')
                h=self.figure2.plot(self.timestep[t], self.ELA, 'ro')
                plt.pause(0.00001)
                l=h.pop(0) #remove last plot entry
                l.remove()                  
        plt.ioff()    

    def finalize(self):
        self.fb.savefig('kwentz_pset5.png', bbox_inches='tight')
 
if __name__== "__main__":
        t_min=-11000 #years
        t_max=100.0 # years
        dt=0.1 #years
        x_min=0.0 #m
        x_max=10000.0 #m
        dx=100.0 #m
        zmax=2750.0 #m
        s=0.15 #slope of glacier elevation linear function--degrees
        h0=100.0 #m
        gamma=0.01 #m/yr
        ELA_mean=2620.0 #m
        A=2.1*10**(-16) #Pa-3yr-1
        density_ice=917.0 #kgm-3
        g=9.81 #ms-2
        C1=0.0012 #mPa-1yr-1 Sliding Coefficient     
        C2=0.003 #erosion coefficient (I multiplied it by 30 to exagerrate erosion)
        
        glacier_model=Glacier()

        #initialize
        glacier_model.start_conds(t_min, t_max, dt, x_min, x_max, dx, zmax, s, h0)
        #run
        glacier_model.run(gamma, dt, dx, C1, C2, s, density_ice,g)
        #finalize
        glacier_model.finalize()

   


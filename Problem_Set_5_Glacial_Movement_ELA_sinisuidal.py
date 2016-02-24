#Glacier Model 
#Created February 16 2016

from __future__ import division
import numpy as np
from matplotlib import pyplot as plt


class Glacier:

    def start_conds(self, t_min, t_max, dt, x_min, x_max, dx, zmax, s, h0):
        self.timestep = np.arange(t_min, t_max+dt, dt)
        self.spacestep = np.arange(x_min, x_max, dx) 
        self.zb = zmax-s*(self.spacestep)
        self.h=np.zeros(len(self.zb))+h0
        self.z=self.zb+self.h
        self.fb = plt.figure(figsize=(10,8)) #figure blueprint
        self.figure1 = plt.subplot(211)
        self.figure2 = plt.subplot(212)
        #inital plot of ELA line
        self.ELA=ELA_mean+ (Amp*np.sin(2*np.pi*(self.timestep/P)))
        self.figure2.plot(self.timestep, self.ELA, color='k')
        plt.ion()
        plt.show()

    def ELA_func(self, Amp, P, ELA_mean):
       self.ELA = ELA_mean + (Amp*np.sin(2*np.pi*(self.timestep[t]/P)))
       return self.ELA

    def b_func (self, gamma, Amp, P, ELA_mean):
        self.b = gamma*(self.z-glacier_model.ELA_func(Amp, P, ELA_mean))
        return self.b

    def slope_func (self, dx):
        self.slope = np.diff(self.z)/dx
        return self.slope

    def havg_func(self):
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

    def run (self, gamma, Amp, P, ELA_mean, dt, dx):
        global t  #make global so that it can be used in ELA_func as index
        for t in range(len(self.timestep)):
            self.dqdx=glacier_model.dqdx_func(dx)
            self.dhdt=glacier_model.b_func(gamma, Amp, P, ELA_mean)-self.dqdx
            self.h+=self.dhdt*dt
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
                self.figure1.set_ylim(2000,4000) #make sure y axis doesn't change
                self.figure1.set_xlim(0, 15000) #make sure x axis doesn't change
                self.figure1.text(12500, 2900, 'Time [yrs]: %d\nELA=%d m' % (self.timestep[t], self.ELA))
                self.figure1.plot(self.spacestep, self.z, label='Glacier Height', color='b')  
                self.figure1.plot(self.spacestep, self.zb, label='Bedrock Height', color='k')
                self.figure1.plot(self.spacestep, (np.zeros(len(self.spacestep))+self.ELA),  '--r', label='ELA')
                self.figure1.fill_between(self.spacestep, self.zb, self.z, color ='b', interpolate=True)
                self.figure1.fill_between(self.spacestep, 0, self.zb, color='k', interpolate=True)
                self.figure1.legend()
                #figure 2 ELA line
                self.figure2.set_title('ELA Oscillation')
                self.figure2.set_xlabel('Time (yrs)')
                self.figure2.set_ylabel('Elevation(m)')
                h=self.figure2.plot(self.timestep[t], self.ELA, 'ro')
                plt.pause(0.00001)
                l=h.pop(0) #remove last plot entry
                l.remove()                  
        plt.ioff()    

    def finalize(self):
        self.fb.savefig('kwentz_pset5.png', bbox_inches='tight')
 
if __name__== "__main__":
        t_min=0.0 #years
        t_max=400.0 # years
        dt=0.001 #years
        x_min=0.0 #m
        x_max=16000.0 #m
        dx=100.0 #m
        zmax=2550 #m
        s=0.03 #slope of glacier elevation linear function
        h0=0.0 #m
        gamma=0.01 #m/yr
        ELA_mean= 2450 #m
        Amp=25 #m
        P=50 #years
        A=(2.1*10**(-16)) #Pa-3yr-1
        density_ice=917 #kgm-3
        g=9.81 #ms-2
        
        glacier_model=Glacier()

        #initiallize
        glacier_model.start_conds(t_min, t_max, dt, x_min, x_max, dx, zmax, s, h0)
        #run
        glacier_model.run(gamma, Amp, P, ELA_mean, dt, dx)
        #finalize
        glacier_model.finalize()

   


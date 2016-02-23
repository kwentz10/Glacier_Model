#Glacier Model 
#Created February 16 2016

from __future__ import division
import numpy as np
from matplotlib import pyplot as plt


class Glacier:

    def start_conds(self, t_min, t_max, dt, x_min, x_max, dx, zmax, s, h0):
        self.timestep = np.arange(t_min, t_max+dt, dt).tolist() 
        self.spacestep = np.arange(x_min, x_max, dx) 
        self.zb = zmax-s*(self.spacestep)
        self.h=np.zeros(len(self.zb))+h0
        self.z=self.zb+self.h
        self.fb = plt.figure(figsize=(10,8)) #figure blueprint
        self.figure = plt.subplot()
        plt.ion()
        plt.show()

    def b_func (self, gamma, ELA):
        self.b = gamma*(self.z-ELA)
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

    def run (self, gamma, ELA, dt, dx):
        for t in range(len(self.timestep)):
            self.dqdx=glacier_model.dqdx_func(dx)
            self.dhdt=glacier_model.b_func(gamma, ELA)-self.dqdx
            self.h+=self.dhdt*dt
            self.z=self.h+self.zb
            self.z=self.z
            for i in range(0, len(self.z)): #make sure bottom limit of z does not go below zb
                if self.z[i]<self.zb[i]:
                    self.z[i]=self.zb[i]
            if self.timestep[t] % 1 ==0:
                self.figure.clear()
                plt.title('Glacier Accumulation & Ablation Over Time')
                plt.xlabel('Distance (m)')
                plt.ylabel('Elevation (m)')
                self.figure.set_ylim(2200, 3000) #make sure y axis doesn't change
                self.figure.set_xlim(0, 9000) #make sure x axis doesn't change
                plt.text(7500, 2800, 'Time [yrs]: %d\nELA=%d m' % (self.timestep[t], ELA))
                self.figure.plot(self.spacestep, self.z, label='Glacier Height', color='b')  
                self.figure.plot(self.spacestep, self.zb, label='Bedrock Height', color='k')
                self.figure.plot(self.spacestep, (np.zeros(len(self.spacestep))+ELA),  '--r', label='ELA')
                self.figure.fill_between(self.spacestep, self.zb, self.z, color ='b', interpolate=True)
                self.figure.fill_between(self.spacestep, 0, self.zb, color='k', interpolate=True)
                self.figure.legend()
                plt.pause(0.00001)
##            print self.z[0], self.z[len(self.z)-1]
        plt.ioff()    

    def finalize(self):
        self.fb.savefig('kwentz_pset5_linear.png', bbox_inches='tight')
 
if __name__== "__main__":
        t_min=0.0 #years
        t_max=150.0 # years
        dt=0.001 #years
        x_min=0.0 #m
        x_max=10000.0 #m
        dx=1000 #m
        zmax=2550 #m
        s=0.03 #slope of glacier elevation linear function
        h0=50.0 #m
        gamma=0.01 #m/yr
        ELA= 2500 #m
        A=(2.1*10**(-16)) #Pa-3yr-1
        density_ice=917 #kgm-3
        g=9.81 #ms-2
        
        glacier_model=Glacier()

        #initiallize
        glacier_model.start_conds(t_min, t_max, dt, x_min, x_max, dx, zmax, s, h0)
        #run
        glacier_model.run(gamma, ELA, dt, dx)
        #finalize
        glacier_model.finalize()

   


# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 13:57:36 2020

@author: Tobias
"""

import geometry
import em_env 
import simulation
import numpy as np
import sparameter
import skrf as rf
import random
import scipy.io as sio

if __name__=='__main__':    
    dists=[0.73*1e-3]#np.linspace(0.56,0.96,3)*1e-3
    for d in dists:
        cable_length=3
        #twistlength=0.065
        random.seed(0)
        twistlength=[]
        while sum(twistlength) < cable_length:
            twistlength.append(0.065+(0.065*0.3*random.uniform(-1,1))) 
        twistlength=tuple(twistlength)
        freq=np.linspace(1*1e9,9e9, 901)
        eps_r=1.2
        loss_t=0.15*0.001
        #d=0.56*1e-3
        
        em=em_env.EmEnv(freq,eps_rel=eps_r, loss_tan=loss_t) #set em env
        
        #set up geometry
        geo=geometry.Geometry()
        geo.set_length(cable_length)
        # define wires
        wir1=geometry.Wire((0,-d))#, radius=0.125*1e-3,sigma=1e7)
        wir2=geometry.Wire((-d/np.sqrt(2),-d/np.sqrt(2)))#, radius=0.125*1e-3,sigma=1e7)
        wir3=geometry.Wire((d/np.sqrt(2),-d/np.sqrt(2)))#, radius=0.125*1e-3,sigma=1e7)
        wir4=geometry.Wire((0,0))
        wir5=geometry.Wire((0,0.93*1e-3))
        pair=geometry.TwistedPair(wir4,wir5,twistlength, cable_length)
        # add wires
        geo.add_wire(wir1)
        geo.add_wire(wir2)
        geo.add_wire(wir3)
        # geo.add_wire(wir4)
        # geo.add_wire(wir5)
        geo.add_pair(pair)
        geo.plot_geometry(z_max=0.065)
    
    
        
        # set up simulation
        Nsamples=int(cable_length*400)
        sim=simulation.Simulation(geo, em, Nsamples)
        abcd=sim.run()
        
        
        s=sparameter.Sparameter()
        s.create_from_abcd(freq, abcd)
        paras=s.get_s()
        dum=1/np.sqrt(5)
        M=np.array([[dum,dum,dum,dum,dum,0,0,0,0,0],[0,0,0,0,0,dum,dum,dum,dum,dum]])
        s_sw=np.empty((np.size(freq),2,2), dtype=np.complex128)
        for idx in range(np.size(freq)):
            s_parameters=np.reshape(paras[idx,:,:],(10,10))
            s_sw[idx,:,:]=np.dot(np.dot(M,s_parameters),M.T)
        
        #s.plot()
        #z=np.linspace(1400,1800)
        z=[1500]
        for imp in z:
            ntw = rf.Network(frequency=freq, s=s_sw)
            ntw.renormalize(imp)
            s_renorm=sparameter.Sparameter()
            s_renorm.set_freq(freq)
            s_renorm.set_s(ntw.s)
            s_renorm.plot(ports=(2,1))
            s_renorm.plot(ports=(1,1))
        # name=f"Pitch_{twistlength}_d_{d:.5f}"
        # sio.savemat(name+".mat", {'f': freq, 'abcd': abcd})
        # ntw.write_touchstone(name+".s2p", )
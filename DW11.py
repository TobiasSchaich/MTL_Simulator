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

if __name__=='__main__':
    
    cable_length=6
    twistlength=1e-1
    freq=np.linspace(1e9,6e9, 901)
    eps_r=1.17
    loss_t=0.14*0.015
    
    em=em_env.EmEnv(freq,eps_rel=eps_r, loss_tan=loss_t) #set em env
    
    #set up geometry
    geo=geometry.Geometry()
    geo.set_length(cable_length)
    # define wires
    wir1=geometry.Wire((0,-1.4*1e-3), radius=0.125*1e-3,sigma=1e7)
    wir2=geometry.Wire((-0.7*1e-3,-0.7*1e-3), radius=0.125*1e-3,sigma=1e7)
    wir3=geometry.Wire((+0.7*1e-3,-0.7*1e-3), radius=0.125*1e-3,sigma=1e7)
    wir4=geometry.Wire((0,0))
    wir5=geometry.Wire((0,0.93*1e-3))
    pair=geometry.TwistedPair(wir4,wir5,twistlength, cable_length)
    # add wires
    geo.add_wire(wir1)
    geo.add_wire(wir2)
    geo.add_wire(wir3)
    #geo.add_wire(wir4)
    #geo.add_wire(wir5)
    geo.add_pair(pair)
    #geo.plot_geometry()


    
    # set up simulation
    Nsamples=cable_length*200
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
    z=np.linspace(1400,1800)
    z=[1500]
    for imp in z:
        ntw = rf.Network(frequency=freq, s=s_sw)
        ntw.renormalize(imp)
        s_renorm=sparameter.Sparameter()
        s_renorm.set_freq(freq)
        s_renorm.set_s(ntw.s)
        s_renorm.plot(ports=(2,1))
        s_renorm.plot(ports=(1,1))

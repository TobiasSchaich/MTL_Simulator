# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 11:06:02 2020

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
    d=1.*1e-3
    cable_length=0.05
    twistlength=0.05
    # random.seed(0)
    # twistlength=[]
    # while sum(twistlength) < cable_length:
    #     twistlength.append(0.065+(0.065*0.3*random.uniform(-1,1))) 
    # twistlength=tuple(twistlength)
    freq=np.linspace(1*1e9, 10e9, 901)
    eps_r=1
    loss_t=0    
    em=em_env.EmEnv(freq,eps_rel=eps_r, loss_tan=loss_t) #set em env
    
    #set up geometry
    geo=geometry.Geometry()
    geo.set_length(cable_length)
    # define wires
    wir1=geometry.Wire((0,d/2))
    wir2=geometry.Wire((0,-d/2))
    wir3=geometry.Wire((0,d))
    pair=geometry.TwistedPair(wir1,wir2,twistlength, cable_length)
    # add wires
    geo.add_wire(wir3)
    geo.add_pair(pair)
    geo.plot_geometry(z_max=0.05)
    
    # set up simulation
    Nsamples=int(cable_length*2000)
    sim=simulation.Simulation(geo, em, Nsamples)
    abcd=sim.run()
    abcd_orig=abcd
    #sio.savemat('Wire.mat', {'f': freq, 'abcd': abcd_orig})
    abcd=np.linalg.matrix_power(abcd_orig, 20)
    
    
    s=sparameter.Sparameter()
    s.create_from_abcd(freq, abcd)
    paras=s.get_s()
    dum=1/np.sqrt(3)
    M=np.array([[dum,dum,dum,0,0,0],[0,0,0,dum,dum,dum]])
    s_sw=np.empty((np.size(freq),2,2), dtype=np.complex128)
    for idx in range(np.size(freq)):
        s_parameters=np.reshape(paras[idx,:,:],(6,6))
        s_sw[idx,:,:]=np.dot(np.dot(M,s_parameters),M.T)
    
    #s.plot()
    #z=np.linspace(1400,1800)
    z=[1200]
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
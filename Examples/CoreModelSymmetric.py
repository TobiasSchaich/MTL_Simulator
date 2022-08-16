# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 11:06:02 2020

@author: Tobias

Simulates the core model with two straight wires symmetrically distributed around TWP (Scenario 2)
of IEEE Transactions on Microwave Theory and Techniques, vol. 70, no. 5, pp. 2541-2552 using MTL model. 
"""



import geometry
import em_env 
import simulation
import numpy as np
import sparameter
import skrf as rf
import matplotlib.pyplot as plt
#import random
#import scipy.io as sio
from scipy.constants import c

if __name__=='__main__':    
    d=1.5*1e-3
    twistlength=0.026
    cable_length=1
    freq=np.linspace(1e9, 12e9, 1601)
    eps_r=1
    loss_t=0    
    em=em_env.EmEnv(freq,eps_rel=eps_r, loss_tan=loss_t) #set em env
    
    #set up geometry
    geo=geometry.Geometry()
    geo.set_length(twistlength)
    # define wires
    wir1=geometry.CoatedWire(position=(0,0.5*1e-3), radius_coating=0.5*1e-3, ed_coating = 2.7, loss_tan = 0.01)
    wir2=geometry.CoatedWire(position=(0,-0.5*1e-3), radius_coating=0.5*1e-3, ed_coating = 2.7, loss_tan = 0.01)
    wir3=geometry.CoatedWire(position=(-d,0), radius_coating=0.5*1e-3, ed_coating = 2.7, loss_tan = 0.01)
    wir4=geometry.CoatedWire(position=(d,0), radius_coating=0.5*1e-3, ed_coating = 2.7, loss_tan = 0.01) #symmetric
    pair=geometry.TwistedPair(wir1,wir2,twistlength, twistlength)
    # add wires
    geo.add_wire(wir3)
    geo.add_wire(wir4)
    geo.add_pair(pair)
    geo.plot_geometry(z_max=twistlength)
    
    # set up simulation
    lambda_min=c/max(freq)
    Nsamples=int(twistlength*2000)

    sim=simulation.Simulation(geo, em, Nsamples)
    abcd=sim.run()
    
    #post-processing
    s=sparameter.Sparameter()
    s.create_from_abcd(freq, abcd)
    paras=s.get_s()
    ntw_50_short = rf.Network(frequency=freq, s=paras)
    exponent=int(np.ceil(cable_length/twistlength))
    ntw_50_long=ntw_50_short
    for _ in range(int(exponent)):
        ntw_50_long=ntw_50_long**ntw_50_short
    paras=ntw_50_long.s
    SW=1/np.sqrt(4)
    DM=1/np.sqrt(2)
    M=np.array([[SW,SW,SW,SW,0,0,0,0],[0,0,0,0, SW, SW,SW,SW], [DM,-DM,0,0,0,0,0,0], [0,0,0,0,DM,-DM,0,0], [0,0,DM,-DM,0,0,0,0], 
                [0,0,0,0,0,0, DM,-DM],[SW,SW,-SW,-SW,0,0,0,0], [0,0,0,0,SW,SW,-SW,-SW]])
    s_sw=np.empty((np.size(freq),8,8), dtype=np.complex128)
    for idx in range(np.size(freq)):
        s_parameters=np.reshape(paras[idx,:,:],(8,8))
        s_sw[idx,:,:]=np.dot(np.dot(M,s_parameters),M.T)
    ntw = rf.Network(frequency=freq/1e9, s=s_sw)
    # ntw.plot_s_db()
    # plt.show()
    Z=np.tile([4*216,4*216, 5, 5, 5, 5, 5, 5], np.size(freq)) # All shorted except SW
    Z=np.reshape(Z, [np.size(freq), 8])
    ntw.renormalize(Z)
    ntw.plot_s_db(n=0, m=1)
    filename="FourWireCoreModelMTL_Symmetric"+".s8p"
    ntw.write_touchstone(filename)
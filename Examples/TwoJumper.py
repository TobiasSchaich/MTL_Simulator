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
from math import gcd
from scipy.constants import c

if __name__=='__main__':    
    freq=np.linspace(1e9, 16e9, 1601)
    eps_r=1
    loss_t=0    
    em=em_env.EmEnv(freq,eps_rel=eps_r, loss_tan=loss_t) #set em env    
    cable_length=1
    twistlength1=0.012 #must be only to mm accuracy
    twistlength2=0.024 #must be only to mm accuracy
    if twistlength1==twistlength2:
        sim_length=twistlength1
    else:
        div=gcd(int(twistlength1*1e3),int(twistlength2*1e3))
        sim_length=twistlength1*(twistlength2*1e3/div) #length of shortest non-repeating cable length
        if sim_length>=cable_length:
            sim_length=cable_length

    #set up geometry
    geo=geometry.Geometry()
    geo.set_length(sim_length)
    # define wires
    wir1=geometry.CoatedWire(position=(1*1e-3,0.5*1e-3), radius_coating=0.5*1e-3, ed_coating = 2.7, loss_tan = 0.01)
    wir2=geometry.CoatedWire(position=(1*1e-3,-0.5*1e-3), radius_coating=0.5*1e-3, ed_coating = 2.7, loss_tan = 0.01)
    wir3=geometry.CoatedWire(position=(-1*1e-3,0.51*1e-3), radius_coating=0.5*1e-3, ed_coating = 2.7, loss_tan = 0.01)
    wir4=geometry.CoatedWire(position=(-1*1e-3,-0.51*1e-3), radius_coating=0.5*1e-3, ed_coating = 2.7, loss_tan = 0.01)
    pair1=geometry.TwistedPair(wir1,wir2,twistlength1, sim_length)
    pair2=geometry.TwistedPair(wir3,wir4,twistlength2, sim_length)
    # add wires
    geo.add_pair(pair1)
    geo.add_pair(pair2)
    geo.plot_geometry(z_max=twistlength1)
    
    # set up simulation
    Nsamples=int(sim_length*500) #every 2 mm

    sim=simulation.Simulation(geo, em, Nsamples)
    abcd=sim.run()
    
    #post-processing
    s=sparameter.Sparameter()
    s.create_from_abcd(freq, abcd)
    paras=s.get_s()
    ntw_50_short = rf.Network(frequency=freq, s=paras)
    exponent=int(np.ceil(cable_length/sim_length))
    ntw_50_long=ntw_50_short
    for _ in range(int(exponent)):
        ntw_50_long=ntw_50_long**ntw_50_short #combine S-parameters of each segment
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
    Z=np.tile([4*216,4*216, 5, 5, 5, 5, 5, 5], np.size(freq))
    Z=np.reshape(Z, [np.size(freq), 8])
    ntw.renormalize(Z)
    ntw.plot_s_db(n=0, m=1) #plot S21 for SW
    filename="TwoJumper"+".s8p"
    ntw.write_touchstone(filename)
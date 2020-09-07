# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 16:12:46 2020

@author: Tobias
"""

import geometry
import em_env 
import simulation
import numpy as np
import sparameter
import skrf as rf


if __name__=='__main__':
    
    cable_length=1
    # defines wires
    wir1=geometry.CoatedWire((-0.54*1e-3,0))
    # wir1=geometry.Wire((-0.45*1e-3,0))
    wir2=geometry.CoatedWire((0.54*1e-3,0))
    # wir3=geometry.Wire((-0.45*1e-3,2*1e-3))
    #wir4=geometry.Wire((0.45*1e-3,2*1e-3))
    
    #pair up wires
    pair=geometry.TwistedPair(wir1,wir2,1*1e-1, cable_length)
    #pair2=geometry.TwistedPair(wir3,wir4,0.8*1e-1,cable_length)
    
    #define geometry
    geo=geometry.Geometry()
    geo.set_length(cable_length)
    
    #add wires
    #geo.add_pair(pair)
    #geo.add_pair(pair2)
    geo.add_wire(wir1)
    geo.add_wire(wir2)
    #geo.add_wire(wir3)
    #set em env
    
    freq=np.logspace(7,9,200)
    eps_r=1#2.25
    em=em_env.EmEnv(freq,eps_rel=eps_r, loss_tan=.001)

    # set up simulation
    Nsamples=cable_length*100
    sim=simulation.Simulation(geo, em, Nsamples)
    abcd=sim.run()
    #s = rf.a2s(abcd)
    
    s=sparameter.Sparameter()
    s.create_from_abcd(freq, abcd)
    s.plot()
    ntw = rf.Network(frequency=freq, s=s.get_s())
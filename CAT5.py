# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 12:59:52 2020

@author: Tobias
"""

import geometry
import em_env 
import simulation
import numpy as np
import abcd2s

if __name__=='__main__':
    #input
    cable_length=2 #metres
    #freq=np.logspace(6,8.7, 100)
    freq=np.linspace(1e7,1e8,2)
    eps_r=2.25
    loss_tan=0.002
    # defines wires
    wir1=geometry.Wire((-0.5*1e-3,0), radius=0.255*1e-3)
    wir2=geometry.Wire((0.5*1e-3,0), radius=0.255*1e-3)
    wir3=geometry.Wire((-1.5*1e-3,1e-3), radius=0.255*1e-3)
    wir4=geometry.Wire((-0.5*1e-3,1*1e-3), radius=0.255*1e-3)
    wir5=geometry.Wire((0.5*1e-3,1*1e-3), radius=0.255*1e-3)
    wir6=geometry.Wire((1.5*1e-3,1*1e-3), radius=0.255*1e-3)
    wir7=geometry.Wire((-0.5*1e-3,2*1e-3), radius=0.255*1e-3)
    wir8=geometry.Wire((0.5*1e-3,2*1e-3), radius=0.255*1e-3)
    
    #pair up wires
    pair=geometry.TwistedPair(wir1,wir2,1.38*1e-2, cable_length)
    pair2=geometry.TwistedPair(wir3,wir4,1.53*1e-2,cable_length)
    pair3=geometry.TwistedPair(wir5,wir6, 1.78*1e-2,cable_length)
    pair4=geometry.TwistedPair(wir7,wir8, 1.94*1e-2,cable_length)

    #define geometry
    geo=geometry.Geometry()
    geo.set_length(cable_length)
    
    #add wires
    geo.add_pair(pair)
    geo.add_pair(pair2)
    geo.add_pair(pair3)
    geo.add_pair(pair4)
    
    #set em env

    em=em_env.EmEnv(freq,eps_rel=eps_r, loss_tan=loss_tan)

    # set up simulation
    sim=simulation.Simulation(geo, em, 400*cable_length)
    abcd=sim.run()
    s=abcd2s.abcd2s(abcd)
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 08:12:30 2020

@author: Tobias
"""

import math
import numpy as np
import matplotlib.pyplot as plt

class Wire: 
    '''Wire Class describing straight Wire in Coordinate system (2D)'''
    
    def __init__(self, position: tuple = (0,0),radius: float = 0.25*1e-3, sigma: float = 5.8*1e7):
        self._position=position
        self._radius=radius
        self._sigma=sigma
        self._wiretype='uncoated' #read only
        
    def get_radius(self):
        return self._radius

    def get_position(self):
        return self._position
    
    def get_sigma(self):
        return self._sigma
    
    def get_wiretype(self):
        return self._wiretype
     
class CoatedWire(Wire):
    
    def __init__(self, position = (0,0),radius = 0.25*1e-3, sigma = 5.8*1e7, radius_coating = 0.54*1e-3, 
                 ed_coating = 3, loss_tan = 0):
        super().__init__(position,radius, sigma)
        self._wiretype='coated'
        self._radius_coating=radius_coating
        self._ed_coating=ed_coating
        self._loss_tan=loss_tan
        
    def get_radius_coating(self):
        return self._radius_coating
    
    def get_ed_coating(self): 
        return self._ed_coating
    
    def get_loss_tan(self):
        return self._loss_tan
         
class TwistedPair:
    """ Creates a Twisted Pair in 3-Dimensional Space"""
    
    def __init__(self, wire1, wire2, pitch: tuple, length:float = 1.):
        if not all((isinstance(wire1, Wire), isinstance(wire2, Wire))):
            raise TypeError("wire1 and wire2 objects must be of type Wire")
        #Pitch can be Tuple with different pitch lengths if alternating along cable
        if isinstance(pitch, float) or isinstance(pitch, int):  
            self._pitch=(pitch,)
        elif isinstance(pitch, tuple):
            self._pitch=pitch
        self._wire1=wire1
        self._wire2=wire2
        self._length=length
        self._centerpoint=self.calculate_centerpoint()
        self._phi=self.calculate_phi()
        self._separation=self.calculate_separation()
                
    def calculate_separation(self):
        x1,y1=self._wire1.get_position()
        x2,y2=self._wire2.get_position()
        return math.sqrt((x1-x2)**2+(y1-y2)**2)
        
    def calculate_phi(self):
        x1,y1=self._wire1.get_position()
        x2,y2=self._wire2.get_position()
        xcenter,ycenter=self._centerpoint
        phi=(math.atan2(y1-ycenter,x1-xcenter), math.atan2(y2-ycenter,x2-xcenter))
        return phi
        
    def calculate_centerpoint(self):
        (x1,y1)=self._wire1.get_position()
        (x2,y2)=self._wire2.get_position()
        centerpoint=((x1+x2)/2, (y1+y2)/2)
        return centerpoint
    
    def get_length(self):
        return self._length
    
    def get_wires(self):
        return (self._wire1, self._wire2)
    
    def set_length(self, length):
        if isinstance(length, float) or isinstance(length, int):
            if length>0:
                self._length=float(length)
            else: 
                raise ValueError("Length has to be positive integer or float greater than 0.")
        else:
            raise TypeError("Length has to be of type int or float.")
    
    def get_xy(self, z):
        if z>self._length or z<0:
            raise ValueError("z must be between 0 and length of pair")
            
        z_working=z
        idx=0
        while True:
            if (z_working-self._pitch[idx%len(self._pitch)])<=0:
                break
            z_working=z_working-self._pitch[idx%len(self._pitch)]
            idx+=1
            idx=idx%len(self._pitch)
        sca_pitch=self._pitch[idx]
        x=self._get_x(sca_pitch,z_working)
        y=self._get_y(sca_pitch, z_working)
        return [x,y]   
        
    def _get_x(self, sca_pitch: float, z: float):
        x_center = self._centerpoint[0]
        x1 = self._separation/2*math.cos(2*math.pi/sca_pitch*z+self._phi[0])
        x2 = self._separation/2*math.cos(2*math.pi/sca_pitch*z+self._phi[1])
        return [x1+x_center, x2+x_center]
    
    def _get_y(self, sca_pitch: float, z: float):
        y_center=self._centerpoint[1]
        y1=self._separation/2*math.sin(2*math.pi/sca_pitch*z+self._phi[0])
        y2=self._separation/2*math.sin(2*math.pi/sca_pitch*z+self._phi[1])
        return [y1+y_center, y2+y_center]
    
    def plotpair(self, elevation=45, angle=45, Nsamples=1000):
        z=np.linspace(0,self._length, Nsamples)
        x1=np.zeros(Nsamples)
        x2=np.zeros(Nsamples)
        y1=np.zeros(Nsamples)
        y2=np.zeros(Nsamples)
        for i in range(Nsamples):
            x1[i],x2[i],y1[i],y2[i]=np.ravel(np.array(self.get_xy(z[i])))
        fig=plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(x1,y1,z)
        ax.plot(x2,y2,z)
        ax.view_init(elevation, angle)

class Geometry:
    """ Container for multiple wires, twisted pairs which need to be simulated"""
    
    def __init__(self, length=1.):
        self.set_length(length)
        self._num_pairs=0
        self._num_single_wires=0
        self._num_wires_total=0
        self._pairs=[]
        self._wires=[]
        
    def set_length(self, length):
        if isinstance(length, float) or isinstance(length, int):
            if length>0:
                self._length=float(length)
            else:
                raise ValueError("Length must have positive value.")
        else:
            raise TypeError("Length of geometry can only be float or integer.")
            
    def get_length(self): 
        return self._length
    
    def get_num_wires_total(self):
        return self._num_wires_total
        
    def get_num_single_wires(self):
        return self._num_single_wires
    
    def get_num_pairs(self):
        return self._num_pairs
    
    def get_all_wires(self):
        '''Returns List of all defined wires in geometry'''
        wires=[]
        for idx in range(self.get_num_single_wires()):
            wires.append(self.get_wire(idx))
        for idx in range(self.get_num_pairs()):
            wires.append(self.get_pair(idx).get_wires()[0])
            wires.append(self.get_pair(idx).get_wires()[1])
        return wires
    
    def get_all_positions(self, z):
        """Returns list of numpy arrays of shape (2,) carrying x and y position of wires at z.Returns first single wires then pairs"""
        pos_wires=[np.array(self.get_wire(idx).get_position()) for idx in range(self.get_num_single_wires())]
        pos_pairs=[np.transpose(self.get_pair(idx).get_xy(z)) for idx in range(self.get_num_pairs())]
        pos_pairs=[np.hsplit(np.ravel(pos_pairs[idx]),2) for idx in range(self.get_num_pairs())]
        pos_pairs=[val for wir in pos_pairs for val in wir] # flatten list
        pos=pos_wires+pos_pairs
        return pos
            
    
    def get_radii(self):
        """Returns two Arrays containing radii of wires and pairs which are empty if objects are defined"""
        rad1=[]
        rad2=[]
        for idx in range(self.get_num_single_wires()):
            rad1.append(self.get_wire(idx).get_radius())
        for idx in range(self.get_num_pairs()):
            wires=self.get_pair(idx).get_wires()
            rad2.append(wires[0].get_radius())
            rad2.append(wires[1].get_radius())
        return [rad1,rad2]
    
    def get_sigmas(self):
        """Returns two Arrays containing radii of wires and pairs which are empty if objects are defined"""
        sig1=[]
        sig2=[]
        for idx in range(self.get_num_single_wires()):
            sig1.append(self.get_wire(idx).get_sigma())
        for idx in range(self.get_num_pairs()):
            wires=self.get_pair(idx).get_wires()
            sig2.append(wires[0].get_sigma())
            sig2.append(wires[1].get_sigma())
        return [sig1,sig2]
    
    def add_wire(self, wire1: Wire):
        if isinstance(wire1, Wire):
            self._wires.append(wire1)
            self._num_single_wires+=1
            self._num_wires_total+=1
        else: 
            raise TypeError("addWire method only takes wire object as input.")
            
    def get_wire(self, idx):
        return self._wires[idx]
    
    def get_pair(self, idx): 
        return self._pairs[idx]
            
    def add_pair(self, pair: TwistedPair):
        if isinstance(pair, TwistedPair):
            if self._chk_length(pair):
                self._pairs.append(pair)
                self._num_wires_total+=2 #2 wires in pair
                self._num_pairs+=1
            else:
                raise ValueError("Length of Pair not compatible with length of geometry object.")
        else: 
            raise TypeError("add_pair method only takes TwistedPair object as input.")
    
    def _chk_length(self, pair: TwistedPair): 
        """Checks if pair length fits with geometry length"""
        pairlength=pair.get_length()
        if pairlength==self._length:
            return True
        else: 
            return False
        
    def plot_geometry(self, z_max=1,zsamples=1000, elev=45, angle=45):
        if z_max>self.get_length():
            raise ValueError("maximum z value zmax larger than length of cable")
        z=np.linspace(0,z_max, zsamples)
        ax = plt.axes(projection='3d')
        x0=np.zeros(zsamples)
        x1=np.zeros(zsamples)
        y0=np.zeros(zsamples)
        y1=np.zeros(zsamples)
        for pair in self._pairs:
            for i in range(zsamples):
                x0[i],x1[i], y0[i],y1[i]=np.ravel(np.array(pair.get_xy(z[i])))
            ax.plot(np.copy(x0),np.copy(y0),z)
            ax.plot(np.copy(x1),np.copy(y1),z)
        for wir in self._wires:
            wir_pos=np.array(wir.get_position())
            x=np.repeat(wir_pos[0],zsamples)
            y=np.repeat(wir_pos[0],zsamples)
            ax.plot(x, y, z)
        ax.view_init(elev, angle)    
        plt.show()
            
if __name__=='__main__':
    wir1=CoatedWire((1,2))
    wir2=Wire((2,-3))
    wir3=Wire((5,6))
    wir4=Wire((7,8))
    pair=TwistedPair(wir1,wir2,0.5, 5)
    pair2=TwistedPair(wir3,wir4,1,5)
    geo=Geometry()
    geo.set_length(5)
    geo.add_pair(pair)
    geo.add_pair(pair2)
    geo.plot_geometry(z_max=2, elev=45, angle=135)
    print(type(geo.get_all_wires()[0]))
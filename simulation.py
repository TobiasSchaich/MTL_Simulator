# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 13:41:02 2020

@author: Tobias
clear"""
import geometry
import em_env
import numpy as np
from scipy.special import hankel1, yn, jv, kn
from scipy.constants import mu_0, epsilon_0, c, pi
from scipy.linalg import expm



class Simulation:
    """Simulation environment combining EM Environment and Geometry"""
    def __init__(self, geo, em, zsamples):
        self.set_geo(geo)
        self.set_em(em)
        self.set_zsamples(zsamples)
        
        
        
        
        
    def set_geo(self, geo):
        if isinstance(geo, geometry.Geometry):
            self._geo=geo
        else:
            raise TypeError("geo must be of Type geometry.Geometry")
            
    def get_geo(self):
        return self._geo
    
    def set_em(self, em):
        if isinstance(em, em_env.EmEnv):
            self._em=em
        else: 
            raise TypeError("em_env must be of type em_env.EmEnv")
            
    def get_em(self):
        return self._em
    
    def set_zsamples(self, zsamples):
        if isinstance(zsamples, int) and zsamples>0:
            self._zsamples=zsamples
        else:
            raise TypeError("zsamples must be integer")
            
    def get_zsamples(self):
        return self._zsamples
    
    def goubau(self): 
        '''calculates propagation constant for SW on coated wire'''
        ##################Fixed variables######################################
        tol=1e-10       # tolerance for itteration
        Nmax = 4e4      # max no. iterations
        r=.1            # relaxation parameter for itteration
        ##################get environment and cable characteristics###################################
        freq=np.array(self.get_em().get_freq())
        ed=self.get_em().get_eps_rel()
        w = 2*pi*freq;
        k0 = w/c;           
        radii=[wir.get_radius() for wir in self.get_geo().get_all_wires()]
        radii_coating=[wir.get_radius_coating() for wir in self.get_geo().get_all_wires()]
        eds_coating=[wir.get_ed_coating() for wir in self.get_geo().get_all_wires()]
        if all([radii[0]==val for val in radii]):
             radius=radii[0]
        else:
            raise NotImplementedError("Wires with unequal radii not yet implemented")
        if all([radii_coating[0]==val for val in radii_coating]):
            radius_coating=radii_coating[0]
        else:
            raise NotImplementedError("Wires with unequal coatings not yet implemented")
        if all([eds_coating[0]==val for val in eds_coating]):
            ed_coating=eds_coating[0]
        else:
            raise NotImplementedError("Wires with unequal coatings dielectrics not yet implemented")   
        if ed_coating<ed:
            raise ValueError("Environment must have lower dielectric constant than coating for a bound SW solution")
        ###################Itteration###########################################
        be0 = 0.999*k0*np.sqrt(ed_coating) #start value for iteration
        N = 1;
        be = be0;       # initialize iteration
        while True:    
            h = np.sqrt(k0**2*ed_coating - be**2);
            g = np.sqrt(be**2 - k0**2*ed);
            z0=yn(0,h*radius)*jv(0,h*radius_coating)-yn(0,h*radius_coating)*jv(0,h*radius)
            z1=yn(0,h*radius)*jv(1,h*radius_coating)-yn(1,h*radius_coating)*jv(0,h*radius)
            fun = -(ed/ed_coating)*kn(1,g*radius_coating)/kn(0,g*radius_coating)*(z0/z1)
            bnew = r*k0*np.sqrt((ed+ed_coating*(fun**2))/(1+(fun**2))) + (1-r) * be;
            if np.linalg.norm(be-bnew)<tol:
                break     # exit loop           
            N = N + 1;
            be = bnew;
            if N>=Nmax:
                raise ValueError(f"goubau did not converge after {Nmax} itterations")       
        be = bnew;
        gd = np.sqrt((k0**2)*ed_coating-be**2)
        g = np.sqrt(be**2-(k0**2)*ed);
        return np.array([be,gd,g])

    
    def sommer(self): 
        '''calculates propagation constant for SW on coated wire'''
        def J01(z):
            P0 = 1 - 9/2./(8*z)**2 
            Q0 = -1./(8*z)
            P1 = 1 + 15/2./(8*z)**2 
            Q1 = 3./(8*z)        
            return 1j*(P0 + 1j*Q0)/(P1 + 1j*Q1);
    
        freq=self.get_em().get_freq()
        ed=self.get_em().get_eps_rel()
        tol=1e-10# iteration tolerance
        radii=[wir.get_radius() for wir in self.get_geo().get_all_wires()]
        sigs=[wir.get_sigma() for wir in self.get_geo().get_all_wires()]
        if all([radii[0]==val for val in radii]):
             radius=radii[0]
        else:
            raise NotImplementedError("Wires with unequal radii not yet implemented")
        if all([sigs[0]==val for val in sigs]):
            sigma=sigs[0]
        else:
            raise NotImplementedError("Wires with unequal sigma not yet implemented")     
        w=2*pi*freq   
        k0=w/c
        ec =  1- 1j * sigma/w/epsilon_0;
        be=0.9*k0
        ga = np.sqrt(ed*k0**2 - be**2);   
        N=1 #counter
        while True: #iterative method to approximate gamma
            gc = np.sqrt(ga**2 + k0**2.*(ec-ed));       
            gnew = hankel1(1,ga*radius)/hankel1(0,ga*radius)*ed*gc/ec*J01(gc*radius);  
            if np.linalg.norm(ga-gnew)<tol:
                break          
            N+=1        
            ga = gnew;       
            if N>1e4:
                raise ValueError("sommer did not converge")    
        be=np.sqrt(ed*k0**2-ga**2)
        return np.array([be,ga])
    
    def _calculate_rclg_coated(self, z): 
        be,gd,g=self.goubau() #Note that g is real instead of complex for increased computational accuracy (ga=1j*g)
        num_cond=self.get_geo().get_num_wires_total()   
        num_freq=np.size(self.get_em().get_freq())
        #all radii are equal for now
        radius_coating=self.get_geo().get_all_wires()[0].get_radius_coating()
        radius=self.get_geo().get_all_wires()[0].get_radius()
        sigma=self.get_geo().get_all_wires()[0].get_sigma()        
        z0 = lambda r : yn(0,gd*radius)*jv(0,gd*r)-yn(0,gd*r)*jv(0,gd*radius) #define z0,z1 function
        z1 = lambda r : yn(0,gd*radius)*jv(1,gd*r)-yn(1,gd*r)*jv(0,gd*radius)
        # get wire positions
        pos_wires=[np.array(self.get_geo().get_wire(idx).get_position()) for idx in range(self.get_geo().get_num_single_wires())]
        pos_pairs=[np.transpose(self.get_geo().get_pair(idx).get_xy(z)) for idx in range(self.get_geo().get_num_pairs())]
        pos_pairs=[np.hsplit(np.ravel(pos_pairs[idx]),2) for idx in range(self.get_geo().get_num_pairs())] # get position as xy array
        pos_pairs=[val for wir in pos_pairs for val in wir] # flatten list
        pos=pos_wires+pos_pairs
        #calculate Inductance and Resistance Matrix
        L_arr=np.tile(np.eye(self.get_geo().get_num_wires_total()), (num_freq,1,1))#init, create 3D unity matrix
        L_self=-mu_0/(2*pi*radius)*(z0(radius_coating)/(gd*z1(radius)) - z1(radius_coating)*kn(0, g*radius_coating)/(g*z1(radius)*kn(1,g*radius_coating)))#values for each frequency
        L_arr=L_arr*np.reshape(L_self, (num_freq,1,1)) #L_arr is 3D array with axis=0 - frequency samples, and axis2,3 being wire samples        
        for ind1 in range(num_cond):
             for ind2 in range((ind1+1),num_cond):
                 L_arr[:,ind1,ind2]=mu_0/(2*pi*g*radius)*(z1(radius_coating)*kn(0,g*np.linalg.norm(pos[ind1]-pos[ind2]))/(z1(radius)*kn(1,g*radius_coating)));
                 L_arr[:,ind2,ind1]=L_arr[:,ind1,ind2]
        delta=np.sqrt(1/(pi*self.get_em().get_freq()*mu_0*sigma));
        R_arr=np.tile(np.eye(num_cond),(num_freq,1,1))*np.reshape(1/(2*pi*radius*delta*sigma), (num_freq,1,1));
        #calculate Capacitance and conductance Matrix 
        ed_air=self.get_em().get_eps_rel()
        loss_tan_coating=self.get_geo().get_all_wires()[0].get_loss_tan()
        if loss_tan_coating!=0:
            ed_coating=self.get_geo().get_all_wires()[0].get_ed_coating()*(1-1j*loss_tan_coating)
        else:
            ed_coating=self.get_geo().get_all_wires()[0].get_ed_coating()
        P_arr=np.tile(np.eye(self.get_geo().get_num_wires_total()), (num_freq,1,1))#init, create 3D unity matrix
        P_self=-z1(radius_coating)/(2*pi*radius*epsilon_0*z1(radius)) * (z0(radius_coating)/
             (ed_coating*gd*z1(radius_coating))-kn(0,g*radius_coating)/(ed_air*g*kn(1, g*radius_coating))) #values for each frequency
        P_arr=P_arr*np.reshape(P_self, (num_freq,1,1)) #L_arr is 3D array with axis=0 - frequency samples, and axis2,3 being wire samples        
        for ind1 in range(num_cond):
             for ind2 in range((ind1+1),num_cond):
                 P_arr[:,ind1,ind2]=1/(2*pi*epsilon_0*g*radius)*(z1(radius_coating)*kn(0,g*np.linalg.norm(pos[ind1]-pos[ind2]))/(z1(radius)*kn(1,g*radius_coating)));
                 P_arr[:,ind2,ind1]=P_arr[:,ind1,ind2]
        C=np.linalg.inv(P_arr)
        C_arr=np.real(C)
        G_arr=-2*pi*np.reshape(self.get_em().get_freq(), (num_freq,1,1))*np.imag(C)
        return R_arr, C_arr, L_arr, G_arr                
        
    def _calculate_rclg_uncoated(self, z): 
        be,ga=self.sommer()
        num_cond=self.get_geo().get_num_wires_total()   
        num_freq=np.size(self.get_em().get_freq())
        #all radii are equal for now
        radius=self.get_geo().get_all_wires()[0].get_radius()
        sigma=self.get_geo().get_all_wires()[0].get_sigma()
        # get wire positions
        pos_wires=[np.array(self.get_geo().get_wire(idx).get_position()) for idx in range(self.get_geo().get_num_single_wires())]
        pos_pairs=[np.transpose(self.get_geo().get_pair(idx).get_xy(z)) for idx in range(self.get_geo().get_num_pairs())]
        pos_pairs=[np.hsplit(np.ravel(pos_pairs[idx]),2) for idx in range(self.get_geo().get_num_pairs())]
        pos_pairs=[val for wir in pos_pairs for val in wir] # flatten list
        pos=pos_wires+pos_pairs
     
        L_arr=np.tile(np.eye(self.get_geo().get_num_wires_total()), (num_freq,1,1))#init
        L_self=(mu_0/(2*pi*radius))*np.real(1j/ga**2*hankel1(0,ga*radius))/np.real(1j/ga*hankel1(1,ga*radius));#values for each frequency
        L_arr=L_arr*np.reshape(L_self, (num_freq,1,1))
        #L_arr is 3D array with axis=0 - frequency samples, and axis2,3 being wire samples   
     
        for ind1 in range(num_cond):
             for ind2 in range((ind1+1),num_cond):
                 L_arr[:,ind1,ind2]=(mu_0/(2*pi*radius))*(np.real(1j/ga**2*hankel1(0,ga*np.linalg.norm(pos[ind1]-pos[ind2])))/np.real(1j/ga*hankel1(1,ga*radius)));
                 L_arr[:,ind2,ind1]=L_arr[:,ind1,ind2]
        #P=1/(epsilon_0*self.get_em().get_eps_rel()*mu_0)*L_arr;
        #C=np.linalg.inv(P);
        L_inv=np.linalg.inv(L_arr)
        C_arr=(epsilon_0*self.get_em().get_eps_rel()*mu_0)*L_inv
        G_arr=np.reshape(mu_0*self.get_em().get_sigma_diel(), (num_freq,1,1))*L_inv; 
        delta=np.sqrt(1/(pi*self.get_em().get_freq()*mu_0*sigma));
        R_arr=np.tile(np.eye(num_cond),(num_freq,1,1))*np.reshape(1/(2*pi*radius*delta*sigma), (num_freq,1,1));
        return R_arr, C_arr, L_arr, G_arr
        
    def _get_mode(self):
        if all([wir.get_wiretype()=='coated' for wir in self.get_geo().get_all_wires()]):
            mode='coated'
        elif all([wir.get_wiretype()=='uncoated' for wir in self.get_geo().get_all_wires()]):
            mode='uncoated'
        else:
            mode='mixed'
            raise NotImplementedError("Mixed Geometries with uncoated and coated wires not yet implemented.")
        return mode
    
    def _chk_geo(self):
        if self.get_geo().get_num_wires_total()==0:
            raise ValueError("No wires defined in geometry.")
        else: 
            print("Geometry OK")
    
    def run(self):
        #check mode
        self._chk_geo()
        mode=self._get_mode() # coated,uncoated,mixed if geometry consists of respective wires types
        if mode=='uncoated':
            self._rclg_function=self._calculate_rclg_uncoated 
        elif mode=='coated':
            self._rclg_function=self._calculate_rclg_coated 
        elif mode=='mixed':
            raise NotImplementedError("Coated wires solution not yet implemented")
        else:
            raise ValueError("Mode should be of type uncoated, coated or mixed only")
        abcd=self._run_in_mode()    #pass correct function for given mode
        return abcd     
            
    def _run_in_mode(self):
        f=self.get_em().get_freq()
        num_freq=np.size(f)
        f=np.reshape(f,(num_freq,1,1))
        num_wires=self.get_geo().get_num_wires_total()

        abcd=np.tile(np.eye(2*num_wires,2*num_wires), (num_freq, 1, 1))
        abcd = abcd.astype(np.complex128, copy=False)
        z_values,dz=np.linspace(0, self.get_geo().get_length(), self.get_zsamples(), endpoint=False, retstep=True)
        for idx in range(self.get_zsamples()):
            z=z_values[idx]
            r_arr,c_arr,l_arr, g_arr=self._rclg_function(z)
            z_arr=r_arr+1j*2*pi*f*l_arr
            y_arr=g_arr+1j*2*pi*f*c_arr
            A=np.concatenate((np.concatenate((np.zeros((num_freq,num_wires, num_wires)), -z_arr), axis=2),
                              np.concatenate((-y_arr, np.zeros((num_freq,num_wires, num_wires))),axis=2)), axis=1)   
            for _ in range(num_freq):
                abcd[_,:,:]=np.dot(abcd[_,:,:],expm(-A[_,:,:]*dz))
        return abcd
    


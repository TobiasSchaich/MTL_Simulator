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


    
def J01(z):
    P0 = 1 - 9/2./(8*z)**2 
    Q0 = -1./(8*z)
    P1 = 1 + 15/2./(8*z)**2 
    Q1 = 3./(8*z)        
    return 1j*(P0 + 1j*Q0)/(P1 + 1j*Q1);


class Simulation:
    """Simulation environment combining EM Environment and Geometry"""
    def __init__(self, geo, em, zsamples):
        self.set_geo(geo)
        self.set_em(em)
        self.set_zsamples(zsamples)
        self._identical_wires_flag=False
        self._clg_function=None
        
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
    
    def run(self):
        #check mode
        self._chk_geo()
        mode=self._get_mode() # coated,uncoated,mixed if geometry consists of respective wires types
        if mode=='uncoated':
            self._clg_function=self._calculate_clg_uncoated 
        elif mode=='coated':
            self._clg_function=self._calculate_clg_coated 
        elif mode=='mixed':
            print("Warning: Experimental Mode - not yet checked for errors")
            raise NotImplementedError("Coated wires solution not yet implemented")
        else:
            raise ValueError("Mode should be of type uncoated, coated or mixed only")
        abcd=self._run_in_mode()    #pass correct function for given mode
        return abcd     
    
    def goubau(self): 
        '''calculates propagation constant for SW on coated wire'''
        ##################Fixed variables######################################
        tol=1e-10       # tolerance for itteration
        Nmax = 4e4      # max no. iterations
        r=.1            # relaxation parameter for itteration
        # variables via em_env class
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
        #Itteration 
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


    def _sommer(self): 
        '''calculates propagation constant for SW on uncoated wire'''
        if self._identical_wires_flag:
            radius=self.get_geo().get_all_wires()[0].get_radius()
            sigma=self.get_geo().get_all_wires()[0].get_sigma()
            be,ga=self._single_wire_sommer(sigma, radius)
            return np.array([be,ga])
        else: 
            num_wires=self.get_geo().get_num_wires_total()
            radii=[wir.get_radius() for wir in self.get_geo().get_all_wires()]
            sigs=[wir.get_sigma() for wir in self.get_geo().get_all_wires()]
            num_freq=np.size(self.get_em().get_freq())
            results=np.empty((2,num_freq,num_wires), dtype=np.complex128) # axis0=be,ga, axis1=freq points, axis2=wire
            for idx in range(num_wires):
                be,ga=self._single_wire_sommer(sigs[idx], radii[idx])
                results[0,:,idx]=be
                results[1,:,idx]=ga
            return results
    
    def _single_wire_sommer(self,sigma, radius):
        freq=self.get_em().get_freq()
        ed=self.get_em().get_eps_rel()
        tol=1e-10# iteration tolerance
        w=2*pi*freq   
        k0=w/c
        ec =  1- 1j * sigma/w/epsilon_0;
        be=0.9*k0*ed
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
    
    def _clg_identical_coated(self, z): 
        ed_air=self.get_em().get_eps_rel()
        ed_coating=self.get_geo().get_all_wires()[0].get_ed_coating_complex()
        be,gd,g=self.goubau() #Note that g is real instead of complex for increased computational accuracy (ga=1j*g)
        num_cond=self.get_geo().get_num_wires_total()   
        num_freq=np.size(self.get_em().get_freq())
        radius_coating=self.get_geo().get_all_wires()[0].get_radius_coating()
        radius=self.get_geo().get_all_wires()[0].get_radius()   
        pos=self.get_geo().get_all_positions(z)
        z0 = lambda r : yn(0,gd*radius)*jv(0,gd*r)-yn(0,gd*r)*jv(0,gd*radius) #define z0,z1 function
        z1 = lambda r : yn(0,gd*radius)*jv(1,gd*r)-yn(1,gd*r)*jv(0,gd*radius)
        #calculate Inductance Matrix
        L_arr=np.tile(np.eye(self.get_geo().get_num_wires_total()), (num_freq,1,1))#init, create 3D unity matrix
        L_self=-mu_0/(2*pi*radius)*z0(radius_coating)/(gd*z1(radius))*(1+
                                    gd**2*ed_air/(g**2*np.real(ed_coating))) #values for each frequency
        L_arr=L_arr*np.reshape(L_self, (num_freq,1,1)) #L_arr is 3D array with axis=0 - frequency samples, and axis2,3 being wire samples        
        for ind1 in range(num_cond):
              for ind2 in range((ind1+1),num_cond):
                  dist=np.linalg.norm(pos[ind1]-pos[ind2])
                  L_arr[:,ind1,ind2]=mu_0/(2*pi*g*radius)*(z1(radius_coating)*kn(0,g*dist)/(z1(radius)*kn(1,g*radius_coating)));
                  L_arr[:,ind2,ind1]=L_arr[:,ind1,ind2]
        #calculate Capacitance and conductance Matrix 
        P_arr=np.tile(np.eye(self.get_geo().get_num_wires_total()), (num_freq,1,1))#init, create 3D unity matrix
        P_self=-z1(radius_coating)/(2*pi*radius*epsilon_0*z1(radius)) * (z0(radius_coating)/
              (ed_coating*gd*z1(radius_coating))-kn(0,g*radius_coating)/(ed_air*g*kn(1, g*radius_coating))) #values for each frequency
        P_arr=P_arr*np.reshape(P_self, (num_freq,1,1)) #L_arr is 3D array with axis=0 - frequency samples, and axis2,3 being wire samples        
        for ind1 in range(num_cond):
              for ind2 in range((ind1+1),num_cond):
                  dist=np.linalg.norm(pos[ind1]-pos[ind2])
                  P_arr[:,ind1,ind2]=1/(2*pi*epsilon_0*g*radius)*(z1(radius_coating)*kn(0,g*dist)/(z1(radius)*kn(1,g*radius_coating)));
                  P_arr[:,ind2,ind1]=P_arr[:,ind1,ind2]
        C=np.linalg.inv(P_arr) # Extract complex C-matrix
        C[0,:,:]=np.linalg.inv(P_arr[0,:,:])#fix for strange error in matrix inversion on my desktop only
        C_arr=np.real(C) #capacitance matrix is real of complex C-Matrix
        G_arr=-2*pi*np.reshape(self.get_em().get_freq(), (num_freq,1,1))*np.imag(C)
        return C_arr, L_arr, G_arr                
    
    def _calculate_clg_coated(self, z): 
        if self._identical_wires_flag:
            return self._clg_identical_coated(z)
        else:
            raise NotImplementedError("Coated wires of different radii not yet implemented")
            return self._clg_unidentical_uncoated(z)
        
    def _calculate_clg_uncoated(self, z): 
        if self._identical_wires_flag:
            return self._clg_identical_uncoated(z)
        else:
            return self._clg_unidentical_uncoated(z)
        
    def _clg_identical_uncoated(self, z):
        """Calculates CLG matrices for setup with all wires equal"""
        be,ga=self._sommer()
        eps_rel=self.get_em().get_eps_rel()
        num_cond=self.get_geo().get_num_wires_total()   
        num_freq=np.size(self.get_em().get_freq())
        radius=self.get_geo().get_all_wires()[0].get_radius()
        sigma=self.get_geo().get_all_wires()[0].get_sigma()
        pos=self.get_geo().get_all_positions(z) #wire positions as list of (2,) numpy arrays carrying x and y information
        #calculate clg matrix
        L_arr=np.tile(np.eye(self.get_geo().get_num_wires_total()), (num_freq,1,1))#init
        L_self=(mu_0/(2*pi*ga*radius))*hankel1(0,ga*radius)/hankel1(1,ga*radius);#values for each frequency
        L_arr=L_arr*np.reshape(L_self, (num_freq,1,1))
        #L_arr is 3D array with axis=0 - frequency samples, and axis2,3 being wire samples   
        for ind1 in range(num_cond):
             for ind2 in range((ind1+1),num_cond):
                 dist=np.linalg.norm(pos[ind1]-pos[ind2])
                 L_arr[:,ind1,ind2]=(mu_0/(2*pi*ga*radius))*hankel1(0,ga*dist)/hankel1(1,ga*radius)
                 L_arr[:,ind2,ind1]=L_arr[:,ind1,ind2]
        P_arr=1/(epsilon_0*eps_rel*mu_0)*L_arr
        C_arr=np.linalg.inv(P_arr)
        C_arr[0,:,:]=np.linalg.inv(P_arr[0,:,:])#fix for strange error in matrix inversion on my desktop only
        G_arr=np.reshape(self.get_em().get_sigma_diel(), (num_freq,1,1))/(epsilon_0*eps_rel)*C_arr
        return C_arr, L_arr, G_arr
    
    def _clg_unidentical_uncoated(self, z):
        """Calculates CLG for setup with wires of different radius or conductivity"""
        propa=self._sommer() # propagation array axis0-be,ga, axis1-frequency, axis2-wire
        num_cond=self.get_geo().get_num_wires_total()   
        num_freq=np.size(self.get_em().get_freq())
        eps_rel=self.get_em().get_eps_rel()
        sigma_diel=self.get_em().get_sigma_diel()
        radii=[wir.get_radius() for wir in self.get_geo().get_all_wires()]
        pos=self.get_geo().get_all_positions(z) #wire positions as list of (2,) numpy arrays carrying x and y information
     
        L_arr=np.empty((num_freq, num_cond,num_cond), dtype=np.complex128)#init
        #L_arr is 3D array with axis=0 - frequency samples, and axis2,3 being wire samples   
        for ind1 in range(num_cond):
              for ind2 in range(num_cond):
                  #be=propa[0,:,ind2] # not needed 
                  ga=propa[1,:,ind2]
                  radius=radii[ind1]
                  if ind1!=ind2:
                      L_arr[:,ind1,ind2]=(mu_0/(2*pi*ga*radius))*hankel1(0,ga*np.linalg.norm(pos[ind1]-pos[ind2]))/hankel1(1,ga*radius);
                  else:
                      L_arr[:,ind1,ind2]=(mu_0/(2*pi*ga*radius))*hankel1(0,ga*radius)/hankel1(1,ga*radius);
        #P=1/(epsilon_0*self.get_em().get_eps_rel()*mu_0)*L_arr;
        #C=np.linalg.inv(P);
        L_inv=np.linalg.inv(L_arr)
        L_inv[0,:,:]=np.linalg.inv(L_arr[0,:,:])#fix for strange error in matrix inversion
        C_arr=(epsilon_0*eps_rel*mu_0)*L_inv
        G_arr=np.reshape(mu_0*sigma_diel, (num_freq,1,1))*L_inv; 
        return C_arr, L_arr, G_arr
    
    def _calculate_r_matrix(self):
        """Calculates resistance matrix for coated and uncoated wires"""
        sigs=[wir.get_sigma() for wir in self.get_geo().get_all_wires()]
        radii=[wir.get_radius() for wir in self.get_geo().get_all_wires()]
        freqs=self.get_em().get_freq()
        num_freq=np.size(freqs)
        num_cond=self.get_geo().get_num_wires_total()
        if self._identical_wires_flag:
            sigma=sigs[0] #identical sigma for all wires
            radius=radii[0]
            delta=np.sqrt(1/(pi*freqs*mu_0*sigma))
            identity_3D=np.tile(np.eye(num_cond),(num_freq,1,1))
            resistances=np.reshape(1/(2*pi*radius*delta*sigma), (num_freq,1,1))
            R_arr=identity_3D*resistances
        else:
            sigs_tiled=np.tile(sigs,(num_freq,1)) # creates 2D array axis0-freq, axis 1 - conductivities
            radii_tiled=np.tile(radii,(num_freq,1)) # creates 2D array axis0-freq, axis 1 - radii
            delta=sigs_tiled*np.reshape(freqs,(num_freq,1))
            delta=np.sqrt(1/(pi*mu_0*delta));
            resistances=1/(2*pi*radii_tiled*delta*sigs_tiled)
            R_arr=np.empty((num_freq,num_cond,num_cond))
            for ind in range(num_freq):
                R_arr[ind,:,:]=np.diag(resistances[ind,:])
        return R_arr
        
    def _get_mode(self):
        if all([wir.get_wiretype()=='coated' for wir in self.get_geo().get_all_wires()]):
            mode='coated'
            #raise NotImplementedError(" Geometries with coated wires not yet implemented.")
        elif all([wir.get_wiretype()=='uncoated' for wir in self.get_geo().get_all_wires()]):
            mode='uncoated'
        else:
            mode='mixed'
            raise NotImplementedError("Mixed Geometries with uncoated and coated wires not yet implemented.")
        return mode
    
    def _chk_identical_wires(self):
        mode=self._get_mode()
        if mode=='mixed': #mix of coated and uncoated wires
            self._identical_wires_flag=False
            return
        radii=[wir.get_radius() for wir in self.get_geo().get_all_wires()]
        sigs=[wir.get_sigma() for wir in self.get_geo().get_all_wires()]
        ident_radii_flag=all(r==radii[0] for r in radii)
        ident_sigma_flag=all(sig==sigs[0] for sig in sigs)
        if ident_radii_flag and ident_sigma_flag:
            self._identical_wires_flag=True
            print("All wires identical. Optimised solver will be used.")
        else:
            self._identical_wires_flag=False
    
    def _chk_geo(self):
        if self.get_geo().get_num_wires_total()==0:
            raise ValueError("No wires defined in geometry.")
        self._chk_identical_wires()            
        print("Geometry OK")
    

            
    def _run_in_mode(self):
        f=self.get_em().get_freq()
        num_freq=np.size(f)
        f=np.reshape(f,(num_freq,1,1))
        num_wires=self.get_geo().get_num_wires_total()
        
        abcd=np.tile(np.eye(2*num_wires,2*num_wires), (num_freq, 1, 1))
        abcd = abcd.astype(np.complex128, copy=False)
        z_values,dz=np.linspace(0, self.get_geo().get_length(), self.get_zsamples(), endpoint=False, retstep=True)
        r_arr=self._calculate_r_matrix() #independent of z
        for idx in range(self.get_zsamples()):
            z=z_values[idx]
            if idx % 10==0: #every 10 steps
                print(f"Calculating for z={z}")
            c_arr,l_arr, g_arr=self._clg_function(z)
            z_arr=r_arr+1j*2*pi*f*l_arr
            y_arr=g_arr+1j*2*pi*f*c_arr
            A=np.concatenate((np.concatenate((np.zeros((num_freq,num_wires, num_wires)), -z_arr), axis=2),
                              np.concatenate((-y_arr, np.zeros((num_freq,num_wires, num_wires))),axis=2)), axis=1)   
            for _ in range(num_freq):
                abcd[_,:,:]=np.dot(abcd[_,:,:],expm(-A[_,:,:]*dz))
        return abcd
    


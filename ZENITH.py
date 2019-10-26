import numpy as np
import scipy as sp
import scipy.linalg
import scipy.signal
import matplotlib as mpl
import matplotlib.pyplot as plt

from AeroDB import *

class ZENITH(object):
    def __init__(self, aerodb, x0=None, rho=None):
        self.aerodb = aerodb
        self.x0 = x0
        self.rho = rho
        self.x_l = [self.x0]
        self.t_l = [0]
        
    @property
    def aerodb(self):
        return self._aerodb

    @aerodb.setter
    def aerodb(self, aerodb):
        self._aerodb = aerodb
        
    @property
    def x(self):
        return self._x

    @x.setter
    def x(self, x):
        self._x = x
        
    @property
    def xE(self):
        R = 6378.15e3
        x1, x2, x3 = self.pos
        xE = np.array([
            [R*np.arctan(np.sqrt(x1**2 + x3**2)/(R+x2))*(x1/np.sqrt(x1**2 + x3**2))],
            [np.sqrt((R+x2)**2 + (x1**2 + x3**2) - R)],
            [R*np.arctan(np.sqrt(x1**2 + x3**2)/(R+x2))*(x3/np.sqrt(x1**2 + x3**2))]
        ])
        
        return xE
        
    @property
    def pos(self):
        return self.x[:3]

    @property
    def vel(self):
        return self.x[3:6]

    @property
    def ang(self):
        return self.x[6:9]

    @property
    def angvel(self):
        return self.x[9:]
        
    @property
    def x0(self):
        return self._x0

    @x0.setter
    def x0(self, x0):
        if x0 is None:
            x0 = np.zeros((12,1))
        else:
            pass
        
        try:
            x0.shape
        except:
            raise TypeError("Initial state x0 must be a numpy array.")
        else:
            if x0.shape[0] != 1:
                x0 = x0.reshape((x0.shape[0],1))
            else:
                pass

            if x0.shape[0] != 12:
                raise ValueError("Initial state x0 must a full state (12-row) vector.")
            else:
                self.x = x0
                self._x0 = x0
                
    def set_rho(self, rho):
        if rho is None:
            self._rho = lambda h : 1.225
        elif callable(rho):
            self._rho = rho
        elif isinstance(rho, float):
            self._rho = lambda h : rho
        else:
            raise TypeError("Invalid rho type")
            
    def set_I(self, I, t_I):
        Idot = np.gradient(I, t_I)

        self._Ifunct = lambda t : np.interp(t, t_I, I)
        self._Idotfunct = lambda t : np.interp(t, t_I, Idot)
        
    def get_I(self, t):
        return self._Ifunct(t)

    def get_Idot(self, t):
        return self._Idotfunct(t)

    def set_Ip(self, Ip, t_Ip):
        Ipdot = np.gradient(Ip, t_Ip)

        self._Ipfunct = lambda t : np.interp(t, t_Ip, Ip)
        self._Ipdotfunct = lambda t : np.interp(t, t_Ip, Ipdot)
        
    def get_Ip(self, t):
        return self._Ipfunct(t)

    def get_Ipdot(self, t):
        return self._Ipdotfunct(t)
            
    def yaw(self, Vw=None):
        if Vw is None:
            Vw = np.zeros((3,1))
        else:
            pass
        
        posDOTvelw = self.pos.flatten().dot((self.vel+Vw).flatten())
        PosMULVelW = np.linalg.norm(self.pos)*np.linalg.norm(self.vel+Vw)
        
        yaw = np.arccos(posDOTvelw/PosMULVelW)
        
        if np.isnan(yaw):
            yaw = 0
        
        return yaw
        
    def rho(self, h=0):
        return self._rho(h)
                
    def xdot(self, t, x):
        self.x = x
        self.x_l.append(x) 
        self.t_l.append(t)
        
        z = self.x[2]
        rho = self.rho(z)
        I = self.get_I(t)
        Ip = self.get_Ip(t)
        IDIVIp = I/Ip
        
        self.aerodb
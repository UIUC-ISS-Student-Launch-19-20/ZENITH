class AeroDB(object):
    def __init__(self, aerodb=None, cal=1, area=None):
        if isinstance(aerodb, AeroDB):
            self.clone(aerodb)
        else:
            self.clear()
        
        if cal is None and area is None:
            raise ValueError("Provide at least one of caliber and area as dimensional input.")
        elif area is None:
            self.cal = cal
            self.area = area
        else:
            self.area = area
            self.cal = cal
            
        
    def clear(self):
        self.cal = 1
        self.area = None
        self.KD0 = 0
        self.KDa = 0
        
    def C2K(self, C):
        return C*self.area/(2*self.cal**2)
    
    def K2C(self, K):
        return (2*K*self.cal**2)/self.area
        
    @property
    def area(self):
        return self._area
    
    @area.setter
    def area(self, area):
        if area is None:
            self._area = np.sqrt(np.pi)*self.cal/2
        elif np.isscalar(area):
            area = np.abs(area)
            self._area = area
        else:
            raise TypeError("Area must be a scalar.")
        
    @property
    def cal(self):
        return self._cal
    
    @cal.setter
    def cal(self, cal):
        if cal is None:
            self._cal = 2*self.area/np.sqrt(np.pi)
        elif np.isscalar(cal):
            cal = np.abs(cal)
            self._cal = cal
        else:
            raise TypeError("Caliber must be a scalar.")
        
    @property
    def KD0(self):
        return self._KD0
    
    @KD0.setter
    def KD0(self, KD0):
        if np.isscalar(KD0):
            KD0 = np.abs(KD0)
            self._KD0 = KD0
        else:
            raise TypeError("Zero-lift drag force coefficient KD0 (AB) must be a scalar.")
            
    @property
    def CD0(self):
        return self.K2C(self.KD0)
    
    @CD0.setter
    def CD0(self, CD0):
        if np.isscalar(CD0):
            CD0 = np.abs(CD0)
            self.KD0 = self.C2K(CD0)
        else:
            raise TypeError("Zero-lift drag force coefficient CD0 (AD) must be a scalar.")
    
    @property
    def KDa(self):
        return self._KDa
    
    @KDa.setter
    def KDa(self, KDa):
        if np.isscalar(KDa):
            KDa = np.abs(KDa)
            self._KDa = KDa
        else:
            raise TypeError("Angle-dep. drag force coefficient KDa (AB) must be a scalar.")
            
    @property
    def CDa(self):
        return self.K2C(self.KDa)
    
    @CDa.setter
    def CDa(self, CDa):
        if np.isscalar(CDa):
            CDa = np.abs(CDa)
            self.KDa = self.C2K(CDa)
        else:
            raise TypeError("Angle-dep. drag force coefficient CDa (AD) must be a scalar.")
            
    @property
    def KDa2(self):
        return self._KDa2
    
    @KDa.setter
    def KDa2(self, KDa2):
        if np.isscalar(KDa2):
            KDa2 = np.abs(KDa2)
            self._KDa2 = KDa2
        else:
            raise TypeError("Square angle-dep. drag force coefficient KDa2 (AB) must be a scalar.")
            
    @property
    def CDa2(self):
        return self.K2C(self.KDa2)
    
    @CDa2.setter
    def CDa2(self, CDa2):
        if np.isscalar(CDa2):
            CDa2 = np.abs(CDa2)
            self.KDa2 = self.C2K(CDa2)
        else:
            raise TypeError("Square angle-dep. drag force coefficient CDa2 (AD) must be a scalar.")
            
    @property
    def KA(self):
        return self._KA
    
    @KA.setter
    def KA(self, KA):
        if np.isscalar(KA):
            KA = np.abs(KA)
            self._KA = KA
        else:
            raise TypeError("Spin damping moment coefficient KA (AB) must be a scalar.")
            
    @property
    def CA(self):
        return self.K2C(self.KA)
    
    @CA.setter
    def CA(self, CA):
        if np.isscalar(CA):
            CA = np.abs(CA)
            self.KA = self.C2K(CA)
        else:
            raise TypeError("Spin damping moment coefficient CA (AD) must be a scalar.")
            
    @property
    def KE(self):
        return self._KE
    
    @KE.setter
    def KE(self, KE):
        if np.isscalar(KE):
            KE = np.abs(KE)
            self._KE = KE
        else:
            raise TypeError("Fin cant moment coefficient KE (AB) must be a scalar.")
            
    @property
    def CE(self):
        return self.K2C(self.KE)
    
    @CE.setter
    def CE(self, CE):
        if np.isscalar(CE):
            CE = np.abs(CE)
            self.KE = self.C2K(CE)
        else:
            raise TypeError("Fin cant moment coefficient CE (AD) must be a scalar.")
            
    @property
    def KL(self):
        return self._KL
    
    @KL.setter
    def KL(self, KL):
        if np.isscalar(KL):
            KL = np.abs(KL)
            self._KL = KL
        else:
            raise TypeError("Lift force coefficient KL (AB) must be a scalar.")
            
    @property
    def CL(self):
        return self.K2C(self.KL)
    
    @CL.setter
    def CL(self, CL):
        if np.isscalar(CL):
            CL = np.abs(CL)
            self.KL = self.C2K(CL)
        else:
            raise TypeError("Lift force coefficient CL (AD) must be a scalar.")
            
    @property
    def KM(self):
        return self._KM
    
    @KM.setter
    def KM(self, KM):
        if np.isscalar(KM):
            KM = np.abs(KM)
            self._KM = KM
        else:
            raise TypeError("Overturning moment coefficient KM (AB) must be a scalar.")
            
    @property
    def CM(self):
        return self.K2C(self.KM)
    
    @CM.setter
    def CM(self, CM):
        if np.isscalar(CM):
            CM = np.abs(CM)
            self.KM = self.C2K(CM)
        else:
            raise TypeError("Overturning moment coefficient CM (AD) must be a scalar.")
            
    @property
    def KF(self):
        return self._KF
    
    @KF.setter
    def KF(self, KF):
        if np.isscalar(KF):
            KF = np.abs(KF)
            self._KF = KF
        else:
            raise TypeError("Magnus force coefficient KF (AB) must be a scalar.")
            
    @property
    def CF(self):
        return self.K2C(self.KF)
    
    @CF.setter
    def CF(self, CF):
        if np.isscalar(CF):
            CF = np.abs(CF)
            self.KF = self.C2K(CF)
        else:
            raise TypeError("Magnus force coefficient CF (AD) must be a scalar.")
            
    @property
    def KT(self):
        return self._KT
    
    @KT.setter
    def KT(self, KT):
        if np.isscalar(KT):
            KT = np.abs(KT)
            self._KT = KT
        else:
            raise TypeError("Magnus moment coefficient KT (AB) must be a scalar.")
            
    @property
    def CT(self):
        return self.K2C(self.KT)
    
    @CT.setter
    def CT(self, CT):
        if np.isscalar(CT):
            CT = np.abs(CT)
            self.KT = self.C2K(CT)
        else:
            raise TypeError("Magnus moment coefficient CT (AD) must be a scalar.")
            
    @property
    def KS(self):
        return self._KS
    
    @KS.setter
    def KS(self, KS):
        if np.isscalar(KS):
            KS = np.abs(KS)
            self._KS = KS
        else:
            raise TypeError("Pitching force coefficient KS (AB) must be a scalar.")
            
    @property
    def CS(self):
        return self.K2C(self.KS)
    
    @CS.setter
    def CS(self, CS):
        if np.isscalar(CS):
            CS = np.abs(CS)
            self.KS = self.C2K(CS)
        else:
            raise TypeError("Pitching force coefficient CS (AD) must be a scalar.")
            
    @property
    def KH(self):
        return self._KH
    
    @KH.setter
    def KH(self, KH):
        if np.isscalar(KH):
            KH = KH
            self._KH = KH
        else:
            raise TypeError("Damping moment coefficient KH (AB) must be a scalar.")
            
    @property
    def CH(self):
        return self.K2C(self.KH)
    
    @CH.setter
    def CH(self, CH):
        if np.isscalar(CH):
            CH = CH
            self.KH = self.C2K(CH)
        else:
            raise TypeError("Damping moment coefficient CH (AD) must be a scalar.")
            
    @property
    def KXF(self):
        return self._KXF
    
    @KXF.setter
    def KXF(self, KXF):
        if np.isscalar(KXF):
            KXF = np.abs(KXF)
            self._KXF = KXF
        else:
            raise TypeError("Magnus cross force coefficient KXF (AB) must be a scalar.")
            
    @property
    def CXF(self):
        return self.K2C(self.KXF)
    
    @CXF.setter
    def CXF(self, CXF):
        if np.isscalar(CXF):
            CXF = np.abs(CXF)
            self.KXF = self.C2K(CXF)
        else:
            raise TypeError("Magnus cross force coefficient CXF (AD) must be a scalar.")
            
    @property
    def KXT(self):
        return self._KXT
    
    @KXT.setter
    def KXT(self, KXT):
        if np.isscalar(KXT):
            KXT = np.abs(KXT)
            self._KXT = KXT
        else:
            raise TypeError("Magnus cross moment coefficient KXT (AB) must be a scalar.")
            
    @property
    def CXT(self):
        return self.K2C(self.KXT)
    
    @CXT.setter
    def CXT(self, CXT):
        if np.isscalar(CXT):
            CXT = np.abs(CXT)
            self.KXT = self.C2K(CXT)
        else:
            raise TypeError("Magnus cross moment coefficient CXT (AD) must be a scalar.")
            
    def yaw(self, x, v):
        xDOTv = x.flatten().dot(v.flatten())
        XMULV = np.linalg.norm(x)*np.linalg.norm(v)
        
        yaw = np.arccos(xDOTv/XMULV)
        
        if np.isnan(yaw):
            yaw = 0
        
        return yaw
    
    def Rtilde(self, x, h, IDIVIp):
        hDOTx = h.flatten().dot(x.flatten())
        
        return IDIVIp*hDOTx
            
    def dragForce(self, x, v, rho=1.225):
        V = np.linalg.norm(v)
        yaw = self.yaw(x, v)
        
        return -rho*self.cal**2 * (self.KD0 + self.KDa*yaw + self.KDa2*yaw**2)*V*v
    
    def spindampingMoment(self, x, v, h, IDIVIp, rho=1.225):
        V = np.linalg.norm(v)
        Rtilde = self.Rtilde(x, h, IDIVIp)
        
        return -rho*self.cal**4 * self.KA * Rtilde * V * x
    
    def fincantMoment(self, x, v, eps=0, rho=1.225):
        V = np.linalg.norm(v)
        
        return -rho*self.cal**3 * self.KE * eps * V**2 * x
    
    def liftForce(self, x, v, rho=1.225):
        V = np.linalg.norm(v)
        vDOTx = v.flatten().dot(x.flatten())
        
        return rho*self.cal**2 * self.KL * (V**2 * x - vDOTx*v)
    
    def overturningMoment(self, x, v, rho=1.225):
        V = np.linalg.norm(v)
        vCRSx = np.cross(v.flatten(), x.flatten())
        
        return rho*self.cal**3 * self.KM * V * vDOTx
    
    def magnusForce(self, x, v, h, IDIVIp, rho=1.225):
        xCRSv = np.cross(x.flatten(), v.flatten())
        Rtilde = self.Rtilde(x, h, IDIVIp)
        
        return rho*self.cal**3 * self.KF * Rtilde * xCRSv
    
    def magnusMoment(self, x, v, h, IDIVIp, rho=1.225):
        vDOTx = v.flatten().dot(x.flatten())
        Rtilde = self.Rtilde(x, h, IDIVIp)
        
        return rho*self.cal**4 * self.KT * Rtilde * (xDOTv*x - v)
    
    def pitchingForce(self, x, v, h, rho=1.225):
        V = np.linalg.norm(v)
        hCRSx = np.cross(h.flatten(), x.flatten())
        
        return -rho*self.cal**3 * self.KS * V * hCRSx
    
    def dampingMoment(self, x, v, h, rho=1.225):
        V = np.linalg.norm(v)
        hDOTx = h.flatten().dot(x.flatten())
        
        return -rho*self.cal**4 * self.KH * V * (h - hDOTx*x)
    
    def magnuscrossForce(self, x, v, h, IDIVIp, rho=1.225):
        V = np.linalg.norm(v)
        hDOTx = h.flatten().dot(x.flatten())
        Rtilde = self.Rtilde(x, h, IDIVIp)
        
        return rho*self.cal**4 * self.KXF * Rtilde * V * (h - hDOTx*x)
    
    def magnuscrossMoment(self, x, v, h, IDIVIp, rho=1.225):
        V = np.linalg.norm(v)
        hCRSx = np.cross(h.flatten(), x.flatten())
        Rtilde = self.Rtilde(x, h, IDIVIp)
        
        return -rho*self.cal**5 * self.KXT * Rtilde * hCRSx
    
    
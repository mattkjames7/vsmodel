import numpy as np


def _F(r,phi,Kp):
	alpha = _alpha(phi,Kp)
	Rs = _Rs(phi,Kp)
	F = 0.5 + (1.0/np.pi)*np.arctan((2.0/alpha)*(r-Rs))

	return F
	
def _dGdphi(phi):
	
	Am = [0.53,0.37,0.1]
	Bm = [0.0,0.21,-0.1]
	
	phi0 = np.pi/2.0
	
	G = 0.0
	dphi = phi - phi0
	for m in range(0,3):
		mdphi = m*dphi
		G += (-m*Am[m]*np.sin(mdphi) + m*Bm[m]*np.cos(mdphi))
	
	return G
	
def _G(phi):
	Am = [0.53,0.37,0.1]
	Bm = [0.0,0.21,-0.1]
	
	phi0 = np.pi/2.0
	
	G = 0.0
	dphi = phi - phi0
	for m in range(0,3):
		mdphi = m*dphi
		G += (Am[m]*np.cos(mdphi) + Bm[m]*np.sin(mdphi))
	
	return G	
	
def _Vs(Kp):
	
	return 0.75*Kp**2
	
def _alpha(phi,Kp):
	
	alpha = 0.15 + (2.55-0.27*Kp)*(1.0 + np.cos(phi - (7*np.pi/12.0)))
	
	return alpha
	
def _R0(Kp):
	
	R0 = 4.4 - 0.6*(Kp - 5.0)
	
	return R0
	
def _Rs(phi,Kp):
	
	R0 = _R0(Kp)
	kappa = 0.14
	beta = 0.97
	
	Rs = R0*(((1 + beta)/(1 + beta*np.cos(phi - np.pi)))**kappa)

	return Rs
	
def _dRsdphi(phi,Kp):
	kappa = 0.14
	beta = 0.97	
	R0 = _R0(Kp)
	p = 1 + beta
	q = 1.0 + beta*np.cos(phi - np.pi)
	dqdphi = -beta*np.sin(phi-np.pi)
	
	dRsdphi = -(R0/q)*dqdphi*kappa*((p/q)**kappa)
	return dRsdphi
	
def _dalphadphi(phi,Kp):
	'''
	The derivative of alpha w.r.t. phi
	
	'''
	dadp = -(2.55 - 0.27*Kp)*np.sin(phi - (7*np.pi/12.0))
	return dadp
	
def _dFdr(r,phi,Kp):
	'''
	Derivative of F(r,phi) w.r.t. r
	
	'''
	alpha = _alpha(phi,Kp)
	Rs = _Rs(phi,Kp)
	f = (2.0/alpha)*(r - Rs)
	datanfdf = 1.0/(1.0 + f**2)
	dfdr = 2.0/alpha
	dFdr = (1.0/np.pi)*datanfdf*dfdr
	return dFdr
	
def _Er(r,phi,Kp):
	'''
	r-component of the electric field.
	
	'''
	
	Vs = _Vs(Kp)
	G = _G(phi)
	dFdr = _dFdr(r,phi,Kp)
	
	Er = Vs*G*dFdr
	
	return Er
	
def _dFdphi(r,phi,Kp):
	
	alpha = _alpha(phi,Kp)
	Rs = _Rs(phi,Kp)
	f = (2.0/alpha)*(r - Rs)
	dFdf = 1.0/(1.0 + f**2)

	dalphadphi = _dalphadphi(phi,Kp)
	dgdphi = (-2.0/alpha)*dalphadphi
	
	g = (2/alpha)
	dRsdphi = _dRsdphi(phi,Kp)

	h = (r - Rs)

	dfdphi =  dgdphi*h - g*dRsdphi
	dFdphi = (1.0/np.pi)*dFdf*dfdphi
	
	return dFdphi
	
def _Ep(r,phi,Kp):
	'''
	phi-component of the electric field.
	
	'''
	
	Vs = _Vs(Kp)
	G = _G(phi)
	F = _F(r,phi,Kp)
	dGdphi = _dGdphi(phi)
	dFdphi = _dFdphi(r,phi,Kp)
	
	Ep = (Vs/r)*(dGdphi*F + dFdphi*G)
	
	return Ep


def SAPSPotential(r,phi,Kp):
	'''
	Calculate the SAPS component of the potential as described by
	Goldstein et al 2005 (doi:10.1029/2005JA011135) in cylindrical
	coordinates.
	
	Inputs
	======
	r : float
		Radial distance from the centre of the Earth in Re
	phi : float
		azimuthal angle, 0 at noon. (radians)
	Kp : float
		Kp index
		
	Returns
	=======
	U : float
		Scalar potential in kV.
		
	'''	
	Vs = _Vs(Kp)
	F = _F(r,phi,Kp)
	G = _G(phi)
	
	U = -Vs*F*G
	return U

def SAPSPotentialCart(x,y,Kp):
	'''
	Calculate the SAPS component of the potential as described by
	Goldstein et al 2005 (doi:10.1029/2005JA011135) in cylindrical
	coordinates.
	
	Inputs
	======
	x : float
		Position(s) in x SM direction
	y : float
		Position(s) in y SM direction
	Kp : float
		Kp index
		
	Returns
	=======
	U : float
		Scalar potential in kV.
		
	'''	
	#calculate the cylindrical coordinates
	r = np.sqrt(x**2 + y**2)
	phi = np.arctan2(y,x)
	
	U = SAPSPotential(r,phi,Kp)
	return U

def SAPSModel(r,phi,Kp):
	'''
	Calculate the SAPS component of the potential as described by
	Goldstein et al 2005 (doi:10.1029/2005JA011135) in cylindrical
	coordinates.
	
	Inputs
	======
	r : float
		Radial distance from the centre of the Earth in Re
	phi : float
		azimuthal angle, 0 at noon. (radians)
	Kp : float
		Kp index
		
	Returns
	=======
	Er, Ephi, Ez : float
		Electric field vector in cylindrical coordinates and units 
		of mV/m.
	
	
	'''
	#get each component
	Er = _Er(r,phi,Kp)
	Ep = _Ep(r,phi,Kp)
	Ez = np.zeros(r.shape,dtype='float32')

	#convert from kV/Re to mV/m
	C = 1000000.0/6.378e6
	
	return (Er*C,Ep*C,Ez)
	
def SAPSModelCart(x,y,Kp):
	'''
	Calculate the SAPS component of the potential as described by
	Goldstein et al 2005 (doi:10.1029/2005JA011135) in Cartesian 
	coordinates.
	
	Inputs
	======
	x : float
		Position(s) in x SM direction
	y : float
		Position(s) in y SM direction
	Kp : float
		Kp index
		
	Returns
	=======
	Er, Ephi, Ez : float
		Electric field vector in cylindrical coordinates and units 
		of mV/m.
	
	
	'''
	#calculate the cylindrical coordinates
	r = np.sqrt(x**2 + y**2)
	phi = np.arctan2(y,x)

	#call the cylindrical version of the model
	Er,Ep,Ez = SAPSModel(r,phi,Kp)
	
	#convert to cartesian components
	Ex = Er*np.cos(phi) - Ep*np.sin(phi)
	Ey = Er*np.sin(phi) + Ep*np.cos(phi)

	return Ex,Ey,Ez

import numpy as np


def MCPotential(r,phi,Kp):
	'''
	Maynard-Chen potential in kV.
	
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
	

	#constant related to the convection E field strength
	b = 0.045/((1.0 - 0.159*Kp + 0.0093*Kp**2)**3)
	
	#the shielding parameter
	gamma = 2.0
	
	#the potential
	U = - b*(r**gamma)*np.sin(phi)
	
	
	return U	
	

def MCPotentialCart(x,y,Kp):
	'''
	Maynard-Chen potential in kV.
	
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
	
	#calculate the polar coordinates
	r = np.sqrt(x**2 + y**2)
	phi = np.arctan2(y,x)
	
	#get the potential
	U = MCPotential(r,phi,Kp)
	
	return U


def MCModel(r,phi,Kp):
	'''
	The convection electric field model using the Kp dependence 
	described by Maynard and Chen 1975
	
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
	
	#convert from kV/Re to mV/m
	C = 1000000.0/6.378e6
	
	#constant related to the convection E field strength
	b = 0.045/((1.0 - 0.159*Kp + 0.0093*Kp**2)**3)
	
	#calculate the components in cylindrical coordinates
	gamma = 2.0
	Er = gamma*b*(r**(gamma-1))*np.sin(phi)
	Ep = b*(r**(gamma-1))*np.cos(phi)
	Ez = np.zeros(np.shape(r),dtype='float32')
	
	return (Er*C,Ep*C,Ez*C)
	
def MCModelCart(x,y,Kp):
	'''
	The convection electric field model using the Kp dependence 
	described by Maynard and Chen 1975
	
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
	Ex, Ey, Ez : float
		Electric field vector in SM coordinates and units 
		of mV/m.
	
	
	'''
	#calculate the cylindrical coordinates
	r = np.sqrt(x**2 + y**2)
	phi = np.arctan2(y,x)

	#call the cylindrical version of the model
	Er,Ep,Ez = MCModel(r,phi,Kp)
	
	#convert to cartesian components
	Ex = Er*np.cos(phi) - Ep*np.sin(phi)
	Ey = Er*np.sin(phi) + Ep*np.cos(phi)

	return Ex,Ey,Ez

import numpy as np

def CorotPotential(r):
	'''
	Corotation potential in kV.

	Inputs
	======
	r : float
		Radial distance from the centre of the Earth in Re

	Returns
	=======
	U : float
		Scalar potential in kV.
	
	'''
	#corotation constant in kV Re
	a = 92.4
	
	#the potential
	U = -(a/r)	
	
	return U	


def CorotPotentialCart(x,y):
	'''
	Corotation potential in kV.
	
	Inputs
	======
	x : float
		Position(s) in x SM direction
	y : float
		Position(s) in y SM direction

	Returns
	=======
	U : float
		Scalar potential in kV.
	
	'''
	
	#calculate the polar coordinates
	r = np.sqrt(x**2 + y**2)
	phi = np.arctan2(y,x)
	
	#get the potential
	U = CorotPotential(r)
	
	return U
	

def CorotModel(r):
	'''
	Corotational part of the Vollan-Stern model in cylindrical
	coordinates.
	
	Inputs
	======
	r : float
		Radial distance from the centre of the Earth in Re

	Returns
	=======
	Er, Ephi, Ez : float
		Electric field vector in cylindrical coordinates and units 
		of mV/m.
	
	
	'''
	#corotation constant in kV Re
	a = 92.4
	
	#get the radial component
	Er = -a/(r**2)
	
	#other components are 0
	Ep = np.zeros(np.shape(r),dtype='float32')
	Ez = np.zeros(np.shape(r),dtype='float32')
	
	#convert from kV/Re to mV/m
	C = 1000000.0/6.378e6
		
	return (Er*C,Ep*C,Ez*C)

def CorotModelCart(x,y):
	'''
	Cartesian version of the corotational component of the Volland-Stern
	model.
	
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
		Electric field vector in cartesian coordinates and units 
		of mV/m.
	
	'''
	
	#calculate the cylindrical coordinates
	r = np.sqrt(x**2 + y**2)
	phi = np.arctan2(y,x)
	
	#get the E field
	Er,Ep,Ez = CorotModel(r)
	
	#convert to cartesian components
	Ex = Er*np.cos(phi) - Ep*np.sin(phi)
	Ey = Er*np.sin(phi) + Ep*np.cos(phi)

	return Ex,Ey,Ez

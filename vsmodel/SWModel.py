import numpy as np

def SWPotential(r,phi,Esw=None,Vsw=-400.0,Bz=0.0):
	'''
	Solar wind convection component of the potential as proposed by 
	Goldstein et al 2005.
	
	Inputs
	======
	r : float
		Radial distance from the centre of the Earth in Re
	phi : float
		azimuthal angle, 0 at noon. (radians)
	Esw : float
		Electric field due to the solar wind propagating past the 
		magnetosphere in mV/m.
	Vsw : float
		x-component (positive towards the Sun) of the solar wind 
		velocity in km/s.
	Bz : float
		Z-component of the IMF in nT..
		
	Returns
	=======
	U : float
		Potential in kV
		
	'''
	#check if Esw is provided or whether we use separate Vsw and Bz
	if Esw is None:
		Esw = (np.float64(Bz*1e-9)*np.float64(Vsw*1e6)).clip(min=0.1)
	
	#convert from kV/Re to mV/m
	C = 1000000.0/6.378e6
			
	#convert E to kV/Re
	Esw = Esw/C
	
	#calculate A
	A = 0.12*Esw/6.6	

	#the shielding parameter
	gamma = 2.0
	
	#the potential
	U = - A*(r**gamma)*np.sin(phi)	
	
	return U	

def SWPotentialCart(x,y,Esw=None,Vsw=-400.0,Bz=0.0):
	'''
	Solar wind convection component of the potential as proposed by 
	Goldstein et al 2005.
	
	Inputs
	======
	x : float
		Position(s) in x SM direction
	y : float
		Position(s) in y SM direction
	Esw : float
		Electric field due to the solar wind propagating past the 
		magnetosphere in mV/m.
	Vsw : float
		x-component (positive towards the Sun) of the solar wind 
		velocity in km/s.
	Bz : float
		Z-component of the IMF in nT.
		
	Returns
	=======
	U : float
		Potential in kV
		
	'''
	
	#calculate the cylindrical coordinates
	r = np.sqrt(x**2 + y**2)
	phi = np.arctan2(y,x)

	U = SWPotential(r,phi,Esw,Vsw,Bz)
	
	return U

def SWModel(r,phi,Esw=None,Vsw=-400.0,Bz=0.0):
	'''
	Solar wind convection component as proposed by Goldstein et al 2005.
	
	Inputs
	======
	r : float
		Radial distance from the centre of the Earth in Re
	phi : float
		azimuthal angle, 0 at noon. (radians)
	Esw : float
		Electric field due to the solar wind propagating past the 
		magnetosphere in mV/m.
	Vsw : float
		x-component (positive towards the Sun) of the solar wind 
		velocity in km/s.
	Bz : float
		Z-component of the IMF in nT.
		
	Returns
	=======
	Er, Ephi, Ez : float
		Electric field vector in cylindrical coordinates and units 
		of mV/m.
	
	
	'''
	
	#convert from kV/Re to mV/m
	C = 1000000.0/6.378e6
	
	#check if Esw is provided or whether we use separate Vsw and Bz
	if Esw is None:
		Esw = (np.float64(Bz*1e-9)*np.float64(Vsw*1e6)).clip(min=0.1)
	
	#convert E to kV/Re
	Esw = Esw/C
	
	#calculate A
	A = 0.12*Esw/6.6
	
	#calculate the components in cylindrical coordinates
	gamma = 2.0
	Er = gamma*A*(r**(gamma-1))*np.sin(phi)
	Ep = A*(r**(gamma-1))*np.cos(phi)
	Ez = np.zeros(np.shape(r),dtype='float32')
	
	return (Er*C,Ep*C,Ez*C)
	
def SWModelCart(x,y,Esw=None,Vsw=-400.0,Bz=0.0):
	'''
	Solar wind convection component as proposed by Goldstein et al 2005.
	
	Inputs
	======
	x : float
		Position(s) in x SM direction
	y : float
		Position(s) in y SM direction
	Esw : float
		Electric field due to the solar wind propagating past the 
		magnetosphere in mV/m.
	Vsw : float
		x-component (positive towards the Sun) of the solar wind 
		velocity in km/s.
	Bz : float
		Z-component of the IMF in nT..
		
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
	Er,Ep,Ez = SWModel(r,phi,Esw,Vsw,Bz)
	
	#convert to cartesian components
	Ex = Er*np.cos(phi) - Ep*np.sin(phi)
	Ey = Er*np.sin(phi) + Ep*np.cos(phi)

	return Ex,Ey,Ez

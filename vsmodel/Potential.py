import numpy as np
from .CorotModel import CorotPotential
from .MCModel import MCPotential
from .SAPSModel import SAPSPotential 
from .SWModel import SWPotential

def Potential(r,phi,Kp=1.0,Esw=None,Vsw=None,Bz=None):
	'''
	Get the Volland-Stern electric field model (see 
	doi:10.1029/JA078i001p00171 and doi:10.1029/JA080i004p00595) using 
	the corotational electric field (as defined in 
	doi:10.1002/2017JA024558) combined with either the Maynard and Chen
	1975 (doi:10.1029/JA080i007p01009) or the Goldstein et al 2005
	(doi:10.1029/2005JA011135) convection electric field contribution.
	
	NOTE:
	To use the Maynard and Chen version, set Esw, Bz and Vsw keywords
	equal to None (this is their default). Otherwise, either set Esw or
	both Bz and Vsw to use the Goldstein version.
	
	Inputs
	======
	r : float
		Radial distance from the centre of the Earth in Re.
	phi : float
		azimuthal angle, 0 at noon. (radians)
	Kp : float
		Kp index.
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
		Potential in kV.
	
	
	'''
	
	
	#check which convection model we are going to use
	if (not Esw is None) or (not Vsw is None and not Bz is None):
		# in this case use the Goldstein et al version of the VS model
		
		U0 = SAPSPotential(r,phi,Kp)
		U1 = SWPotential(r,phi,Esw,Vsw,Bz)
		
		Uc = U0 + U1
	else:
		#otherwise use the Maynard and Chen version
		Uc = MCPotential(r,phi,Kp)
		
	#get the corotation E field
	Ur = CorotPotential(r)
	
	#now combine to get the output field
	U = Ur + Uc
	
	return U


def PotentialCart(x,y,Kp=1.0,Esw=None,Vsw=None,Bz=None):
	'''
	Get the Volland-Stern electric field model (see 
	doi:10.1029/JA078i001p00171 and doi:10.1029/JA080i004p00595) using 
	the corotational electric field (as defined in 
	doi:10.1002/2017JA024558) combined with either the Maynard and Chen
	1975 (doi:10.1029/JA080i007p01009) or the Goldstein et al 2005
	(doi:10.1029/2005JA011135) convection electric field contribution.
	
	NOTE:
	To use the Maynard and Chen version, set Esw, Bz and Vsw keywords
	equal to None (this is their default). Otherwise, either set Esw or
	both Bz and Vsw to use the Goldstein version.

	Inputs
	======
	x : float
		Position(s) in x SM direction
	y : float
		Position(s) in y SM direction
	Kp : float
		Kp index.
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
		Potential in kV.
	'''
	#calculate the cylindrical coordinates
	r = np.sqrt(x**2 + y**2)
	phi = np.arctan2(y,x)

	#call the cylindrical version of the model
	U = Potential(r,phi,Kp,Esw,Vsw,Bz)
	
	return U

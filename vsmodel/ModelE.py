import numpy as np
from .CorotModel import CorotModel
from .SWModel import SWModel
from .SAPSModel import SAPSModel
from .MCModel import MCModel

def ModelE(r,phi,Kp=1.0,Esw=None,Vsw=None,Bz=None):
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
	Er, Ep, Ez : float
		Electric field vector in cylindrical coordinates and units 
		of mV/m.
	
	
	'''
	
	
	#check which convection model we are going to use
	if (not Esw is None) or (not Vsw is None and not Bz is None):
		# in this case use the Goldstein et al version of the VS model
		
		Er0,Ep0,Ez0 = SAPSModel(r,phi,Kp)
		Er1,Ep1,Ez1 = SWModel(r,phi,Esw,Vsw,Bz)
		
		Erc = Er0 + Er1
		Epc = Ep0 + Ep1
		Ezc = Ez0 + Ez1
	else:
		#otherwise use the Maynard and Chen version
		Erc,Epc,Ezc = MCModel(r,phi,Kp)
		
	#get the corotation E field
	Err,Epr,Ezr = CorotModel(r)
	
	#now combine to get the output field
	Er = Erc + Err
	Ep = Epc + Epr
	Ez = Ezc + Ezr
	
	return Er,Ep,Ez


	
def ModelECart(x,y,Kp=1.0,Esw=None,Vsw=None,Bz=None):
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
	Ex, Ey, Ez : float
		Electric field vector in SM coordinates and units 
		of mV/m.	
	'''
	#calculate the cylindrical coordinates
	r = np.sqrt(x**2 + y**2)
	phi = np.arctan2(y,x)

	#call the cylindrical version of the model
	Er,Ep,Ez = ModelE(r,phi,Kp,Esw,Vsw,Bz)
	
	#convert to cartesian components
	Ex = Er*np.cos(phi) - Ep*np.sin(phi)
	Ey = Er*np.sin(phi) + Ep*np.cos(phi)

	return Ex,Ey,Ez

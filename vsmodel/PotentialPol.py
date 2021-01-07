import numpy as np


def PotentialPol(r,phi,Kp):
	'''
	An attempt to construct the Volland-Stern model, see
	doi:10.1029/JA078i001p00171 and doi:10.1029/JA080i004p00595
	
	This model is created using the following description of the model:
	doi: 10.1002/2017JA024558
	
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
	
	#corotation constant in kV Re
	a = 92.4
	
	#constant related to the convection E field strength
	b = 0.045/((1.0 - 0.159*Kp + 0.0093*Kp**2)**3)
	
	#the shielding parameter
	gamma = 2.0
	
	#the potential
	U = -(a/r) - b*(r**gamma)*np.sin(phi)
	
	
	return U

import numpy as np
from .PotentialPol import PotentialPol

def PotentialCart(x,y,Kp):
	'''
	Cartesian version of an attempt to construct the Volland-Stern model, 
	see	doi:10.1029/JA078i001p00171 and doi:10.1029/JA080i004p00595
	
	This model is created using the following description of the model:
	doi: 10.1002/2017JA024558
	
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
	U = PotentialPol(r,phi,Kp)
	
	return U

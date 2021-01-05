import numpy as np


def ModelPol(r,phi,Kp):
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
	Er, Etheta, Ephi : float
		Electric field vector in spherical polar coordinates and units 
		of mV/m.
	
	'''
	
	#corotation constant in kV Re
	a = 92.4
	
	#constant related to the convection E field strength
	b = 0.045/((1.0 - 0.159*Kp + 0.0093*Kp**2)**3)
	
	#the shielding parameter
	gamma = 2.0
	
	#dU/dr
	dUdr = a*(r**-2) - gamma*b*(r**(gamma - 1.0))*np.sin(phi)
	
	#dU/dphi
	dUdphi = -b*(r**gamma)*np.cos(phi)
	
	#dU/dtheta
	dUdtheta = np.zeros(dUdr.shape,dtype='float32')
	
	#Er
	Er = -dUdr
	
	#Etheta = -(1/r)*(dU/dtheta)
	Etheta = -dUdtheta
	
	#Ephi = -(1/(r*sin(theta)))*(dU/dphi)
	Ephi = -(1.0/r)*dUdphi
	
	#convert from kV/Re to mV/m
	C = 1000000.0/6.378e6
	
	
	return (Er*C,Etheta*C,Ephi*C)

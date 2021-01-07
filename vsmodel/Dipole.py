import numpy as np


def GetDipole(Beq=-31200.0):
	'''
	Create a function which provides magnetic field vector at a given positional vector.
	
	Args:
		Beq: Magnetic field strength at the equator in nT
		
	Returns:
		Function which takes positional vector p and returns a field vector.
	'''
	def DipoleField(x,y,z):
		
		#get the shape of x
		sh = np.shape(x)
		_x = np.array(x).flatten()
		_y = np.array(y).flatten()
		_z = np.array(z).flatten()
		
		#p must be a vector, or array of vectors which are in Rp
		mu = np.array([0.0,0.0,Beq]) #nT R_p^3
		M = np.linalg.norm(mu)
		mhat = np.array([mu/M])
		
		#get R and its unit vector
		R = np.sqrt(_x**2 + _y**2 + _z**2)
		if np.size(-x) == 1:
			rhat = np.array([[_x/R],[_y/R],[_z/R]]).T
		else:		
			rhat = np.array([_x/R,_y/R,_z/R]).T
		
		#now calcualte the field components
		mdotr = np.sum(mhat*rhat,axis=1)
		B = M*((3.0*mdotr*rhat.T - mhat.T)/R**3.0)
		
		Bx,By,Bz = B[0].reshape(sh),B[1].reshape(sh),B[2].reshape(sh)
		return Bx,By,Bz
		
	
	return DipoleField

import numpy as np

def VExB(Ex,Ey,Ez,Bx,By,Bz):
	'''
	Calculate the ExB velocity based on the electric and magnetic fields
	
	Inputs
	======
	Ex : float
		x-component of the electric field, V/m
	Ey : float
		y-component of the electric field, V/m
	Ez : float
		z-component of the electric field, V/m
	Bx : float
		x-component of the magnetic field, T
	By : float
		y-component of the magnetic field, T
	Bz : float
		z-component of the magnetic field, T
	
	Returns
	=======
	Vx : float
		x-component of the drift velocity
	Vy : float
		y-component of the drift velocity
	Vz : float
		z-component of the drift velocity
	
	'''
	
	#square of the field magnitude
	B2 = Bx**2 + By**2 + Bz**2
	
	#V = (E x B)/B^2
	Vx = Ey*Bz - Ez*By
	Vy = Ez*Bx - Ex*Bz
	Vz = Ex*By - Ey*Bx
	
	Vx /= B2
	Vy /= B2
	Vz /= B2
	
	return Vx,Vy,Vz
	

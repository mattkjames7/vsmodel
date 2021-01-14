import numpy as np
from .ModelE import ModelECart
from .Dipole import GetDipole
from .VExB import VExB

def VExBModel(x,y,z,Kp,Esw=None,Vsw=None,Bz=None):
	'''
	Estimate the ExB drift using a dipole magnetic field with the
	Volland-Stern E field model.
	
	'''

	#get the dipole field function
	dip = GetDipole()
	
	#work out B
	Bx,By,Bz = dip(x,y,z)
	
	#convert to T
	Bx *= 1e-9
	By *= 1e-9
	Bz *= 1e-9

	#get the E field
	Ex,Ey,Ez = ModelECart(x,y,Kp,Esw,Vsw,Bz)
	
	#convert from mV/m to V/m
	Ex *= 1e-3
	Ey *= 1e-3
	Ez *= 1e-3
	
	#calculate the ExB drift velocity vector
	Vx,Vy,Vz = VExB(Ex,Ey,Ez,Bx,By,Bz)
	
	return Vx,Vy,Vz

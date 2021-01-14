import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from .PlotPlanet import PlotPlanetXY
from .ModelE import ModelECart
from .Potential import PotentialCart
from .VExBModel import VExBModel
from .Dipole import GetDipole
from .CorotModel import CorotModelCart,CorotPotentialCart
from .MCModel import MCModelCart,MCPotentialCart
from .SAPSModel import SAPSModelCart,SAPSPotentialCart
from .SWModel import SWModelCart,SWPotentialCart
	
zlabs = {	'U':		'Potential (kV)',
			'E':		'Electric Field, $|\mathbf{E}|$, (mV m$^{-1}$)',
			'Ex':		'Electric Field, $E_x$, (mV m$^{-1}$)',
			'Ey':		'Electric Field, $E_y$, (mV m$^{-1}$)',
			'V':		'Velocity, $|\mathbf{V}|$, (m s$^{-1}$)',
			'Vx':		'Velocity, $V_x$, (m s$^{-1}$)',
			'Vy':		'Velocity, $V_y$, (m s$^{-1}$)',
			'B':		'Magnetic Field, $|\mathbf{B}|$, (nT)',
			'Bx':		'Magnetic Field, $B_x$, (nT)',
			'By':		'Magnetic Field, $B_y$, (nT)',
			'Bz':		'Magnetic Field, $B_z$, (nT)',
			'U_SAPS': 	'SAPS Potential (kV)',
			'U_Corot': 	'Corotation Potential (kV)',
			'U_MC':		'Maynard and Chen 1975 Potential (kV)',
			'U_SW':		'SW Convection Potential (kV)',
			'E_SAPS':	'SAPS Electric Field, $|\mathbf{E}|$, (mV m$^{-1}$)',
			'Ex_SAPS':	'SAPS Electric Field, $E_x$, (mV m$^{-1}$)',
			'Ey_SAPS':	'SAPS Electric Field, $E_y$, (mV m$^{-1}$)',
			'E_Corot':	'SAPS Electric Field, $|\mathbf{E}|$, (mV m$^{-1}$)',
			'Ex_Corot':	'SAPS Electric Field, $E_x$, (mV m$^{-1}$)',
			'Ey_Corot':	'SAPS Electric Field, $E_y$, (mV m$^{-1}$)',
			'E_MC':		'Maynard and Chen 1975 Electric Field, $|\mathbf{E}|$, (mV m$^{-1}$)',
			'Ex_MC':	'Maynard and Chen 1975 Electric Field, $E_x$, (mV m$^{-1}$)',
			'Ey_MC':	'Maynard and Chen 1975 Electric Field, $E_y$, (mV m$^{-1}$)',
			'E_SW':		'SW Convection Electric Field, $|\mathbf{E}|$, (mV m$^{-1}$)',
			'Ex_SW':	'SW Convection Electric Field, $E_x$, (mV m$^{-1}$)',
			'Ey_SW':	'SW Convection Electric Field, $E_y$, (mV m$^{-1}$)',
			}
			
def PlotModelEq(Model='V',Rmax=10.0,dR=0.1,Kp=1.0,Esw=None,Vsw=None,Bz=None,fig=None,maps=[1,1,0,0],
		ShowContour=True,zlog=True,cmap='gnuplot',scale=None,Verbose=False,fmt='%1.0f',nlvl=10):
	'''
	Plots various models in the equatorial plane.

	NOTE:
	To use the Maynard and Chen version, set Esw, Bz and Vsw keywords
	equal to None (this is their default). Otherwise, either set Esw or
	both Bz and Vsw to use the Goldstein version.
	
	Inputs
	======
	Model : str
		'U' - potential
		'Ex'|'Ey'|'E' - Electric field
		'Vx'|'Vy'|'V' - ExB drift
		'Bx'|'By'|'Bz'|'B' - Magnetic dipole field
	Rmax : float
		Limits of the plot
	dR : float
		Size of bins on plot
	Kp : float
		Kp index
	Esw : float
		Electric field due to the solar wind propagating past the 
		magnetosphere in mV/m.
	Vsw : float
		x-component (positive towards the Sun) of the solar wind 
		velocity in km/s.
	Bz : float
		Z-component of the IMF in nT.
	fig : object
		This should be either None (to create a new figure), a matplotlib.pyplot
		object (to add a new set of axes to an existing plot) or an 
		instance of matplotlib.pyplot.Axes in order to add to an existing
		set of axes.
	maps : list/tuple
		This will determine the position of the subplot on the figure.
	zlog : bool
		If True - the color scale will be logarithmic
	cmap : str
		String to say which colormap to use.
	scale : None or list/tuple
		If None - scale will be generated automatically, otherwise set 
		scale = [minimum,maximum]
	Verbose : bool
		If True, model calculation progress will be displayed	
	
	
	
	'''
	
	#create a grid of points to calculate the model at
	nR = 2*np.int32(Rmax/dR)
	xe = np.linspace(-Rmax,Rmax,nR+1)
	ye = xe
	xc = xe[:-1] + 0.5*dR
	yc = ye[:-1] + 0.5*dR
	xe,ye = np.meshgrid(xe,ye)
	xc,yc = np.meshgrid(xc,yc)
	zc = np.zeros(xc.shape)
	R = np.sqrt(xc**2 + yc**2 + zc**2)
	
	#calculate the model
	if Model == 'U':
		data = PotentialCart(xc,yc,Kp,Esw,Vsw,Bz)
	elif Model in ['E','Ex','Ey']:
		data = ModelECart(xc,yc,Kp,Esw,Vsw,Bz)
	elif Model in ['B','Bx','By','Bz']:
		dip = GetDipole()
		data = dip(xc,yc,zc)
	elif Model in ['V','Vx','Vy']:
		data = VExBModel(xc,yc,zc,Kp,Esw,Vsw,Bz)
	elif 'E' in Model and 'SAPS' in Model:
		data = SAPSModelCart(xc,yc,Kp)
	elif 'E' in Model and 'Corot' in Model:
		data = CorotModelCart(xc,yc)
	elif 'E' in Model and 'MC' in Model:
		data = MCModelCart(xc,yc,Kp)
	elif 'E' in Model and 'SW' in Model:
		data = SWModelCart(xc,yc,Esw,Vsw,Bz)
	elif 'U' in Model and 'SAPS' in Model:
		data = SAPSPotentialCart(xc,yc,Kp)
	elif 'U' in Model and 'Corot' in Model:
		data = CorotPotentialCart(xc,yc)
	elif 'U' in Model and 'MC' in Model:
		data = MCPotentialCart(xc,yc,Kp)
	elif 'U' in Model and 'SW' in Model:
		data = SWPotentialCart(xc,yc,Esw,Vsw,Bz)
	
	if 'x' in Model:
		data = data[0]
	elif 'y' in Model:
		data = data[1]
	elif 'z' in Model:
		data = data[2]
	
	if not ('U' in Model or 'x' in Model or 'y' in Model or 'z' in Model):
		data = np.sqrt(data[0]**2 + data[1]**2 + data[2]**2)
	
	#remove stuff from within the planet
	bad = np.where(R < 1.0)
	data[bad] = np.nan
	
	#get the scale limits
	if scale is None:
		if zlog:
			scale = [np.nanmin(data[data > 0]),np.nanmax(data)]
		else:
			scale = [np.nanmin(data),np.nanmax(data)]
			
	zlabel = zlabs[Model]
	if zlog:
		norm = colors.LogNorm(vmin=scale[0],vmax=scale[1])
		lvl = 10**np.linspace(np.log10(scale[0]),np.log10(scale[1]),nlvl)
	else:
		norm = colors.Normalize(vmin=scale[0],vmax=scale[1])	
		lvl = np.linspace(scale[0],scale[1],nlvl)
	

	

	#create the plot window
	if fig is None:
		fig = plt
		fig.figure()
	if hasattr(fig,'Axes'):
		ax = fig.subplot2grid((maps[1],maps[0]),(maps[3],maps[2]))
		ax.set_xlim([Rmax,-Rmax])
		ax.set_ylim([-Rmax,Rmax])
		ax.set_aspect(1.0)
	else:
		ax = fig		
		
	#plot the mesh
	sm = ax.pcolormesh(ye,xe,data,cmap=cmap,norm=norm,)

	if ShowContour:
		cs = ax.contour(yc,xc,data,lvl,cmap='Greys',norm=norm)
		ax.clabel(cs, inline=1, fontsize=10,fmt=fmt)
		
	#plot the planet
	PlotPlanetXY(ax)


	#sort the axis labels
	ax.set_ylabel('$x_{SM}$ ($R_E$)')
	ax.set_xlabel('$y_{SM}$ ($R_E$)')	
	
	

	#colorbar
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="5%", pad=0.05)
	cbar = plt.colorbar(sm,cax=cax) 
	cbar.set_label(zlabel)

	return ax

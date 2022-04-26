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
from .VExB import VExB
	
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
			'V_SAPS':	'SAPS Velocity, $|\mathbf{V}|$, (m s$^{-1}$)',
			'Vx_SAPS':	'SAPS Velocity, $V_x$, (m s$^{-1}$)',
			'Vy_SAPS':	'SAPS Velocity, $V_y$, (m s$^{-1}$)',
			'V_Corot':	'Corotation Velocity, $|\mathbf{V}|$, (m s$^{-1}$)',
			'Vx_Corot':	'Corotation Velocity, $V_x$, (m s$^{-1}$)',
			'Vy_Corot':	'Corotation Velocity, $V_y$, (m s$^{-1}$)',
			'V_MC':		'Maynard and Chen 1975 Velocity, $|\mathbf{V}|$, (m s$^{-1}$)',
			'Vx_MC':	'Maynard and Chen 1975 Velocity, $V_x$, (m s$^{-1}$)',
			'Vy_MC':	'Maynard and Chen 1975 Velocity, $V_y$, (m s$^{-1}$)',
			'V_SW':		'SW Convection Velocity, $|\mathbf{V}|$, (m s$^{-1}$)',
			'Vx_SW':	'SW Convection Velocity, $V_x$, (m s$^{-1}$)',
			'Vy_SW':	'SW Convection Velocity, $V_y$, (m s$^{-1}$)',
			}

#A function wrapper for the potential
def _PotentialCart(x,y,Kp,Esw,Vsw,Bz,Emin,Escale):
	return PotentialCart(x,y,Kp,Esw,Vsw,Bz,Emin,Escale)
			
#List the models 0: Full; 1: Corot; 2: MC; 3: SW; 4: SAPS
UModels = [PotentialCart,CorotPotentialCart,MCPotentialCart,SWPotentialCart,SAPSPotentialCart]
EModels = [ModelECart,CorotModelCart,MCModelCart,SWModelCart,SAPSModelCart]
BModel = GetDipole()

#map the model string to the index
modelmap = {'Full' : 0,
			'Corot' : 1,
			'MC' : 2,
			'SW' : 3,
			'SAPS' : 4}


def PlotModelEq(Model='V',Rmax=10.0,dR=0.1,Kp=1.0,Esw=None,Vsw=None,
		Bz=None,Emin=0.25,Escale=0.2,fig=None,maps=[1,1,0,0],
		ShowContour=True,zlog=True,cmap='gnuplot',scale=None,
		Verbose=False,fmt='%1.0f',nlvl=10):
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
		The above strings can also be combined with one of the following
		substrings in order to look at a specific model component: 
		'_Corot'|'_MC'|'_SW'|'_SAPS' - corresponding to the corotational,
		Maynard and Chen 1975, Solar wind driven convection and SAPS
		components. e.g. 'Ex_SAPS' will present the x-component of the 
		electric field due to SAPS.
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
	
	
	#work out which model to use
	FullModel = not '_' in Model
	if FullModel:
		Mtype = Model
		Mcomp = 'Full'
	else:
		Mtype,Mcomp = Model.split('_')
	Mindex = modelmap[Mcomp]
	if Mcomp == 'Full':
		args = (xc,yc,Kp,Esw,Vsw,Bz,Emin,Escale)
	elif Mcomp in ['SAPS','MC']:
		args = (xc,yc,Kp)
	elif Mcomp == 'SW':
		args = (xc,yc,Esw,Vsw,Bz,Emin,Escale)
	else:
		args = (xc,yc)
	
	useEModel = 'E' in Model or 'V' in Model
	useBModel = 'B' in Model or 'V' in Model
	if useEModel:
		EModel = EModels[Mindex]
		E = EModel(*args)
	if useBModel:
		B = BModel(xc,yc,zc)
	
	if 'U' in Model:
		UModel = UModels[Mindex]
		data = UModel(*args)
	else:
		if 'B' in Mtype:
			data = B
		elif 'E' in Mtype:
			data = E
		else:
			data = VExB(E[0]*1e-3,E[1]*1e-3,E[2]*1e-3,B[0]*1e-9,B[1]*1e-9,B[2]*1e-9)
	
		if 'x' in Model:
			data = data[0]
		elif 'y' in Model:
			data = data[1]
		elif 'z' in Model:
			data = data[2]
		else:
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

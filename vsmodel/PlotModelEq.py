import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from .PlotPlanet import PlotPlanetXY
from .ModelCart import ModelCart
from .PotentialCart import PotentialCart
from .VExBModel import VExBModel
from .Dipole import GetDipole

zlabs = {	'U':	'Potential (V)',
			'E':	'Electric Field, $|\mathbf{E}|$, (mV m$^{-1}$)',
			'Ex':	'Electric Field, $E_x$, (mV m$^{-1}$)',
			'Ey':	'Electric Field, $E_y$, (mV m$^{-1}$)',
			'V':	'Velocity, $|\mathbf{V}|$, (m s$^{-1}$)',
			'Vx':	'Velocity, $V_x$, (m s$^{-1}$)',
			'Vy':	'Velocity, $V_y$, (m s$^{-1}$)',
			'B':	'Magnetic Field, $|\mathbf{B}|$, (nT)',
			'Bx':	'Magnetic Field, $B_x$, (nT)',
			'By':	'Magnetic Field, $B_y$, (nT)',
			'Bz':	'Magnetic Field, $B_z$, (nT)',
			}
			
def PlotModelEq(Model='V',Rmax=10.0,dR=0.1,Kp=1.0,fig=None,
		maps=[1,1,0,0],zlog=True,cmap='gnuplot',scale=None,Verbose=False):
	'''
	Plots various models in the equatorial plane.
	
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
		data = PotentialCart(xc,yc,Kp)
	elif 'E' in Model:
		data = ModelCart(xc,yc,Kp)
	elif 'B' in Model:
		dip = GetDipole()
		data = dip(xc,yc,zc)
	elif 'V' in Model:
		data = VExBModel(xc,yc,zc,Kp)
	
	if 'x' in Model:
		data = data[0]
	elif 'y' in Model:
		data = data[1]
	elif 'z' in Model:
		data = data[2]
	
	if Model in ['E','B','V']:
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
	else:
		norm = colors.Normalize(vmin=scale[0],vmax=scale[1])	
	
	

	

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
